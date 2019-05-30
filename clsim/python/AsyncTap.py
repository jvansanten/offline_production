#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: AsyncTap.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file AsyncTap.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from __future__ import print_function

try:
   import zmq
   has_zmq = True
except ImportError:
   has_zmq = False
has_zmq = False

import multiprocessing
import tempfile
import signal
import sys
import time

from icecube import icetray, dataclasses, dataio
from I3Tray import *


# run a writer tray, getting frames from a queue (this runs as a subprocess)
def RunAsyncTray(source,sink,prefix,streams,segment,segmentArgs,identifier,debug):
    
    # just pushes frame on a queue
    class AsyncReceiver(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddParameter("Debug", "Output some status information for debugging", False)
            self.AddParameter("Source", "", None)
            self.AddParameter("Prefix", "", None)
            self.AddParameter("Identifier", "", None)
            self.AddOutBox("OutBox")
            self.nframes = 0
        def Configure(self):
            self.source = self.GetParameter("Source")
            
            self.prefix = self.GetParameter("Prefix")
            self.identifier = self.GetParameter("Identifier")
            if not isinstance(self.identifier, int):
                raise TypeError("Identifier must be an integer")

        def Process(self):
            if self.prefix is not None:
                for frame in dataio.I3File(self.prefix):
                    self.PushFrame(frame)
                self.prefix = None
            
            # request work
            self.source.send_pyobj(self.identifier)
            frame = self.source.recv_pyobj()
            if frame is None:
                icetray.logging.log_debug("({}) got sentinel, requesting suspension".format(self.identifier), unit="clsim.AsyncTap")
                self.RequestSuspension()
                return
            else:
                self.nframes += 1
                self.PushFrame(frame)

        def Finish(self):
            icetray.logging.log_debug("received {} frames".format(self.nframes), unit="clsim.AsyncTap")
    
    class AsyncSender(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddParameter("Sink", "", None)
            self.AddParameter("Streams", "", [icetray.I3Frame.DAQ])
            self.AddOutBox("OutBox")
            self.nframes = 0
        def Configure(self):
            self.sink = self.GetParameter("Sink")
            self.streams = self.GetParameter("Streams")

        def Process(self):
            
            frame = self.PopFrame()
            if frame is not None and frame.Stop in self.streams:
                self.sink.send_pyobj(frame)
                self.nframes += 1
                self.PushFrame(frame)
        
        def Finish(self):
            icetray.logging.log_debug("sent {} frames".format(self.nframes), unit="clsim.AsyncTap")
    
    class RerouteSigIntToDefault(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddOutBox("OutBox")
        def Configure(self):
            # Reroute ctrl+C (sigint) to the default handler.
            # (This needs to be done in a module, as IceTray re-routes
            # the signal in tray.Execute())
            signal.signal(signal.SIGINT, signal.SIG_IGN)
        def DAQ(self, frame):
            self.PushFrame(frame)

    zmq_context = zmq.Context()
    icetray.logging.log_debug("setting up zmq to receive on {}".format(source), unit="clsim.AsyncTap")
    icetray.logging.log_debug("setting up zmq to push on {}".format(sink), unit="clsim.AsyncTap")

    source_socket = zmq_context.socket(zmq.REQ)
    sink_socket = zmq_context.socket(zmq.PUSH)
    sink_socket.setsockopt(zmq.SNDHWM, 1)
    source_socket.connect(source)
    sink_socket.connect(sink)

    try:
        writeTray = I3Tray()
        writeTray.AddModule(AsyncReceiver, "theAsyncReceiver",
            Source=source_socket,
            Prefix=prefix,Debug=debug, Identifier=identifier)
        writeTray.AddModule(RerouteSigIntToDefault, "rerouteSigIntToDefault")
    
        if hasattr(segment, '__i3traysegment__'):
            writeTray.AddSegment(segment,"theSegment", **segmentArgs)
        else:
            writeTray.AddModule(segment,"theModule", **segmentArgs)
    
        writeTray.AddModule(AsyncSender, "theAsyncReturner",
            Sink=sink_socket, Streams=streams)

        icetray.logging.log_info("worker starting..", unit="clsim")
        writeTray.Execute()
        writeTray.Finish()
    finally:
        sink_socket.send_pyobj(identifier)
    icetray.logging.log_debug("({}) finished".format(identifier), unit="clsim.AsyncTap")

class AsyncTap(icetray.I3ConditionalModule):
    """
    Starts N copies of a module or a segment on its own tray in its own process
    and load balances frames from the current tray to the children. Since the
    load balancing destroys the order of the frames anyway, no effort is made
    to preserve the original order when they return to parent tray.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("Segment", "The tray segment to run asynchronously", None)
        self.AddParameter("Args", "A dictionary with keyword arguments for the asynchronous segment.", dict())
        self.AddParameter("NWorkers", "Number of subprocesses to start", 1)
        self.AddParameter("Streams", "Streams to distribute to workers", [icetray.I3Frame.DAQ])
        self.AddParameter("Prefix", "File with frames to prepend to each worker's stream", None)
        self.AddOutBox("OutBox")
        self.nframes = 0
        self.suspensionRequested=False

    def CheckOnChildProcess(self):
        if self.suspensionRequested: return False
        
        # check to see if child process is still alive
        if not all((proc.is_alive for proc in self.processes)):
            icetray.logging.log_error("child tray died unexpectedly. Terminating master tray.", unit="clsim.AsyncTap")
            self.RequestSuspension()
            self.suspensionRequested=True
            return False
        else:
            return True

    def Configure(self):
        self.segment = self.GetParameter("Segment")
        self.args = self.GetParameter("Args")
        nworkers = self.GetParameter("NWorkers")
        assert nworkers > 0

        icetray.logging.log_debug("starting child process..", unit="clsim.AsyncTap")
        
        self.debug = True
        prefix = self.GetParameter("Prefix")
        self.streams = self.GetParameter("Streams")
        source, sink = ('ipc://'+tempfile.mktemp(prefix='clsim-ventilator-') for _ in range(2))
        self.processes = [multiprocessing.Process(target=RunAsyncTray, args=(source,sink,prefix,self.streams,self.segment,self.args,i,True,)) for i in range(nworkers)]
        self.child_ids = set(range(nworkers))

        icetray.logging.log_debug("binding to zmq REP socket {}".format(source), unit="clsim.AsyncTap")
        icetray.logging.log_debug("binding to zmq PULL socket {}".format(sink), unit="clsim.AsyncTap")

        self.zmq_context = zmq.Context()
        self.source = self.zmq_context.socket(zmq.REP)
        self.sink = self.zmq_context.socket(zmq.PULL)
        self.sink.setsockopt(zmq.RCVHWM, 1)
        self.source.bind(source)
        self.sink.bind(sink)
        self.poller = zmq.Poller()
        self.poller.register(self.source, zmq.POLLIN)
        self.poller.register(self.sink, zmq.POLLIN)
        
        icetray.logging.log_debug("zmq socket is bound.", unit="clsim.AsyncTap")

        for proc in self.processes:
            proc.start()

    def shutdown_workers(self):
        while self.child_ids:
            sockets = dict(self.poller.poll())
            if self.sink in sockets:
                self.consume_sink()
            if self.source in sockets:
                msg = self.source.recv_pyobj()
                self.child_ids.add(msg)
                icetray.logging.log_info("shutting down worker {}".format(msg), unit="clsim.AsyncTap")
                self.source.send_pyobj(None)

    def consume_sink(self):
        msg = self.sink.recv_pyobj()
        if isinstance(msg, int):
            icetray.logging.log_info("worker {} finished".format(msg), unit="clsim.AsyncTap")
            self.child_ids.remove(msg)
        else:
            self.PushFrame(msg)

    def Process(self):
        
        frame = self.PopFrame()
        if not frame.Stop in self.streams:
            self.PushFrame(frame)
            return

        while self.child_ids:
            sockets = dict(self.poller.poll())
            if self.sink in sockets:
                self.consume_sink()
            if self.source in sockets:
                msg = self.source.recv_pyobj()
                self.child_ids.add(msg)
                self.source.send_pyobj(frame)
                self.nframes += 1
                break
        else:
            raise RuntimeError("No workers left alive")

    def Finish(self):
        self.shutdown_workers()
        icetray.logging.log_debug("sent {} frames".format(self.nframes), unit="clsim.AsyncTap")

        icetray.logging.log_debug("async module finished. waiting for child tray..", unit="clsim.AsyncTap")
        for proc in self.processes:
            proc.join()
        icetray.logging.log_debug("child tray finished.", unit="clsim.AsyncTap")

