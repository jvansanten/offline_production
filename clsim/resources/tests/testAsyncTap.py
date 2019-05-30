
import sys

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray
import os, time
try:
    from icecube.clsim.AsyncTap import AsyncTap
except ImportError:
    sys.exit(0)

# icetray.logging.set_level_for_unit('clsim.AsyncTap', 'DEBUG')
gcd = os.path.expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz")

@icetray.traysegment
def Tappy(tray, name, raise_error=False):
    class BufferMatic(icetray.I3Module):
        """
        Buffer frames, breaking 1-to-1 reply pattern
        """
        def __init__(self, ctx):
            super(BufferMatic,self).__init__(ctx)
        def Configure(self):
            if raise_error:
                raise ValueError
            pass
            self.frames = []
        def DAQ(self, frame):
            self.frames.append(frame)
        def Finish(self):
            for frame in self.frames:
                self.PushFrame(frame)
    tray.Add(BufferMatic)
    def printy(frame):
        # print("{} in PID {}".format(frame['num'].value, os.getpid()))
        frame['pid'] = icetray.I3Int(os.getpid())
    tray.Add(printy, Streams=[icetray.I3Frame.DAQ])
workers = 4
nframes = 100

def test_with_errors():
    tray = I3Tray()

    tray.Add("I3InfiniteSource", Prefix=gcd)

    def number(frame):
        frame['num'] = icetray.I3Int(number.count)
        number.count += 1
    number.count = 0
    tray.Add(number, Streams=[icetray.I3Frame.DAQ])


    tray.Add(AsyncTap, Segment=Tappy, Prefix=gcd, NWorkers=workers, Args={'raise_error': True}, FrameKey='num')

    try:
        tray.Execute(nframes+3)
        assert False, "tray.Execute() should have raised"
    except RuntimeError:
        pass
    tray.Finish()

def test_num_frames():
    tray = I3Tray()

    tray.Add("I3InfiniteSource", Prefix=gcd)

    def number(frame):
        frame['num'] = icetray.I3Int(number.count)
        number.count += 1
    number.count = 0
    tray.Add(number, Streams=[icetray.I3Frame.DAQ])


    tray.Add(AsyncTap, Segment=Tappy, Prefix=gcd, NWorkers=workers, FrameKey='num')

    def printy(frame):
        printy.workers[frame['pid'].value] += 1
    from collections import defaultdict
    printy.workers = defaultdict(int)
    tray.Add(printy, streams=[icetray.I3Frame.DAQ])

    tray.Execute(nframes+3)
    tray.Finish()
    print(printy.workers)
    assert len(printy.workers) == workers
    assert sum(printy.workers.values()) == nframes

test_num_frames()
test_with_errors()

