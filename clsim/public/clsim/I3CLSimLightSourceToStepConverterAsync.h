/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3CLSimLightSourceToStepConverterAsync.h 179360 2020-03-10 16:07:35Z eganster $
 *
 * @file I3CLSimLightSourceToStepConverterAsync.h
 * @version $Revision: 179360 $
 * @date $Date: 2020-03-10 10:07:35 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMLIGHTSOURCETOSTEPCONVERTERASYNC_H_INCLUDED
#define I3CLSIMLIGHTSOURCETOSTEPCONVERTERASYNC_H_INCLUDED

#include "clsim/I3CLSimLightSourceToStepConverter.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"

I3_FORWARD_DECLARATION(I3CLSimLightSourcePropagator);

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "clsim/I3CLSimQueue.h"

#include <map>
#include <string>

/**
 * @brief A particle-to-step converter using a mix of MC propagators and
 * parameterizations running in a thread.
 */
class I3CLSimLightSourceToStepConverterAsync : public I3CLSimLightSourceToStepConverter
{
public:
    typedef std::tuple<I3CLSimStepSeriesConstPtr, std::vector<uint32_t>, std::map<uint32_t, I3MCTreePtr>, bool> FromGeant4Pair_t;
    
    static const uint32_t default_maxQueueItems;
    
    I3CLSimLightSourceToStepConverterAsync(uint32_t maxQueueItems=default_maxQueueItems);
    virtual ~I3CLSimLightSourceToStepConverterAsync();

    // inherited:
    
    /**
     * Sets the granularity of the bunch size for the
     * return vectors. The number of entries in a vector
     * will always be a multiple of this number.
     * This will stop the worker thread; restart it with
     * a call to Initialize() before enqueuing steps.
     */
    virtual void SetBunchSizeGranularity(uint64_t num);

    /**
     * Sets the maximum bunch size.
     * This will stop the worker thread; restart it with
     * a call to Initialize() before enqueuing steps.
     */
    virtual void SetMaxBunchSize(uint64_t num);

    /**
     * Sets the wavelength bias. Set this to a constant value
     * of 1 if you do not need biased photon generation.
     * The Cherenkov spectrum will be multiplied by this
     * value at each wavelength.
     * This will influence the number of photons produced.
     * This will stop the worker thread; restart it with
     * a call to Initialize() before enqueuing steps.
     */
    virtual void SetWlenBias(I3CLSimFunctionConstPtr wlenBias);

    /**
     * Sets the medium properties.
     * This will stop the worker thread; restart it with
     * a call to Initialize() before enqueuing steps.
     */
    virtual void SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties);

    /**
     * Sets the propagators to use.
     * This will stop the worker thread; restart it with
     * a call to Initialize() before enqueuing steps.
     */
    virtual void SetPropagators(const std::vector<I3CLSimLightSourcePropagatorPtr> &);

    /**
     * Initialize the simulation and start the worker thread
     */
    virtual void Initialize();

    /**
     * Returns true if the worker thread has been started.
     * Never throws.
     */
    virtual bool IsInitialized() const;
    
    /**
     * Adds a new I3Particle to the queue for use as a primary in tracking.
     * The resulting I3CLSimSteps can be retrieved from the
     * I3CLSimLightSourceToStepConverter after some processing time.
     *
     * Enqueuing a particle after calling EnqueueBarrier 
     * will throw if not all steps have been retrieved.
     *
     * Each step produced by this particles will be tagged
     * with the id set by "identifier".
     * 
     * Will throw if not initialized.
     */
    virtual void EnqueueLightSource(const I3CLSimLightSource &lightSource, I3CLSimStepFactoryPtr);
    
    /**
     * Adds a "barrier" to the particle queue. This will keep the
     * simulation running until all steps have been retrieved.
     * New particles can only be added to the queue after all
     * "old" steps have been retrieved.
     * 
     * Will throw if not initialized.
     */
    virtual void EnqueueBarrier();
    
    /**
     * Returns true an enqueued barrier is still active. And active
     * barrier means that no new particles can currently be added
     * to the queue. Steps have to be retrieved using GetConversionResult()
     * until this function returns false.
     * 
     * Will throw if not initialized.
     */
    virtual bool BarrierActive() const;
    
    /**
     * Returns true if more steps are available for the current particle.
     * If the return value is false, the current simulation is finished
     * and a new particle may be set.
     * 
     * Will throw if not initialized.
     */
    virtual bool MoreStepsAvailable() const;

    /**
     * Returns a bunch of steps stored in a vector<I3CLSimStep> or
     * a pair of <uint32_t, I3ParticleConstPtr>.
     * Returned particles were produced by the primary particle.
     * The uint32_t value is the primary's identifier as passed to
     * EnqueueLightSource().
     *
     * Might block if no steps are available.
     * The steps may belong to various particles at the same time.
     * 
     * Will throw if not initialized.
     */
    virtual I3CLSimStepSeriesConstPtr GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout=NAN);
    
    /**
     * Like GetConversionResultWithBarrierInfo(), but also returns
     * a vector of the identifiers that have finished processing.
     */
    virtual std::tuple<I3CLSimStepSeriesConstPtr, std::vector<uint32_t>, std::map<uint32_t, I3MCTreePtr>> GetConversionResultWithBarrierInfoAndMarkers(bool &barrierWasReset, double timeout=NAN);
    
    
private:

    typedef std::pair<I3CLSimStepFactoryPtr, I3CLSimLightSourceConstPtr> ToGeant4Pair_t;

    void WorkerThread();
    void WorkerThread_impl(boost::this_thread::disable_interruption &di);
    void StartThread();
    void StopThread();
    boost::shared_ptr<boost::thread> thread_;
    boost::condition_variable_any threadStarted_cond_;
    boost::mutex threadStarted_mutex_;
    bool threadStarted_;

    mutable boost::mutex barrier_is_enqueued_mutex_;
    bool barrier_is_enqueued_;

    boost::shared_ptr<I3CLSimQueue<ToGeant4Pair_t> > queueToGeant4_;
    boost::shared_ptr<I3CLSimQueue<FromGeant4Pair_t> > queueFromGeant4_;

    bool initialized_;
    uint64_t bunchSizeGranularity_;
    uint64_t maxBunchSize_;
    
    I3CLSimFunctionConstPtr wlenBias_;
    I3CLSimMediumPropertiesConstPtr mediumProperties_;
    
    std::vector<I3CLSimLightSourcePropagatorPtr> propagators_;
    
    SET_LOGGER("I3CLSimLightSourceToStepConverterAsync");
};

I3_POINTER_TYPEDEFS(I3CLSimLightSourceToStepConverterAsync);

#endif //I3CLSIMLIGHTSOURCETOSTEPCONVERTERASYNC_H_INCLUDED
