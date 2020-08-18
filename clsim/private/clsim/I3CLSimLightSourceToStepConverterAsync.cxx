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
 * $Id: I3CLSimLightSourceToStepConverterAsync.cxx 179360 2020-03-10 16:07:35Z eganster $
 *
 * @file I3CLSimLightSourceToStepConverterAsync.cxx
 * @version $Revision: 179360 $
 * @date $Date: 2020-03-10 10:07:35 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#include "clsim/I3CLSimLightSourceToStepConverterAsync.h"
#include "clsim/I3CLSimLightSourcePropagator.h"

#include "clsim/I3CLSimQueue.h"

#include <boost/thread/locks.hpp>
#include <boost/foreach.hpp>

#include <limits>
#include <deque>
#include <boost/tuple/tuple.hpp>
#include <cmath>

#include "clsim/I3CLSimStepStore.h"

// other headers
#include <stdlib.h>


// static definitions

const uint32_t I3CLSimLightSourceToStepConverterAsync::default_maxQueueItems=5;

I3CLSimLightSourceToStepConverterAsync::I3CLSimLightSourceToStepConverterAsync(uint32_t maxQueueItems)
:
queueToGeant4_(new I3CLSimQueue<ToGeant4Pair_t>(maxQueueItems)),
queueFromGeant4_(new I3CLSimQueue<FromGeant4Pair_t>(0)),
initialized_(false),
bunchSizeGranularity_(512),
maxBunchSize_(512000)
{
}

I3CLSimLightSourceToStepConverterAsync::~I3CLSimLightSourceToStepConverterAsync()
{
    StopThread();
}

void I3CLSimLightSourceToStepConverterAsync::StopThread()
{

    if (thread_)
    {
        if (thread_->joinable())
        {
            log_debug("Stopping the worker thread..");

            thread_->interrupt();
            
            thread_->join(); // wait for it indefinitely


            log_debug("Worker thread stopped.");
        }
        
        thread_.reset();
    }
}

void I3CLSimLightSourceToStepConverterAsync::StartThread()
{
    if (thread_)
        I3CLSimLightSourceToStepConverter_exception("Thread is already started!");
    log_debug("Starting the worker thread..");
    threadStarted_=false;
    barrier_is_enqueued_=false;

    thread_ = boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&I3CLSimLightSourceToStepConverterAsync::WorkerThread, this)));

    // wait for startup
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        for (;;)
        {
            if (threadStarted_) break;
            threadStarted_cond_.wait(guard);
        }
    }        
    
    log_debug("Worker thread started.");

}

void I3CLSimLightSourceToStepConverterAsync::Initialize()
{
    if (initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync already initialized!");

    if (!wlenBias_)
        throw I3CLSimLightSourceToStepConverter_exception("WlenBias not set!");

    if (!mediumProperties_)
        throw I3CLSimLightSourceToStepConverter_exception("MediumProperties not set!");
    
    if (bunchSizeGranularity_ > maxBunchSize_)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity must not be greater than MaxBunchSize!");
    
    if (maxBunchSize_%bunchSizeGranularity_ != 0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize is not a multiple of BunchSizeGranularity!");

    // stop the worker thread
    StopThread();

    // distribute medium properties and bias factors
    const I3CLSimLightSourceParameterizationSeries &parameterizations = this->GetLightSourceParameterizationSeries();
    for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimLightSourceParameterization &parameterization = *it;
        if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
        
        parameterization.converter->SetMediumProperties(mediumProperties_);
        parameterization.converter->SetWlenBias(wlenBias_);
        parameterization.converter->SetBunchSizeGranularity(1); // we do not send the bunches directly, the steps are integrated in the step store first, so granularity does not matter
        parameterization.converter->SetMaxBunchSize(maxBunchSize_); // use the same bunch size for the parameterizations
    }

    // now initialize them, skipping those that are already initialized
    for (I3CLSimLightSourceParameterizationSeries::const_iterator it=parameterizations.begin();
         it!=parameterizations.end(); ++it)
    {
        const I3CLSimLightSourceParameterization &parameterization = *it;

        // skip initialized converters
        if (!parameterization.converter->IsInitialized())
            parameterization.converter->Initialize();
    }

    // do the same two-phase initialization for propagators
    for (auto &propagator : propagators_) {
        if (!propagator) log_fatal("Internal error: NULL propagator");
        propagator->SetMediumProperties(mediumProperties_);
        propagator->SetWlenBias(wlenBias_);
    }
    for (auto &propagator : propagators_) {
        if (!propagator->IsInitialized())
            propagator->Initialize();
    }

    // [re]start the worker thread
    StartThread();

    initialized_=true;
}

void I3CLSimLightSourceToStepConverterAsync::WorkerThread()
{
    // do not interrupt this thread by default
    boost::this_thread::disable_interruption di;

    try {
        WorkerThread_impl(di);
    } catch(...) { // any exceptions?
        log_error("Worker thread died unexpectedly..");
        // If we do not throw here, I3CLSimClientModule will wait forever
        throw;
    }
}

void I3CLSimLightSourceToStepConverterAsync::WorkerThread_impl(boost::this_thread::disable_interruption &di)
{
    // set up geant4
    
    // this thing stores all the steps generated, sorted by their number of
    // Cherenkov photons
    I3CLSimStepStore stepStore(0);
    
    // this stores all particles that will be ent to parametrizations
    typedef std::tuple<I3CLSimLightSourceConstPtr, uint32_t, I3CLSimLightSourceParameterization&> sendToParameterizationQueueItem;
    std::deque<sendToParameterizationQueueItem> sendToParameterizationQueue;
    
    // this stores the ids of light sources we are still processing
    std::deque<uint32_t> markers;

    std::map<uint32_t, I3MCTreePtr> particleHistories;

    // notify the main thread that everything is set up
    {
        boost::unique_lock<boost::mutex> guard(threadStarted_mutex_);
        threadStarted_=true;
    }
    threadStarted_cond_.notify_all();
    
    // make a copy of the list of available parameterizations
    I3CLSimLightSourceParameterizationSeries parameterizations = this->GetLightSourceParameterizationSeries();

    // Push steps to output queue
    auto flushStepStore = [&](bool resetBarrier=false) {
        // Flush full-sized bunches
        while (stepStore.size() > maxBunchSize_)
        {
            I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
            stepStore.pop_bunch_to_vector(maxBunchSize_, *steps);
            
            // emit markers for any _previous_ light sources whose output
            // is no longer in the store
            std::vector<uint32_t> finished;
            std::map<uint32_t, I3MCTreePtr> daughters;
            while (!markers.empty() && stepStore.count(markers.front()) == 0) {
                finished.push_back(markers.front());
                auto target = particleHistories.find(markers.front());
                if (target != particleHistories.end()) {
                    daughters.emplace(std::move(*target));
                    particleHistories.erase(target);
                }
                markers.pop_front();
            }
            
            {
                boost::this_thread::restore_interruption ri(di);
                // this blocks if the queue from Geant4 to
                // OpenCL is full.
                queueFromGeant4_->Put(std::make_tuple(steps, std::move(finished), std::move(daughters), false));
            }
        }
        
        // Flush a partial bunch if requested
        if (resetBarrier) {
            
            // when flushing the last steps for an event, it may be necessary
            // to include dummy steps to fill the number of steps to be a multiple
            // of bunchSizeGranularity_. This is a step with weight==0 and numPhotons==0,
            // so it should not contribute to the final results!
            I3CLSimStep NoOpStepTemplate;
            NoOpStepTemplate.SetPos(I3Position(0.,0.,0.));
            NoOpStepTemplate.SetDir(I3Direction(0.,0.,-1.));
            NoOpStepTemplate.SetTime(0.);
            NoOpStepTemplate.SetLength(0.);
            NoOpStepTemplate.SetNumPhotons(0);
            NoOpStepTemplate.SetWeight(0.);
            NoOpStepTemplate.SetBeta(1.);
            
            const std::size_t numStepsWithDummyFill = bunchSizeGranularity_>1?(((stepStore.size()/bunchSizeGranularity_)+1)*bunchSizeGranularity_):stepStore.size();
            
            I3CLSimStepSeriesPtr steps(new I3CLSimStepSeries());
            stepStore.pop_bunch_to_vector(numStepsWithDummyFill, *steps, NoOpStepTemplate);
            
            if (!stepStore.empty())
                log_fatal("Internal logic error. step store should be empty.");
            
            // emit all remaining markers
            {
                boost::this_thread::restore_interruption ri(di);
                queueFromGeant4_->Put(std::make_tuple(steps, std::vector<uint32_t>(markers.begin(), markers.end()), particleHistories, true /* this is the LAST reply before the barrier is reached! */));
            }
            markers.clear();
            particleHistories.clear();
        }
    };
    
    // Add a step to the store, and push bunches if needed
    std::function<void(const I3CLSimStep &step)> emitStep = [&](const I3CLSimStep &step) {
        stepStore.insert_copy(step.GetNumPhotons(), step);
        flushStepStore();
    };
    
    // Call a parameterization to get steps
    auto getStepsFromParameterization = [&](I3CLSimLightSourceParameterization &parameterization, I3CLSimLightSourceConstPtr &lightSource, I3CLSimStepFactoryPtr stepFactory) {
        // call the converter
        if (!parameterization.converter) log_fatal("Internal error: parameteriation has NULL converter");
        if (!parameterization.converter->IsInitialized()) log_fatal("Internal error: parameterization converter is not initialized.");
        if (parameterization.converter->BarrierActive()) log_fatal("Logic error: parameterization converter has active barrier.");
                
        parameterization.converter->EnqueueLightSource(*lightSource, stepFactory);
        parameterization.converter->EnqueueBarrier();
        
        // get steps from the parameterization until the barrier is reached
        for (;;)
        {
            I3CLSimStepSeriesConstPtr res;
            bool barrierHasBeenReached=false;
            
            {
                boost::this_thread::restore_interruption ri(di);
                // this blocks if there are no steps yet and the
                // parameterization code is still working.
                res = parameterization.converter->GetConversionResultWithBarrierInfo(barrierHasBeenReached);
            }
            
            if (!res) {
                log_debug("NULL result from parameterization GetConversionResult(). ignoring.");
            } else {
                // add steps from the parameterization to the step store
                BOOST_FOREACH(const I3CLSimStep &step, *res)
                {
                    stepStore.insert_copy(step.GetNumPhotons(), step);
                }
            }
            
            // we don't need the results anymore
            res.reset();
            
            // push steps out if there are enough of them
            flushStepStore();

            if (barrierHasBeenReached) break; // get out of the loop if the barrier has been reached
        }
    };
    
    // Push a particle into the step generation stack
    std::function<bool(I3CLSimLightSourceConstPtr &, I3CLSimStepFactoryPtr, I3CLSimLightSourcePropagatorPtr)> addLightSource =
        [&](I3CLSimLightSourceConstPtr &lightSource, I3CLSimStepFactoryPtr stepFactory, I3CLSimLightSourcePropagatorPtr source) ->bool
    {
        for (auto &parameterization : parameterizations) {
            if (parameterization.IsValidForLightSource(*lightSource))
            {
                getStepsFromParameterization(parameterization, lightSource, stepFactory);
                return true;
            }
        }
        
        for (auto &propagator : propagators_) {
            if (propagator != source && propagator->IsValidForLightSource(*lightSource))
            {
                // Prevent this propagator from consuming its own output
                namespace ph = std::placeholders;
                I3CLSimLightSourcePropagator::secondary_callback emitSecondary =
                    std::bind<bool>(addLightSource, ph::_1, ph::_2, propagator);
                I3MCTreePtr daughters = propagator->Convert(lightSource, stepFactory, emitSecondary, emitStep);
                if (daughters) {
                    particleHistories[stepFactory->GetLightSourceID()] = daughters;
                }
                return true;
            }
        }
        
        return false;
    };

    // start the main loop
    for (;;)
    {
        I3CLSimLightSourceConstPtr lightSource;
        I3CLSimStepFactoryPtr stepFactory;

        {
            boost::this_thread::restore_interruption ri(di);
            try {
                std::tie(stepFactory, lightSource) = queueToGeant4_->Get();
            }
            catch(boost::thread_interrupted &i)
            {
                log_debug("G4 thread was interrupted. closing.");
                break;
            }
        }
        
        // Flush any full-sized bunches from the store, plus the last partial
        // bunch if the barrier was reached
        try {
            flushStepStore(!lightSource);
        } catch (boost::thread_interrupted &i) {
            log_debug("G4 thread was interrupted. closing.");
            break;
        }
        
        // If the barrier was reached, there's nothing left to do
        if (!lightSource)
            continue;

        if (lightSource->GetType() == I3CLSimLightSource::Unknown)
        {
            log_warn("Ignoring a light source with type \"Unknown\".");
            continue;
        }
        
        if ((lightSource->GetType() == I3CLSimLightSource::Particle) && (lightSource->GetParticle().GetType() == I3Particle::unknown))
        {
            log_warn("Ignoring a particle with type \"unknown\".");
            continue;
        }
        
        try {
            if (!addLightSource(lightSource, stepFactory, NULL))
                log_fatal_stream("No propagator or parameterization can handle this "<<lightSource->GetParticle().GetTypeString());
        } catch (boost::thread_interrupted &i) {
            log_debug("Thread was interrupted. closing.");
            break;
        }
        
        // Light source is eligible for finalization after the next bunch
        markers.push_back(stepFactory->GetLightSourceID());
        
    }

    log_debug("Worker thread terminated.");
}

bool I3CLSimLightSourceToStepConverterAsync::IsInitialized() const
{
    return initialized_;
}

void I3CLSimLightSourceToStepConverterAsync::SetBunchSizeGranularity(uint64_t num)
{
    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("BunchSizeGranularity of 0 is invalid!");

    StopThread();
    bunchSizeGranularity_=num;
    initialized_=false;
}

void I3CLSimLightSourceToStepConverterAsync::SetMaxBunchSize(uint64_t num)
{
    if (num<=0)
        throw I3CLSimLightSourceToStepConverter_exception("MaxBunchSize of 0 is invalid!");

    StopThread();
    maxBunchSize_=num;
    initialized_=false;
}

void I3CLSimLightSourceToStepConverterAsync::SetWlenBias(I3CLSimFunctionConstPtr wlenBias)
{
    StopThread();
    wlenBias_=wlenBias;
    initialized_=false;
}

void I3CLSimLightSourceToStepConverterAsync::SetMediumProperties(I3CLSimMediumPropertiesConstPtr mediumProperties)
{
    StopThread();
    mediumProperties_=mediumProperties;
    initialized_=false;
}

void I3CLSimLightSourceToStepConverterAsync::SetPropagators(const std::vector<I3CLSimLightSourcePropagatorPtr> &propagators)
{
    StopThread();
    propagators_ = propagators;
    initialized_=false;
}

void I3CLSimLightSourceToStepConverterAsync::EnqueueLightSource(const I3CLSimLightSource &lightSource, I3CLSimStepFactoryPtr stepFactory)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimLightSourceToStepConverter_exception("A barrier is enqueued! You must receive all steps before enqueuing a new particle.");
    }
    
    I3CLSimLightSourceConstPtr lightSourceCopy(new I3CLSimLightSource(lightSource));
    queueToGeant4_->Put(std::make_pair(stepFactory, lightSourceCopy));
}

void I3CLSimLightSourceToStepConverterAsync::EnqueueBarrier()
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
        if (barrier_is_enqueued_)
            throw I3CLSimLightSourceToStepConverter_exception("A barrier is already enqueued!");
        
        barrier_is_enqueued_=true;

        // we use a NULL pointer as the barrier
        queueToGeant4_->Put(std::make_pair(I3CLSimStepFactoryPtr(), I3CLSimLightSourceConstPtr()));
    }
}

bool I3CLSimLightSourceToStepConverterAsync::BarrierActive() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync is not initialized!");

    {
        boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);

        return barrier_is_enqueued_;
    }
}

bool I3CLSimLightSourceToStepConverterAsync::MoreStepsAvailable() const
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync is not initialized!");

    return (!queueFromGeant4_->empty());
}

std::tuple<I3CLSimStepSeriesConstPtr, std::vector<uint32_t>, std::map<uint32_t, I3MCTreePtr>>
I3CLSimLightSourceToStepConverterAsync::GetConversionResultWithBarrierInfoAndMarkers(bool &barrierWasReset, double timeout)
{
    if (!initialized_)
        throw I3CLSimLightSourceToStepConverter_exception("I3CLSimLightSourceToStepConverterAsync is not initialized!");

    barrierWasReset=false;
    
    FromGeant4Pair_t ret;
    if (!std::isnan(timeout))
        ret = queueFromGeant4_->Get(timeout/I3Units::second, FromGeant4Pair_t(nullptr, {}, {}, false));
    else
        ret = queueFromGeant4_->Get(); // no timeout
    
    if (std::get<3>(ret))
    {
        {
            boost::unique_lock<boost::mutex> guard(barrier_is_enqueued_mutex_);
            if (!barrier_is_enqueued_)
                log_error("Internal logic error. Barrier is not set as enqueued, yet a finalization message was received.");
            barrierWasReset=true;
            barrier_is_enqueued_=false;
        }
    }
    return std::make_tuple(std::move(std::get<0>(ret)), std::move(std::get<1>(ret)), std::move(std::get<2>(ret)));
}

I3CLSimStepSeriesConstPtr I3CLSimLightSourceToStepConverterAsync::GetConversionResultWithBarrierInfo(bool &barrierWasReset, double timeout)
{
    I3CLSimStepSeriesConstPtr steps;
    
    std::tie(steps, std::ignore, std::ignore) = GetConversionResultWithBarrierInfoAndMarkers(barrierWasReset, timeout);
    
    return steps;
}

