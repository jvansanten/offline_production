
#include "clsim/I3CLSimStepFactory.h"
#include "clsim/random_generator/I3Philox4x32.h"

I3CLSimStepFactory::I3CLSimStepFactory(uint64_t frameKey, uint32_t lightSourceIndex, uint32_t lightSourceIdentifier)
    : lightSourceIdentifier_(lightSourceIdentifier), lightSourceIndex_(lightSourceIndex), stepIndex_(0)
{
    frameKey_[0] = frameKey & ((1ul<<32)-1);
    frameKey_[1] = frameKey >> 32;
    // use stream 1 for sequential operations. this could in principle be
    // extended to multiple streams, given a deterministic partioning scheme
    randomService_ = I3RandomServicePtr(new I3Philox4x32({{lightSourceIndex_,1,0,0}}, {{frameKey_[0], frameKey_[1]}}));
}

I3CLSimStep I3CLSimStepFactory::createStep()
{
    I3CLSimStep step;
    step.rngKey = frameKey_;
    step.rngCounter = {{
        lightSourceIndex_,
        0,
        stepIndex_++,
        0
    }};
    step.SetID(lightSourceIdentifier_);
    return step;
}
