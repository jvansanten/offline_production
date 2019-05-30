
#ifndef I3CLSIMSTEPFACTORY_H_INCLUDED
#define I3CLSIMSTEPFACTORY_H_INCLUDED

#include "clsim/I3CLSimStep.h"
#include "phys-services/I3RandomService.h"

class I3CLSimStepFactory {
public:
    I3CLSimStepFactory(uint64_t frameKey, uint32_t lightSourceIndex, uint32_t lightSourceIdentifier);
    uint32_t GetLightSourceID() const { return lightSourceIdentifier_; }
    I3RandomServicePtr GetRandomStream() { return randomService_; };
    I3CLSimStep createStep();
private:
    /// A globally unique identifier for the source frame (used for RNG key)
    std::array<uint32_t,2> frameKey_;
    /// A locally unique identifier for the light source (used only to tie photons to frames)
    uint32_t lightSourceIdentifier_;
    /// Index of the light source within its frame (used for RNG counter)
    uint32_t lightSourceIndex_;
    /// Index of the step from this light source (used for RNG counter)
    uint32_t stepIndex_;
    I3RandomServicePtr randomService_;
};

I3_POINTER_TYPEDEFS(I3CLSimStepFactory);

#endif // I3CLSIMSTEPFACTORY_H_INCLUDED
