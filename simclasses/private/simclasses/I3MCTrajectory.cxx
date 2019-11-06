
#include "simclasses/I3MCTrajectory.h"
#include <icetray/serialization.h>
#include <serialization/variant.hpp>

I3MCTrajectory::I3MCTrajectory() {
    static_assert(sizeof(I3MCTrajectory) == 88, "I3MCTrajectory has no internal padding");
    static_assert(sizeof(I3MCTrajectory)*2 < sizeof(I3Particle), "I3MCTrajectory is less than half the size of I3Particle");
    static_assert(sizeof(I3MCTrajectory::Checkpoint) == 5*8, "I3MCTrajectory is less than half the size of I3Particle");
}



I3MCTrajectory::operator I3Particle() const
{
    I3Particle p(majorID_, minorID_);
    p.SetPos(GetPos());
    p.SetDir(GetDir());
    p.SetTime(GetTime());
    p.SetLength(GetDisplacement());
    p.SetType(GetType());
    p.SetKineticEnergy(GetKineticEnergy());
    p.SetSpeed(I3Constants::c*GetBeta());
    p.SetLocationType(I3Particle::InActiveVolume);
    p.SetFitStatus(I3Particle::NotSet);
    p.SetShape(GetShape());
    return p;
}

template <class Archive> 
void I3MCTrajectory::Checkpoint::serialize(Archive& ar, unsigned version)
{
    ar & make_nvp("time", time);
    ar & make_nvp("time", time);
    ar & make_nvp("x", x);
    ar & make_nvp("y", y);
    ar & make_nvp("z", z);
}

template <class Archive> 
void I3MCTrajectory::InitialState::serialize(Archive& ar, unsigned version)
{
    ar & make_nvp("zenith", zenith);
    ar & make_nvp("azimuth", azimuth);
}

template <class Archive> 
void I3MCTrajectory::FinalState::serialize(Archive& ar, unsigned version)
{
    ar & make_nvp("zenith", zenith);
    ar & make_nvp("azimuth", azimuth);
} 

template <class Archive> 
void I3MCTrajectory::serialize(Archive& ar, unsigned version)
{
    if (version>i3mctrajectory_version_){
        log_fatal("Attempting to read version %u from file but running "
            "version %u of I3MCTrajectory class.",
            version,i3mctrajectory_version_);
    }

    // ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("majorID", majorID_);
    ar & make_nvp("minorID", minorID_);
    ar & make_nvp("pdgEncoding", pdgEncoding_);
    ar & make_nvp("vertex", position_);
    ar & make_nvp("time", time_);
    ar & make_nvp("energy", energy_);
    ar & make_nvp("state", state_);
} 

I3_SERIALIZABLE(I3MCTrajectory);
