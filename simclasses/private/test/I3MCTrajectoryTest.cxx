
#include <I3Test.h>

#include "simclasses/I3MCTrajectory.h"

#include "dataclasses/physics/I3MCTree.h"
TEST_GROUP(I3MCTrajectory);

TEST(state_machine) {
    I3Position pos(1,2,3);
    I3Direction dir(1,0,0);
    double time = 0;
    double energy = 2*I3Units::MeV;
    I3MCTrajectory track(I3Particle::MuMinus, pos, dir, energy, time);
    {
        ENSURE_EQUAL(track.GetShape(), I3Particle::Null);
        ENSURE(std::isnan(track.GetLength()));
        ENSURE_EQUAL(track.GetNumSteps(), 0u);
        ENSURE_EQUAL(track.GetDir(), dir);
        ENSURE_EQUAL(track.GetKineticEnergy(), energy);
        ENSURE_EQUAL(track.GetTotalEnergy(), energy + I3Particle::GetMassForType(I3Particle::MuMinus));
    }

    {
        I3Position zig(0,1,0);
        I3Position zag(0,0,1);
        double new_energy = 1*I3Units::MeV;
        EXPECT_THROW(track.AddPoint(time, energy, pos + zig), "time offset must be positive");
        ENSURE_EQUAL(track.GetShape(), I3Particle::Null, "shape preserved after exception");
        EXPECT_THROW(track.AddPoint(time+1, energy+1, pos + zig), "particle may not gain energy");
        ENSURE_EQUAL(track.GetShape(), I3Particle::Null, "shape preserved after exception");

        track.AddPoint(1, new_energy, pos + zig);
        ENSURE_EQUAL(track.GetShape(), I3Particle::MCTrack);
        ENSURE_EQUAL(track.GetDisplacement(), 1);
        ENSURE_EQUAL(track.GetNumSteps(), 1u);
        ENSURE_EQUAL(track.GetLength(0), 1);
        ENSURE(std::isnan(track.GetLength(1)));
        ENSURE_EQUAL(track.GetDir(0), I3Direction(zig));
        ENSURE_EQUAL(track.GetDir(1), I3Direction(zig));
        ENSURE_EQUAL(track.GetKineticEnergy(0), energy);
        ENSURE_EQUAL(track.GetKineticEnergy(1), new_energy);
        ENSURE_EQUAL(track.GetTime(0), 0);
        ENSURE_EQUAL(track.GetTime(1), 1);

        track.AddPoint(2, new_energy, pos + zig + zag);
        ENSURE_EQUAL(track.GetShape(), I3Particle::MCTrack);
        ENSURE_EQUAL(track.GetDisplacement(), std::sqrt(2));
        ENSURE_EQUAL(track.GetNumSteps(), 2u);
        ENSURE_EQUAL(track.GetLength(0), 1);
        ENSURE_EQUAL(track.GetLength(1), 1);
        ENSURE(std::isnan(track.GetLength(2)));
        ENSURE_EQUAL(track.GetDir(0), I3Direction(zig));
        ENSURE_EQUAL(track.GetDir(1), I3Direction(zag));
        ENSURE_EQUAL(track.GetKineticEnergy(0), energy);
        ENSURE_EQUAL(track.GetKineticEnergy(1), new_energy);
        ENSURE_EQUAL(track.GetKineticEnergy(2), new_energy);
        ENSURE(track.GetBeta(1) < track.GetBeta(0));

        ENSURE_EQUAL(track.GetTime(0), 0);
        ENSURE_EQUAL(track.GetTime(1), 1);
        ENSURE_EQUAL(track.GetTime(2), 2);
    }

    {
        I3MCTrajectory clone(track);
        ENSURE_EQUAL(clone.GetShape(), I3Particle::MCTrack);
        clone.Clear();
        ENSURE_EQUAL(track.GetShape(), I3Particle::MCTrack);
        ENSURE_EQUAL(clone.GetShape(), I3Particle::Null);
    }

    {
        I3MCTrajectory clone(track);
        ENSURE_EQUAL(clone.GetShape(), I3Particle::MCTrack);
        clone.SetPointlike();
        ENSURE_EQUAL(track.GetShape(), I3Particle::MCTrack);
        ENSURE_EQUAL(clone.GetShape(), I3Particle::Cascade);
        ENSURE_EQUAL(clone.GetDisplacement(), 0);
        ENSURE_EQUAL(clone.GetNumSteps(), 0u);
    }


}

TEST(simplify) {
    I3Position pos(1,2,3);
    I3Direction dir(1,0,0);
    double time = 0;
    double energy = 2*I3Units::GeV;
    I3MCTrajectory track(I3Particle::MuMinus, pos, dir, energy, time);
    I3Position zig(0,1,0);
    I3Position zag(0,0,1);
    track.AddPoint(1, energy/2, pos + zig);
    track.AddPoint(2, energy/3, pos + zag);

    auto simplified = track.Simplify([](I3MCTrajectory::TrajectoryPoint a, I3MCTrajectory::TrajectoryPoint b) {
        return (a.GetPos()-b.GetPos()).Magnitude() > 30;
    });
    ENSURE_EQUAL(track.GetNumSteps(), 2u);
    ENSURE_EQUAL(simplified.GetNumSteps(), 1u);
}
