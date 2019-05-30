
from os.path import expandvars
from icecube import clsim
from icecube import phys_services
from icecube.dataclasses import I3Particle, I3Direction, I3Position


def test_factory():
    """
    I3CLSimStepFactory creates steps with unique stepId
    """
    frameKey = 37 + (42 << 32)
    lightSourceIdentifier = 42
    lightSourceIndex = 3
    factory = clsim.I3CLSimStepFactory(frameKey, lightSourceIndex, lightSourceIdentifier)
    step = factory.createStep()
    assert step.rngKey == [frameKey & (1<<32)-1, frameKey>>32]
    assert step.rngCounter == [lightSourceIndex,0,0,0]
    step = factory.createStep()
    assert step.rngKey == [frameKey & (1<<32)-1, frameKey>>32]
    assert step.rngCounter == [lightSourceIndex,0,1,0]

def test_conversion():
    """
    I3CLSimLightSourceToStepConverterAsync creates steps with unique stepId
    """
    RandomService = phys_services.I3GSLRandomService(0)
    DetectorSettings = clsim.traysegments.common.setupDetector(
        expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
        UseGeant4=True,
    )
    stepGenerator = clsim.I3CLSimLightSourceToStepConverterAsync(1)
    stepGenerator.SetLightSourceParameterizationSeries(DetectorSettings['ParameterizationList'])
    stepGenerator.SetMediumProperties(DetectorSettings['MediumProperties'])
    stepGenerator.SetWlenBias(DetectorSettings['WavelengthGenerationBias'])
    try:
        stepGenerator.SetPropagators([clsim.I3CLSimLightSourcePropagatorGeant4()])
    except AttributeError:
        pass
    stepGenerator.Initialize()

    p = I3Particle(I3Position(0,0,0), I3Direction(1,0,0), 0, I3Particle.Cascade)
    p.type = p.EMinus
    p.energy = 5
    frameKey = 37 + (42 << 32)
    lightSourceIdentifier = 42
    lightSourceIndex = 3
    factory = clsim.I3CLSimStepFactory(frameKey, lightSourceIndex, lightSourceIdentifier)

    stepGenerator.EnqueueLightSource(clsim.I3CLSimLightSource(p), factory)
    stepGenerator.EnqueueBarrier()
    steps = stepGenerator.GetConversionResult()
    step_ids = set()
    for step in steps:
        if step.num > 0:
            assert step.rngKey == [frameKey & (1<<32)-1, frameKey>>32]
            assert step.rngCounter[0] == lightSourceIndex
            step_index = step.rngCounter[2]
            assert step_index not in step_ids
            step_ids.add(step_index)
    assert len(step_ids) > 1

if __name__ == "__main__":
    for name in sorted(filter(lambda k: k.startswith('test'), locals().keys())):
        locals()[name]()
