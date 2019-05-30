#!/usr/bin/env python

"""
Test whether cascade extension can be disabled in the PPC parameterizations
"""

from icecube import icetray, dataclasses, clsim, phys_services

p = dataclasses.I3Particle()
p.pos = dataclasses.I3Position(0,0,0)
p.dir = dataclasses.I3Direction(0,0,1)
p.time = 0
p.energy = 1
p.type = p.EMinus

source = clsim.I3CLSimLightSource(p)

converter = clsim.I3CLSimLightSourceToStepConverterPPC()
converter.SetWlenBias(clsim.GetIceCubeDOMAcceptance())
converter.SetMediumProperties(clsim.MakeIceCubeMediumProperties())
converter.Initialize()

factory = clsim.I3CLSimStepFactory(37, 0, 0)
converter.EnqueueLightSource(source, factory)
steps = converter.GetConversionResult()
assert all( [s.pos.z > 0 for s in steps] ), "Steps are spread along the cascade by default"

converter.SetUseCascadeExtension(False)
converter.EnqueueLightSource(source, factory)
steps = converter.GetConversionResult()
assert all( [s.pos.z == 0 for s in steps] ), "Steps are all at the origin"
