/**
 * @file I3Philox4x32.cxx
 * @copyright (c) 2019 The IceCube Collaboration
 * @author Jakob van Santen
 * @date March 2019
 */

#include "clsim/random_generator/I3Philox4x32.h"

I3Philox4x32::I3Philox4x32(
    I3Philox4x32::BitGeneratorType::ctr_type counter,
    I3Philox4x32::BitGeneratorType::key_type key
) : engine_(counter, key)
{}

I3Philox4x32::~I3Philox4x32(){}
