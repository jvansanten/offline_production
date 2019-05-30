/**
 * @file I3Philox4x32.h
 * @copyright (c) 2019 The IceCube Collaboration
 * @author Jakob van Santen
 * @date March 2019
 */

#ifndef I3PHILOX4X32_H
#define I3PHILOX4X32_H

#include "icetray/I3ServiceFactory.h"
#include "Random123/philox.h"
#include "Random123/MicroURNG.hpp"
#include "phys-services/I3StdRandomEngine.h"

/**
 * @class I3Philox4x32
 * @brief An implementation of the I3RandomService interface using the C++ 
 * random number engine for I3Philox4x32
 *  
 */ 
class I3Philox4x32 : public I3StdRandomEngine<I3Philox4x32>
{
 public:
  typedef r123::MicroURNG<r123::Philox4x32> BitGeneratorType;
  I3Philox4x32(BitGeneratorType::ctr_type, BitGeneratorType::key_type);

  /**
   * destructor
   */
  virtual ~I3Philox4x32();
  
private:
  BitGeneratorType engine_;
public:
  BitGeneratorType& engine() { return engine_; }
  const BitGeneratorType& engine() const { return engine_; }

  SET_LOGGER("I3Philox4x32");
};

I3_POINTER_TYPEDEFS(I3Philox4x32);

#endif //I3PHILOX4X32_H
