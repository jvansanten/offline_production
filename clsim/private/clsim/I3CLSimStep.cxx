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
 * $Id: I3CLSimStep.cxx 177840 2019-12-11 18:50:26Z olivas $
 *
 * @file I3CLSimStep.cxx
 * @version $Revision: 177840 $
 * @date $Date: 2019-12-11 11:50:26 -0700 (Wed, 11 Dec 2019) $
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/I3CLSimStep.h>

#include <boost/static_assert.hpp>
#include <serialization/binary_object.hpp>

using namespace icecube::archive;

namespace {
    const std::size_t blobSizeV0 = 48; // size of our structure in bytes
    const std::size_t blobSizeV1 = 72; // size of our structure in bytes
    struct I3CLSimStepV0 {
        cl_float4 posAndTime;   // x,y,z,time
        cl_float4 dirAndLengthAndBeta; // theta,phi,length,beta
        cl_uint numPhotons;
        cl_float weight;
        cl_uint identifier;
        cl_uchar sourceType;
        cl_uchar dummy1;
        cl_ushort dummy2;
    };
}

I3CLSimStep::~I3CLSimStep() { }

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Waddress-of-packed-member"
#endif
template <class Archive>
void I3CLSimStep::save(Archive &ar, unsigned version) const
{
    ar << make_nvp("x", ((const cl_float *)&posAndTime)[0]);
    ar << make_nvp("y", ((const cl_float *)&posAndTime)[1]);
    ar << make_nvp("z", ((const cl_float *)&posAndTime)[2]);
    ar << make_nvp("time", ((const cl_float *)&posAndTime)[3]);
    
    ar << make_nvp("theta", ((const cl_float *)&dirAndLengthAndBeta)[0]);
    ar << make_nvp("phi", ((const cl_float *)&dirAndLengthAndBeta)[1]);
    ar << make_nvp("length", ((const cl_float *)&dirAndLengthAndBeta)[2]);
    ar << make_nvp("beta", ((const cl_float *)&dirAndLengthAndBeta)[3]);

    ar << make_nvp("num", numPhotons);
    ar << make_nvp("weight", weight);
    ar << make_nvp("id", identifier);
    ar << make_nvp("sourceType", sourceType);
    ar << make_nvp("dummy1", dummy1);
    ar << make_nvp("dummy2", dummy2);

    ar << make_nvp("rngKey", rngKey);
    ar << make_nvp("rngCounter", rngCounter);
}     

template <class Archive>
void I3CLSimStep::load(Archive &ar, unsigned version)
{
    if (version > i3clsimstep_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimStep class.",version,i3clsimstep_version_);

    float temp; uint32_t temp_uint; uint64_t temp_ulong;
    uint16_t temp_uint16; uint8_t temp_uint8;
    ar >> make_nvp("x", temp); ((cl_float *)&posAndTime)[0]=temp;
    ar >> make_nvp("y", temp); ((cl_float *)&posAndTime)[1]=temp;
    ar >> make_nvp("z", temp); ((cl_float *)&posAndTime)[2]=temp;
    ar >> make_nvp("time", temp); ((cl_float *)&posAndTime)[3]=temp;

    ar >> make_nvp("theta", temp); ((cl_float *)&dirAndLengthAndBeta)[0]=temp;
    ar >> make_nvp("phi", temp); ((cl_float *)&dirAndLengthAndBeta)[1]=temp;
    ar >> make_nvp("length", temp); ((cl_float *)&dirAndLengthAndBeta)[2]=temp;
    ar >> make_nvp("beta", temp); ((cl_float *)&dirAndLengthAndBeta)[3]=temp;
    ar >> make_nvp("num", temp_uint); numPhotons=temp_uint;
    ar >> make_nvp("weight", temp); weight=temp;
    ar >> make_nvp("id", temp_uint); identifier=temp_uint;
    ar >> make_nvp("sourceType", temp_uint8); sourceType=temp_uint8;
    ar >> make_nvp("dummy1", temp_uint8); dummy1=temp_uint8;
    ar >> make_nvp("dummy2", temp_uint16); dummy2=temp_uint16;
    if (version > 0) {
        // g++ refuses to take non-const references to packed members
        {
            auto temp = rngKey;
            ar >> make_nvp("rngKey", temp);
            rngKey = temp;
        }
        {
            auto temp = rngCounter;
            ar >> make_nvp("rngCounter", temp);
            rngCounter = temp;
        }
    }

}
#ifdef __clang__
#pragma clang diagnostic pop
#endif

// just save the binary blob for binary archives (internal storage is little-endian)

template <>
void I3CLSimStep::save(portable_binary_oarchive &ar, unsigned version) const
{
    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimStep) == blobSizeV1));

    ar << make_nvp("blob", icecube::serialization::make_binary_object((void *)this, blobSizeV1));
}

template <>
void I3CLSimStep::load(portable_binary_iarchive &ar, unsigned version)
{
    if (version > i3clsimstep_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimStep class.",version,i3clsimstep_version_);

    // check an assumption we will make throughout the code
    BOOST_STATIC_ASSERT((sizeof(I3CLSimStepV0) == blobSizeV0));
    BOOST_STATIC_ASSERT((sizeof(I3CLSimStep) == blobSizeV1));

    if (version == 0) {
        I3CLSimStepV0 temp;
        ar >> make_nvp("blob", icecube::serialization::make_binary_object(&temp, blobSizeV0));
        this->rngKey = {{0,0}};
        this->rngCounter = {{0,0,0,0}};
        this->posAndTime = temp.posAndTime;
        this->dirAndLengthAndBeta = temp.dirAndLengthAndBeta;
        this->numPhotons = temp.numPhotons;
        this->weight = temp.weight;
        this->identifier = temp.identifier;
        this->sourceType = temp.sourceType;
        this->dummy1 = temp.dummy1;
        this->dummy2 = temp.dummy2;
    } else {
        ar >> make_nvp("blob", icecube::serialization::make_binary_object(this, blobSizeV1));
    }
}

template<>
template<>
void I3Vector<I3CLSimStep>::serialize(portable_binary_iarchive &ar, unsigned version)
{
    ar >> make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    unsigned I3CLSimStep_version;
    ar >> make_nvp("I3CLSimStep_version", I3CLSimStep_version);
    if (I3CLSimStep_version != i3clsimstep_version_)
        log_fatal("This reader can only read I3Vector<I3CLSimStep> version %u, but %u was provided.",i3clsimstep_version_,I3CLSimStep_version);
    uint64_t size;
    ar >> make_nvp("num", size);

    this->resize(size);

    // read the binary blob in one go..
    if (I3CLSimStep_version == 0) {
        std::vector<I3CLSimStepV0> temp(size);
        ar >> make_nvp("blob", icecube::serialization::make_binary_object( temp.data(), blobSizeV0*size));
        for (size_t i = 0; i < size; i++) {
            (*this)[i].rngKey = {{0,0}};
            (*this)[i].rngCounter = {{0,0,0,0}};
            (*this)[i].posAndTime = temp[i].posAndTime;
            (*this)[i].dirAndLengthAndBeta = temp[i].dirAndLengthAndBeta;
            (*this)[i].numPhotons = temp[i].numPhotons;
            (*this)[i].weight = temp[i].weight;
            (*this)[i].identifier = temp[i].identifier;
            (*this)[i].sourceType = temp[i].sourceType;
            (*this)[i].dummy1 = temp[i].dummy1;
            (*this)[i].dummy2 = temp[i].dummy2;
        }
    } else {
        ar >> make_nvp("blob", icecube::serialization::make_binary_object( &((*this)[0]), blobSizeV1*size));
    }
}

template<>
template<>
void I3Vector<I3CLSimStep>::serialize(portable_binary_oarchive &ar, unsigned version)
{
    ar << make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar << make_nvp("I3CLSimStep_version", i3clsimstep_version_);
    uint64_t size = this->size();
    ar << make_nvp("num", size);
    ar << make_nvp("blob", icecube::serialization::make_binary_object( &((*this)[0]), blobSizeV1*size ));
}


#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winstantiation-after-specialization"
#endif
I3_SERIALIZABLE(I3CLSimStep);
I3_SERIALIZABLE(I3CLSimStepSeries);

#ifdef __clang__
#pragma clang diagnostic pop
#endif


std::ostream& operator<<(std::ostream& os, const I3CLSimStep& s){
    os << "[I3CLSimStep:\n"
       << "      Position: (" << s.GetPosX() << ',' << s.GetPosY() << ',' << s.GetPosZ() << ")\n"
       << "          Time: " << s.GetTime() << '\n'
       << "     Direction: (" << s.GetDirTheta() << ',' << s.GetDirPhi() << ")\n"
       << "        Length: " << s.GetLength() << '\n'
       << "          Beta: " << s.GetBeta() << '\n'
       << "  Num. Photons: " << s.GetNumPhotons() << '\n'
       << "        Weight: " << s.GetWeight() << '\n'
       << "            ID: " << s.GetID() << '\n'
       << "   Source Type: " << s.GetSourceType() << "\n]";
  return os;
}
