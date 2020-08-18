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
 * $Id: I3PhotonToMCPEConverter.h 179373 2020-03-10 21:25:48Z jvansanten $
 *
 * @file I3PhotonToMCPEConverter.h
 * @version $Revision: 179373 $
 * @date $Date: 2020-03-10 15:25:48 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#ifndef I3PHOTONTOMCPECONVERTER_H_INCLUDED
#define I3PHOTONTOMCPECONVERTER_H_INCLUDED

#include "icetray/I3Module.h"
#include "icetray/I3ConditionalModule.h"

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"
#include "dataclasses/physics/I3MCTree.h"
#include "simclasses/I3MCPE.h"
#include <simclasses/I3ParticleIDMap.hpp>

#include "phys-services/I3RandomService.h"

#include "clsim/function/I3CLSimFunction.h"
#include "clsim/dom/I3CLSimPhotonToMCPEConverter.h"

#include <string>

class I3CompressedPhoton;
class ModuleKey;

/**
 * @brief This module reads I3PhotonSeriesMaps generated
 * by CLSim, applies (D)OM acceptances (wavelength&angular)
 * to the photons and stores the results in an I3MCPESeriesMap.
 *
 */
class I3PhotonToMCPEConverter : public I3ConditionalModule
{
public:
    /**
     * Builds an instance of this class
     */
    I3PhotonToMCPEConverter(const I3Context& ctx);
    
    /**
     * Destroys an instance of this class
     */
    virtual ~I3PhotonToMCPEConverter();
    
    /**
     * This module takes a configuration parameter and so it must be configured.
     */
    virtual void Configure();
    
    /**
     * The module needs to process Physics frames
     */
    void DAQ(I3FramePtr frame);

    virtual void Finish();

    
private:
    // parameters
    
    /// Parameter: A random number generating service (derived from I3RandomService).
    boost::shared_ptr<I3RandomService> randomService_;

    boost::shared_ptr<I3CLSimPhotonToMCPEConverter> mcpeGenerator_;

    /// Parameter: Name of the input I3PhotonSeriesMap frame object. 
    std::string inputPhotonSeriesMapName_;

    /// Parameter: Name of the output I3MCPESeries frame object. 
    std::string outputMCPESeriesMapName_;

    /// Parameter: Compress the output I3MCPEs
    bool mergeHits_;

private:
    // default, assignment, and copy constructor declared private
    I3PhotonToMCPEConverter();
    I3PhotonToMCPEConverter(const I3PhotonToMCPEConverter&);
    I3PhotonToMCPEConverter& operator=(const I3PhotonToMCPEConverter&);
    
    template <typename PhotonMapType>
    std::pair<I3MCPESeriesMapPtr,I3ParticleIDMapPtr> Convert(I3FramePtr frame);
    
    // record some statistics
    uint64_t numGeneratedHits_;
    
    SET_LOGGER("I3PhotonToMCPEConverter");
};

class I3CLSimPhotonToMCPEConverterForDOMs : public I3CLSimPhotonToMCPEConverter {
public:
    I3CLSimPhotonToMCPEConverterForDOMs(boost::shared_ptr<const std::map<OMKey, I3CLSimFunctionConstPtr>>, I3CLSimFunctionConstPtr);
    virtual ~I3CLSimPhotonToMCPEConverterForDOMs();
    virtual boost::optional<std::tuple<OMKey,I3MCPE>> Convert(I3RandomService &randomService, const ModuleKey&, const I3CompressedPhoton &) const override;
private:
    boost::shared_ptr<const std::map<OMKey, I3CLSimFunctionConstPtr>> wavelengthAcceptance_;
    I3CLSimFunctionConstPtr angularAcceptance_;
};

#endif //I3PHOTONTOMCPECONVERTER_H_INCLUDED
