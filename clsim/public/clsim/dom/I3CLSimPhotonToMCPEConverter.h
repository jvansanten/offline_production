
#ifndef CLSIM_I3CLSimPhotonToMCPEConverter_h_INCLUDED
#define CLSIM_I3CLSimPhotonToMCPEConverter_h_INCLUDED

#include <icetray/OMKey.h>
#include <icetray/I3PointerTypedefs.h>
#include <simclasses/I3MCPE.h>
#include <boost/optional.hpp>

class ModuleKey;
class I3CompressedPhoton;
I3_FORWARD_DECLARATION(I3RandomService);

class I3CLSimPhotonToMCPEConverter {
public:
    virtual ~I3CLSimPhotonToMCPEConverter();
    virtual boost::optional<std::tuple<OMKey,I3MCPE>> Convert(I3RandomService &randomService, const ModuleKey&, const I3CompressedPhoton &) const = 0;
    boost::optional<std::tuple<OMKey,I3MCPE>> Convert(I3RandomServicePtr &&randomService, const ModuleKey& key, const I3CompressedPhoton &photon)
    { return Convert(*randomService, key, photon); }
};

#endif
