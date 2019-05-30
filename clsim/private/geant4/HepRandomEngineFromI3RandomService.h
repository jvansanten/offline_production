
#ifndef CLSIM_PHILOX4X32HEPRANDOMENGINE_INCLUDED
#define CLSIM_PHILOX4X32HEPRANDOMENGINE_INCLUDED

#include <CLHEP/Random/RandomEngine.h>
#include "icetray/I3Logging.h"
#include "phys-services/I3RandomService.h"

class HepRandomEngineFromI3RandomService : public CLHEP::HepRandomEngine {
public:
    HepRandomEngineFromI3RandomService(I3RandomServicePtr rng);
    virtual ~HepRandomEngineFromI3RandomService();
    virtual double flat();
    virtual void flatArray(const int size, double *vect);
    virtual std::string name() const { return "Philox4x32"; };
    virtual void setSeed(long seed, int) { log_fatal("unimplemented"); };
    virtual void setSeeds(const long* seeds, int) { log_fatal("unimplemented"); };
    virtual void saveStatus (const char filename[]="Config.conf") const { log_fatal("unimplemented"); };
    virtual void restoreStatus (const char filename[]="Config.conf") { log_fatal("unimplemented"); };
    virtual void showStatus () const { log_fatal("unimplemented"); };
    virtual std::ostream & put (std::ostream &os) const { return os; };
private:
    I3RandomServicePtr randomService_;

};

#endif //CLSIM_PHILOX4X32HEPRANDOMENGINE_INCLUDED