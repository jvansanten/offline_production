
#include <icetray/serialization.h>
#include "simclasses/I3MCTrajectoryTree.h"

namespace TreeBase {
template<>
std::ostream& Tree<I3MCTrajectory,I3ParticleID>::Print(std::ostream& s) const{
  s << "[I3MCTrajectoryTree:\n";
  I3Position pos;
  I3Direction dir;
  for(Tree::const_iterator iter=cbegin(),end=cend(); iter!=end; iter++){
    for(unsigned int d=0; d<depth(*iter)+1; d++)
      s << "  ";
    s << I3ParticleID(*iter).minorID << " " << I3Particle::GetTypeString(iter->GetType()) << " ";
    pos = iter->GetPos();
    dir = iter->GetDir();
    s << "(" << pos.GetX()/I3Units::m << "m, ";
    s << pos.GetY()/I3Units::m << "m, ";
    s << pos.GetZ()/I3Units::m << "m) ";
    s << "(" << dir.GetZenith()/I3Units::degree << "deg, ";
    s << dir.GetAzimuth()/I3Units::degree << "deg) ";
    s << iter->GetTime()/I3Units::ns << "ns ";
    s << iter->GetKineticEnergy()/I3Units::GeV << "GeV ";
    s << iter->GetDisplacement()/I3Units::m << "m\n";
  }
  s << ']';
  return(s);
}
} //namespace TreeBase

I3_SERIALIZABLE(I3MCTrajectoryTree);
