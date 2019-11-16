
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

/// @brief Convert to legacy I3MCTree
I3MCTree I3MCTrajectoryTree::to_I3MCTree() const {
    I3MCTree tree;
    for (auto t = this->cbegin(); t != this->cend(); t++) {
        auto parent = this->parent(t);
        if (parent == this->cend()) {
            tree.insert_after(*t);
        } else {
            // Keep siblings ordered in time
            auto sibling = tree.children(tree.find(*parent));
            while (sibling != tree.end_sibling() && !(t->GetTime() < sibling->GetTime()))
                sibling++;
            if (sibling == tree.end_sibling()) {
                tree.append_child(*parent, *t);
            } else {
                tree.insert(sibling, *t);
            }
        }
        if (t->GetShape() != I3Particle::MCTrack)
            continue;
        // Emulate I3MuonSlicer by marking the track as Dark and adding
        // constant-energy slices as daughters.
        auto p = tree.find(*t);
        p->SetShape(I3Particle::Dark);
        for (I3MCTrajectory::size_type i = 0 ; i < t->GetNumSteps(); i++) {
            auto segment = tree.append_child(p, I3Particle(
                t->GetPos(i),
                t->GetDir(i),
                t->GetTime(i),
                I3Particle::ContainedTrack,
                t->GetLength(i)
            ));
            segment->SetKineticEnergy(t->GetKineticEnergy(i));
            segment->SetSpeed(I3Constants::c*t->GetBeta(i));
            segment->SetLocationType(I3Particle::InActiveVolume);
            segment->SetType(t->GetType());
        }
    }

    return tree;
}

I3_SERIALIZABLE(I3MCTrajectoryTree);
