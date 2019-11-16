#ifndef SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
#define SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
#include "simclasses/I3MCTrajectory.h"
#include "dataclasses/physics/I3MCTree.h"

/**
 * I3MCTree - This goes into the frame and everyone can see it
 */
class I3MCTrajectoryTree : public TreeBase::Tree<I3MCTrajectory,I3ParticleID> {
public:
    using TreeBase::Tree<I3MCTrajectory,I3ParticleID>::Tree;

    explicit operator I3MCTree() const { return to_I3MCTree(); };
    I3MCTree to_I3MCTree() const;

    I3MCTrajectoryTree Clip(const I3Surfaces::Surface &surface) const;
private:
    friend class icecube::serialization::access;
    using TreeBase::Tree<I3MCTrajectory,I3ParticleID>::load;
    using TreeBase::Tree<I3MCTrajectory,I3ParticleID>::save;
    I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(I3MCTrajectoryTree,TreeBase::tree_version_);
I3_POINTER_TYPEDEFS(I3MCTrajectoryTree);
I3_DEFAULT_NAME(I3MCTrajectoryTree);

#endif //SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
