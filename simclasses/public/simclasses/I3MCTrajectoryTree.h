#ifndef SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
#define SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
#include "simclasses/I3MCTrajectory.h"
#include "dataclasses/physics/I3MCTree.h"

/**
 * I3MCTree - This goes into the frame and everyone can see it
 */
typedef TreeBase::Tree<I3MCTrajectory,I3ParticleID> I3MCTrajectoryTree;

I3_CLASS_VERSION(I3MCTrajectoryTree,TreeBase::tree_version_);
I3_POINTER_TYPEDEFS(I3MCTrajectoryTree);
I3_DEFAULT_NAME(I3MCTrajectoryTree);

#endif //SIMCLASSES_I3MCTRAJECTORYTREE_H_INCLUDED
