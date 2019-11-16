

#include <simclasses/I3MCTrajectoryTree.h>

#include <dataclasses/python/tree_indexing_suite.hpp>
#include <icetray/python/stream_to_string.hpp>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/operator_suite.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

void register_I3MCTrajectoryTree()
{
    class_<I3MCTrajectoryTree, I3MCTrajectoryTreePtr, bases<I3FrameObject> >(
        "I3MCTrajectoryTree", "A tree of lightweight particles used for simulation"
    )
        .def(init<const I3MCTrajectoryTree::value_type&>())
        .def(TreeBase::tree_indexing_suite<I3MCTrajectoryTree>())
        // dataclass suites
        .def(copy_suite<I3MCTrajectoryTree>())
        .def(operator_suite<I3MCTrajectoryTree>())
        .def_pickle(boost_serializable_pickle_suite<I3MCTrajectoryTree>())
        .def("__str__", &stream_to_string<I3MCTrajectoryTree>)
        .def("to_I3MCTree", &I3MCTrajectoryTree::to_I3MCTree)
        .def("clip", &I3MCTrajectoryTree::Clip)
    ;
    register_pointer_conversions<I3MCTrajectoryTree>();
}
