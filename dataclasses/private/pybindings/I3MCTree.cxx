//
//   Copyright (c) 2004, 2005, 2006, 2007   Troy D. Straszheim  
//   
//   $Id: I3MCTree.cxx 175121 2019-08-20 16:50:35Z kjmeagher $
//
//   This file is part of IceTray.
//
//   IceTray is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   IceTray is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <vector>
#include <algorithm>
#include <exception>

#include <dataclasses/physics/I3MCTree.h>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/operator_suite.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/stream_to_string.hpp>
#include <dataclasses/python/tree_indexing_suite.hpp>

using namespace boost::python;

void register_I3MCTree()
{
  class_<I3MCTree, bases<I3FrameObject>, I3MCTreePtr>("I3MCTree",
    "A tree of I3Particles used for simulation.")
    // Extra Constructors
    .def(init<const I3Particle&>())

    // I3MCTreeUtils
    .def("get_daughters", &I3MCTreeUtils::GetDaughters, "Get all daughters/children of an I3ParticleID")
    .def("has_parent", &I3MCTreeUtils::HasParent, "Does the I3ParticleID have a parent?")
    .def("has", &I3MCTreeUtils::Has, "Does the I3ParticleID exist in the tree?")
    .def("add_primary", &I3MCTreeUtils::AddPrimary, "Add an I3Particle as a primary (at root-level)")
    .def("get_primary", &I3MCTreeUtils::GetPrimary, "Get the primary that created the I3ParticleID")
    .def("get_primaries", &I3MCTreeUtils::GetPrimaries, "Get a list of all primaries")
    .add_property("primaries", &I3MCTreeUtils::GetPrimaries)
    .def("append_child", &I3MCTreeUtils::AppendChild, "Add a child to an I3ParticleID")
    .def("get_particle", &I3MCTreeUtils::GetParticle, "Get the I3Particle represented by the I3ParticleID")
    .def("dump", &I3MCTreeUtils::Dump, "Return tree as a string")

    // base class methods
    .def(TreeBase::tree_indexing_suite<I3MCTree>())
    // dataclass suites
    .def(copy_suite<I3MCTree>())
    .def(operator_suite<I3MCTree>())
    .def_pickle(boost_serializable_pickle_suite<I3MCTree>())
    .def("__str__", &stream_to_string<I3MCTree>)
  ;
  register_pointer_conversions<I3MCTree>();
}

