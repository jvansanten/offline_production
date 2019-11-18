

#include <simclasses/I3MCTrajectory.h>
#include <boost/preprocessor/seq.hpp>

#include <icetray/python/list_indexing_suite.hpp>

using namespace boost::python;

static std::string 
I3MCTrajectory_prettyprint(const I3MCTrajectory& p)
{
    std::ostringstream oss;
    I3Position pos = p.GetPos();
    I3Direction dir = p.GetDir();
    oss << "[ I3MCTrajectory " << std::endl
        << "              ID : " << I3ParticleID(p) << std::endl
        << "        zen, azi : " << dir.GetZenith()/I3Units::deg  << ", " << dir.GetAzimuth()/I3Units::deg << " deg" << std::endl
        << "           x,y,z : " << pos.GetX() << ", " << pos.GetY() << ", " << pos.GetZ() << std::endl
        << "            time : " << p.GetTime() << std::endl
        << "     kin. energy : " << p.GetKineticEnergy() << std::endl
        << "          length : " << p.GetDisplacement() << std::endl
        << "            type : " << I3Particle::GetTypeString(p.GetType()) << std::endl
        << "           shape : " << I3Particle::GetShapeString(p.GetShape()) << std::endl
        << "           steps : " << p.GetNumSteps() << std::endl
        << "]" ;
    return oss.str();
}

static I3MCTrajectory Simplify(const I3MCTrajectory &t, object predicate)
{
    return t.Simplify(predicate);
}

template <typename Target>
static object Clip(const I3MCTrajectory &t, const Target &target)
{
    object result;
    auto optional = t.Clip(target);
    if (optional) {
        result = object(*optional);
    }
    return result;
}

void register_I3MCTrajectory()
{
    {
        scope outer = 
        class_<I3MCTrajectory, I3MCTrajectoryPtr>(
            "I3MCTrajectory",
            init<
                I3Particle::ParticleType,
                const I3Position&,
                const I3Direction&,
                double,
                double
            >((arg("type"),"pos","dir","energy","time"))
        )
        .add_property("type", &I3MCTrajectory::GetType, &I3MCTrajectory::SetType)
        .add_property("shape", &I3MCTrajectory::GetShape)
        .def("__len__", &I3MCTrajectory::GetNumSteps)
        .def("__str__", &I3MCTrajectory_prettyprint)
        .def("__eq__", &I3MCTrajectory::operator==)
        .def("GetKineticEnergy", &I3MCTrajectory::GetKineticEnergy, arg("index")=0)
        .def("GetTime", &I3MCTrajectory::GetTime, arg("index")=0)
        .def("GetPos", &I3MCTrajectory::GetPos, arg("index")=0)
        .def("GetDir", &I3MCTrajectory::GetDir, arg("index")=0)
        .def("Simplify", &Simplify, arg("predicate"))
        .def("Clip", &Clip<I3Surfaces::Surface>, arg("surface"))
        .def("Clip", &Clip<I3TimeWindow>, arg("interval"))
        ;

        class_<I3MCTrajectory::TrajectoryPoint>("TrajectoryPoint", no_init)
            .add_property("kinetic_energy", &I3MCTrajectory::TrajectoryPoint::GetKineticEnergy)
            .add_property("time", &I3MCTrajectory::TrajectoryPoint::GetTime)
            .add_property("pos", make_function(&I3MCTrajectory::TrajectoryPoint::GetPos, return_value_policy<copy_const_reference>()))
        ;
    }


    implicitly_convertible<I3MCTrajectory,I3ParticleID>();

    class_<std::vector<I3MCTrajectory>, boost::shared_ptr<std::vector<I3MCTrajectory>>>("I3MCTrajectorySeries")
    .def(list_indexing_suite<std::vector<I3MCTrajectory>>())
    ;
}
