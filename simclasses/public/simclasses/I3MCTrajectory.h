#ifndef SIMCLASSES_I3MCTRAJECTORY_H_INCLUDED
#define SIMCLASSES_I3MCTRAJECTORY_H_INCLUDED
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3ParticleID.h"
#include "dataclasses/physics/detail/I3MCTree_fwd.h"

#include <boost/variant.hpp>

static const unsigned i3mctrajectory_version_ = 1;

/// @brief A lightweight representation of particle's path through the active
/// volume
///
/// A simulated particle initially has an identity, a particle type, a vertex
/// position and time, an energy, and a direction. The energy, time, position,
/// and direction may change as the particle propagates, but the identity and
/// and type do not. In the absence of a strong magnetic field, long-ranged
/// particles tend to travel on approximately straight lines between stochastic
/// energy losses. This leads to a natural representation of a trajectory as a
/// series of checkpoints (time,energy,x,y,z).
///
/// An I3MCTrajectory can represent one of three states:
/// - An initial state: The state at construction, or after a call to Clear().
///   It has a position, time, initial direction, and kinetic energy.
///   GetDisplacement() returns NAN, and GetShape() returns Null.
/// - A pointlike particle: The state after a call to SetPointlike().
///   GetDisplacement() returns 0, and GetShape() returns Cascade.
/// - A trajectory: The state after one or more calls to AddPoint(). In this
///   state, the trajectory has one or more checkpoints where the time, energy,
///   and position of the particle are known. In the segments between these
///   checkpoints, the particle momentum is constant. GetDisplacement() returns
///   the distance from the vertex to the last checkpoint, and GetShape()
///   returns MCTrack.
///
/// I3MCTrajectory is designed to be as small as possible: 88 bytes on the
/// stack, plus 40 on the heap for each segment. By contrast, I3Particles are
/// 200 bytes apiece (all on the stack), mostly due to bloated direction and
/// position classes, extra fields that are not generally useful for particles
/// in simulation, and internal padding.
///
/// The trajectory can be compressed further by reducing the number of
/// checkpoints, e.g. adding a checkpoint only when the energy changes by more
/// than 5%, the checkpoint is more than 20 m from the last one, the direction
/// deviates by more than 1 degree, etc.
class I3MCTrajectory {
public:
    /// @brief Initialize a trajectory
    /// The trajectory has a type, position, time, direction, and kinetic
    /// energy, but no length. GetShape() will return `I3Particle::Null`
    I3MCTrajectory(I3Particle::ParticleType type,
        const I3Position &pos, const I3Direction &dir,
        double energy, double time)
        : pdgEncoding_(type),
          position_({pos.GetX(),pos.GetY(),pos.GetZ()}),
          time_(time),
          energy_(energy),
          state_(InitialState({dir.GetZenith(),dir.GetAzimuth()}))
    {
        auto id = I3ParticleID::create();
        majorID_ = id.majorID;
        minorID_ = id.minorID;
    }

    /// \cond
    typedef size_t size_type;
    /// \endcond

    /// \cond
    operator I3ParticleID() const { return I3ParticleID(majorID_,minorID_); }
    /// \endcond

    I3Particle::ParticleType GetType() const { return I3Particle::ParticleType(pdgEncoding_); } 

    /// @brief Implicit conversion to I3Particle
    ///
    ///   - ID, Pos, Dir, KineticEnergy, Type are set from their corresponding
    ///     getters
    ///   - Length is set from GetDisplacment()
    ///   - Speed is set from GetBeta()
    ///   - LocationType is InActiveVolume
    ///   - FitStatus is NotSet
    operator I3Particle() const;

    /// Get the shape of the particle
    /// Null: has not yet been propagated
    /// MCTrack: has been propagated and has a trajectory
    /// Cascade: has been proapgated, but has not trajectory
    I3Particle::ParticleShape GetShape() const
    {
        return boost::apply_visitor(visitors::GetShape(), state_);
    }

    /// @brief Get the number of segments in the trajectory
    ///
    /// @returns a valid argument to any of the Get methods that take
    /// an index. It is >= 1 if GetShape() returns `I3Particle::MCTrack`, and 0
    /// otherwise.
    size_type GetNumSteps() const { return boost::apply_visitor(visitors::GetNumSteps(), state_); }

    /// @brief Get the index of the track segment at the given time
    ///
    /// @returns a valid argument to any of the Get methods that take
    /// an index. It is 0 if GetShape() returns something other than
    /// `I3Particle::MCTrack`.
    size_type GetIndexForTime(double time) const
    {
        return boost::apply_visitor(visitors::GetIndexForTime(time-time_), state_);
    }

    /// @param[in] index index of a segment
    /// @returns position of the beginning of the segment in I3Units
    I3Position GetPos(size_type index=0) const
    {
        return boost::apply_visitor(visitors::GetPos(*this, index), state_);
    }

    /// @param[in] index index of a segment
    /// @returns time of the beginning of the segment in I3Units
    double GetTime(size_type index=0) const
    {
        return boost::apply_visitor(visitors::GetTime(*this, index), state_);
    }

    /// @param[in] index index of a segment
    /// @returns direction at the beginning of a segment
    I3Direction GetDir(size_type index=0) const
    {
        return boost::apply_visitor(visitors::GetDir(*this, index), state_);
    }

    /// @param[in] index index of a segment
    /// @returns kinetic energy at the beginning of a segment in I3Units
    double GetKineticEnergy(size_type index=0) const
    {
        return boost::apply_visitor(visitors::GetEnergy(*this, index), state_);
    }
    
    /// @returns The particle's mass in I3Units
    double GetMass() const
    {
        return I3Particle::GetMassForType(GetType());
    }

    /// @param[in] index index of the segment
    /// @returns The Lorentz beta of the segment
    double GetBeta(size_type index=0) const
    {
        double beta_squared = 1.-1./std::pow(1. + GetKineticEnergy(index)/GetMass(),2);
        if (beta_squared > 0) {
            return std::sqrt(beta_squared);
        } else {
            return 0.;
        }
    }

    /// @param[in] index index of the segment
    /// @returns The length of the segment in I3Units. If `index` >=
    /// GetNumSteps(), returns NAN
    double GetLength(size_type index=0) const
    {
        return boost::apply_visitor(visitors::GetLength(*this, index), state_);
    }

    /// @returns distance from vertex to endpoint in I3Units.
    double GetDisplacement() const
    {
        return boost::apply_visitor(visitors::GetDisplacement(), state_);
    }

    /// @brief Set vertex position
    /// @param[in] pos The new vertex position
    /// Since the trajectory points are relative to the vertex, this moves them
    /// as well.
    void SetPos(I3Position &&pos) { SetPos(pos); }
    /// @copydoc SetPos
    void SetPos(const I3Position &pos) {
        position_[0] = pos.GetX();
        position_[1] = pos.GetY();
        position_[2] = pos.GetZ();
    }

    /// @brief Set vertex time
    /// @param[in] time The new vertex time
    /// Since the trajectory times are relative to the vertex time, this shifts
    /// them as well.
    void SetTime(double time) { time_ = time; }

    /// @brief Set the particle type
    void SetType(I3Particle::ParticleType type) { pdgEncoding_ = type; }
    /// @copydoc SetType
    void SetPdgEncoding(int32_t type) { pdgEncoding_ = type; }

    /// @brief Add a point to the trajectory
    ///
    /// @param[in] time   time, on the same scale as the vertex time
    /// @param[in] energy kinetic energy
    /// @param[in] pos    kosition, in same coordinate system as the vertex
    ///
    /// If Shape == MCTrack, this appends a point to the end of the trajectory.
    /// Otherwise, it sets the shape to MCTrack and inserts the point as the
    /// first.
    void AddPoint(double time, double energy, const I3Position &pos)
    {
        double reltime = time - time_;
        I3Position displacement = pos - GetVertexPosition();
        return boost::apply_visitor(visitors::AddPoint(*this,reltime,energy,displacement), state_);
    }
    /// @brief Mark as a point-like trajectory.
    ///
    /// Point-like trajectories a direction, but zero length. This changes the
    /// shape to Cascade.
    void SetPointlike()
    {
        return boost::apply_visitor(visitors::SetPointlike(*this), state_);
    }
    /// @brief Reset the trajectory to an initial state
    /// 
    /// Initial trajectories have a direction, but _no_ length. This removes
    /// any trajectory points that may have been added previously.
    void Clear()
    {
        return boost::apply_visitor(visitors::Clear(*this), state_);
    }

private:
    /// \cond
    I3MCTrajectory();
    I3Position GetVertexPosition() const { return I3Position(position_[0],position_[1],position_[2]); }

    struct Checkpoint {
        double time;
        double energy;
        double x,y,z;
        I3Position GetPos() const { return I3Position(x,y,z); }
        template <class Archive> void serialize(Archive& ar, unsigned version);
    };
    struct InitialState {
        double zenith, azimuth;
        template <class Archive> void serialize(Archive& ar, unsigned version);
    };
    static_assert(2*sizeof(InitialState) < sizeof(I3Direction), "I3Direction is bloated");
    struct FinalState {
        double zenith, azimuth;
        template <class Archive> void serialize(Archive& ar, unsigned version);
    };

    // NB: we store the components of I3ParticleID separately to prevent
    // the compiler from adding 4 bytes of padding before and after pdgEncoding_
    // Inherting from I3FrameObject would also add 8 bytes (4 for the virtual
    // table pointer, and 4 more to pad out for the alignment of majorID_)
    uint64_t majorID_;
    uint32_t minorID_;
    int32_t pdgEncoding_;
    /// @brief Vertex position
    /// NB: I3Position contains representations of itself in 3 different
    /// coordinate systems. We store the Cartesian 3-vector here to save space.
    std::array<double,3> position_;
    static_assert(2*sizeof(position_) < sizeof(I3Position), "I3Position is bloated");
    /// Vertex time
    double time_;
    /// Initial kinetic energy
    double energy_;

    boost::variant<InitialState,FinalState,std::vector<Checkpoint> > state_;
    static_assert(sizeof(state_) == 32, "state variant is small enough");

    struct visitors {

        struct GetNumSteps : public boost::static_visitor<size_type> {
            result_type operator()(const std::vector<Checkpoint> &v) const { return v.size(); }
            result_type operator()(const InitialState &v) const { return 0; }
            result_type operator()(const FinalState &v) const { return 0; }
        };
        struct GetDisplacement : public boost::static_visitor<double> {
            result_type operator()(const std::vector<Checkpoint> &v) const { return v.back().GetPos().Magnitude(); }
            result_type operator()(const InitialState &v) const { return NAN; }
            result_type operator()(const FinalState &v) const { return 0; }
        };
        struct GetShape : public boost::static_visitor<I3Particle::ParticleShape> {
            result_type operator()(const std::vector<Checkpoint> &v) const { return I3Particle::MCTrack; }
            result_type operator()(const InitialState &v) const { return I3Particle::Null; }
            result_type operator()(const FinalState &v) const { return I3Particle::Cascade; }
        };
        struct GetPos : public boost::static_visitor<I3Position> {
            const I3MCTrajectory &self_;
            size_type index_;
            GetPos(const I3MCTrajectory &self, size_type i=0) : self_(self), index_(i) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                if (index_ == 0) {
                    return self_.GetVertexPosition();
                } else {
                    return self_.GetVertexPosition() + v[std::min(index_,v.size())-1].GetPos();
                }
            }
            result_type operator()(const InitialState &v) const { return self_.GetVertexPosition(); }
            result_type operator()(const FinalState &v) const { return self_.GetVertexPosition(); }
        };
        struct GetLength : public boost::static_visitor<double> {
            const I3MCTrajectory &self_;
            size_type index_;
            GetLength(const I3MCTrajectory &self, size_type i=0) : self_(self), index_(i) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                if (index_ >= v.size()) {
                    // We can't know where the segment ends
                    return NAN;
                } else if (index_ < 1) {
                    // NB: checkpoint positions are relative to the vertex
                    return v[0].GetPos().Magnitude();
                } else {
                    size_type i = std::min(index_,v.size()-1);
                    return (v[i].GetPos()-v[i-1].GetPos()).Magnitude();
                }
            }
            result_type operator()(const InitialState &v) const { return NAN; }
            result_type operator()(const FinalState &v) const { return 0; }
        };
        struct GetDir : public boost::static_visitor<I3Direction> {
            const I3MCTrajectory &self_;
            size_type index_;
            GetDir(const I3MCTrajectory &self, size_type i=0) : self_(self), index_(i) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                if (index_ < 1 || v.size() < 2) {
                    return I3Direction(v[0].GetPos());
                } else {
                    size_type i = std::min(index_,v.size()-1);
                    return I3Direction(v[i].GetPos()-v[i-1].GetPos());
                }
            }
            result_type operator()(const InitialState &v) const { return I3Direction(v.zenith, v.azimuth); }
            result_type operator()(const FinalState &v) const { return I3Direction(v.zenith, v.azimuth); }
        };
        struct GetTime : public boost::static_visitor<double> {
            const I3MCTrajectory &self_;
            size_type index_;
            GetTime(const I3MCTrajectory &self, size_type i=0) : self_(self), index_(i) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                if (index_ == 0) {
                    return self_.time_;
                } else {
                    return self_.time_+v[std::min(index_,v.size())-1].time;
                }
            }
            result_type operator()(const InitialState &v) const { return self_.time_; }
            result_type operator()(const FinalState &v) const { return self_.time_; }
        };
        struct Clear : public boost::static_visitor<> {
            I3MCTrajectory &self_;
            Clear(I3MCTrajectory &self) : self_(self) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                I3Direction dir = self_.GetDir();
                self_.state_ = InitialState({dir.GetZenith(), dir.GetAzimuth()});
            }
            result_type operator()(const InitialState &v) const { }
            result_type operator()(const FinalState &v) const {
                self_.state_ = InitialState({v.zenith, v.azimuth});
            }
        };
        struct SetPointlike : public boost::static_visitor<> {
            I3MCTrajectory &self_;
            SetPointlike(I3MCTrajectory &self) : self_(self) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                I3Direction dir = self_.GetDir();
                self_.state_ = FinalState({dir.GetZenith(), dir.GetAzimuth()});
            }
            result_type operator()(InitialState &v) const {
                self_.state_ = FinalState({v.zenith, v.azimuth});
            }
            result_type operator()(const FinalState &v) const { }
        };
        struct AddPoint : public boost::static_visitor<> {
            I3MCTrajectory &self_;
            double reltime_;
            double energy_;
            const I3Position &displacement_;
            AddPoint(I3MCTrajectory &self, double reltime, double energy, const I3Position &displacement)
                : self_(self), reltime_(reltime), energy_(energy), displacement_(displacement) {};
            result_type operator()(std::vector<Checkpoint> &v) const {
                if (!(reltime_ > v.back().time)) {
                    throw std::domain_error("Checkpoint times must be strictly increasing");
                }
                if (!(energy_ <= v.back().energy)) {
                    throw std::domain_error("Checkpoint energy must be <= previous checkpoint");
                }
                v.push_back(
                    {
                        reltime_,
                        energy_,
                        displacement_.GetX(),
                        displacement_.GetY(),
                        displacement_.GetZ()
                    });
            }
            result_type operator()(InitialState &v) const { InitializeTrajectory(); }
            result_type operator()(FinalState &v) const { InitializeTrajectory(); }
            result_type InitializeTrajectory() const
            {
                if (!(reltime_ > 0)) {
                    throw std::domain_error("Checkpoint times must be strictly increasing");
                }
                if (!(energy_ <= self_.energy_)) {
                    throw std::domain_error("Checkpoint energy must be <= previous checkpoint");
                }
                self_.state_ = std::vector<Checkpoint>(
                    1,
                    {
                        reltime_,
                        energy_,
                        displacement_.GetX(),
                        displacement_.GetY(),
                        displacement_.GetZ()
                    });
            }
        };
        struct GetIndexForTime : public boost::static_visitor<size_type> {
            double reltime_;
            GetIndexForTime(double reltime) : reltime_(reltime) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                return std::distance(v.begin(),
                    std::lower_bound(v.begin(), v.end(), reltime_,
                    [](const Checkpoint &segment, double time) { return segment.time < time; }));
            }
            result_type operator()(const InitialState &v) const { return 0; }
            result_type operator()(const FinalState &v) const { return 0; }
        };
        struct GetEnergy : public boost::static_visitor<double> {
            const I3MCTrajectory &self_;
            size_type index_;
            GetEnergy(const I3MCTrajectory &self, size_type i=0) : self_(self), index_(i) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                if (index_ == 0) {
                    return self_.energy_;
                } else {
                    return v[std::min(index_,v.size())-1].energy;
                }
            }
            result_type operator()(const InitialState &v) const { return self_.energy_; }
            result_type operator()(const FinalState &v) const { return self_.energy_; }
        };

    };

    friend class icecube::serialization::access;
    friend class TreeBase::Tree<I3MCTrajectory,I3ParticleID>;
    friend class TreeBase::TreeNode<I3MCTrajectory>;
    template <class Archive> void serialize(Archive& ar, unsigned version);
    /// \endcond
};

std::ostream& operator<<(std::ostream& oss, const I3MCTrajectory& t);

I3_POINTER_TYPEDEFS(I3MCTrajectory);

I3_CLASS_VERSION(I3MCTrajectory,i3mctrajectory_version_);

#endif //SIMCLASSES_I3MCTRAJECTORY_H_INCLUDED
