#ifndef SIMCLASSES_I3MCTRAJECTORY_H_INCLUDED
#define SIMCLASSES_I3MCTRAJECTORY_H_INCLUDED

#include <functional>

#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3ParticleID.h"
#include "dataclasses/physics/detail/I3MCTree_fwd.h"
#include "phys-services/surfaces/Surface.h"
#include "dataclasses/I3TimeWindow.h"

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
/// - A pointlike trajectory: The state after a call to SetPointlike().
///   GetDisplacement() returns 0, and GetShape() returns Cascade.
/// - A ranged trajectory: The state after one or more calls to AddPoint(). In
///   this state, the trajectory has one or more checkpoints where the time, energy,
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
          energy_(std::max(0., energy - GetMass())),
          state_(InitialState({dir.GetZenith(),dir.GetAzimuth()}))
    {
        auto id = I3ParticleID::create();
        majorID_ = id.majorID;
        minorID_ = id.minorID;
    }

    /// \cond
    typedef size_t size_type;
    /// \endcond

    class TrajectoryPoint {
    public:
        const I3Position& GetPos() const { return pos_; }
        double GetTime() const { return time_; }
        double GetKineticEnergy() const { return energy_; }
    private:
        friend class I3MCTrajectory;
        TrajectoryPoint() {}
        TrajectoryPoint(double time, double energy, const I3Position &pos)
            : time_(time), energy_(energy), pos_(pos) {}
        double time_, energy_;
        I3Position pos_;
    };

    /// \cond
    operator I3ParticleID() const { return I3ParticleID(majorID_,minorID_); }
    /// \endcond

    bool operator==(const I3MCTrajectory& rhs) const {
      return ( majorID_ == rhs.majorID_ &&
          minorID_ == rhs.minorID_ && 
          pdgEncoding_ == rhs.pdgEncoding_ &&
          CompareFloatingPoint::Compare_NanEqual(position_[0], rhs.position_[0]) &&
          CompareFloatingPoint::Compare_NanEqual(position_[1], rhs.position_[1]) &&
          CompareFloatingPoint::Compare_NanEqual(position_[2], rhs.position_[2]) &&
          CompareFloatingPoint::Compare_NanEqual(time_, rhs.time_) &&
          CompareFloatingPoint::Compare_NanEqual(energy_, rhs.energy_) &&
          state_ == rhs.state_
      );
    }
    bool operator!=(const I3MCTrajectory& rhs) const {
      return !(*this == rhs);
    }

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
    /// NB: pseudoparticles in the PDG code range above 2e6 are considered massless
    double GetMass() const
    {
        return (std::abs(GetType()) > 2000000000) ? 0 : I3Particle::GetMassForType(GetType());
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
    /// @brief Clip the trajectory to the given bounding volume
    ///
    /// @returns the portions of the trajectory that are inside the volume, if any
    boost::optional<I3MCTrajectory> Clip(const I3Surfaces::Surface &surface) const
    {
        boost::optional<I3MCTrajectory> result;
        boost::optional<I3TimeWindow> interval = FindCrossingTimes(surface);
        if (interval) {
            result = Clip(*interval);
        }
        return result;
    }
    /// @brief Clip the trajectory to the given time window
    ///
    /// @returns the portions of the trajectory that are within the given time window, if any
    boost::optional<I3MCTrajectory> Clip(const I3TimeWindow &interval) const
    {
        return boost::apply_visitor(visitors::Clip(*this, interval), state_);
    }
    /// @brief Find the times where a trajectory crosses a closed surface
    ///
    /// @returns the time window during which the trajectory is inside the volume
    boost::optional<I3TimeWindow> FindCrossingTimes(const I3Surfaces::Surface &surface) const
    {
        return boost::apply_visitor(visitors::FindCrossingTimes(*this, surface), state_);
    }
    /// @brief Remove intermediate trajectory points
    /// @tparam BinaryPredicate a binary predicate with a `bool operator(const &TrajectoryPoint, const &TrajectoryPoint)`. If the predicate returns `true`, the right-hand point is kept, otherwise it is removed.
    /// @returns a simplified I3MCTrajectory containing the vertex and endpoint, plus any intermediate points that pass the predicate
    template <class BinaryPredicate>
    I3MCTrajectory Simplify(BinaryPredicate&& pred) const
    {
        return boost::apply_visitor(visitors::make_Simplify(*this, pred), state_);
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
        bool operator==(const Checkpoint& rhs) const {
          return (
              CompareFloatingPoint::Compare_NanEqual(time, rhs.time) &&
              CompareFloatingPoint::Compare_NanEqual(energy, rhs.energy) &&
              CompareFloatingPoint::Compare_NanEqual(x, rhs.x) &&
              CompareFloatingPoint::Compare_NanEqual(y, rhs.y) &&
              CompareFloatingPoint::Compare_NanEqual(z, rhs.z)
          );
        }
    };
    struct InitialState {
        double zenith, azimuth;
        template <class Archive> void serialize(Archive& ar, unsigned version);
        bool operator==(const InitialState& rhs) const {
          return (
              CompareFloatingPoint::Compare_NanEqual(zenith, rhs.zenith) &&
              CompareFloatingPoint::Compare_NanEqual(azimuth, rhs.azimuth)
          );
        }
    };
    static_assert(2*sizeof(InitialState) < sizeof(I3Direction), "I3Direction is bloated");
    struct FinalState {
        double zenith, azimuth;
        template <class Archive> void serialize(Archive& ar, unsigned version);
        bool operator==(const FinalState& rhs) const {
          return (
              CompareFloatingPoint::Compare_NanEqual(zenith, rhs.zenith) &&
              CompareFloatingPoint::Compare_NanEqual(azimuth, rhs.azimuth)
          );
        }
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
        struct Clip : public boost::static_visitor<boost::optional<I3MCTrajectory>> {
            const I3MCTrajectory &self_;
            const I3TimeWindow &interval_;
            Clip(const I3MCTrajectory &self, const I3TimeWindow &interval) : self_(self), interval_(interval) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                result_type result;
                if (interval_.Contains(self_.GetTime(0)) || interval_.Contains(self_.GetTime(v.size()+1))) {
                    auto first = self_.GetIndexForTime(interval_.GetStart());
                    auto last = self_.GetIndexForTime(interval_.GetStop());
                    double first_trim = std::max(0., interval_.GetStart()-self_.GetTime(first));
                    double last_trim = std::max(0., self_.GetTime(last)-interval_.GetStop());

                    auto &clipped = *result;
                    // copy constant properties
                    clipped.majorID_ = self_.majorID_;
                    clipped.minorID_ = self_.minorID_;
                    clipped.SetType(self_.GetType());
                    clipped.energy_ = std::max(0., self_.GetKineticEnergy(first)-self_.GetMass());
                    // Shift vertex to the volume border
                    clipped.SetPos(self_.GetPos(first) + (first_trim*I3Constants::c*self_.GetBeta(first))*self_.GetDir(first));
                    clipped.SetTime(self_.GetTime(first) + first_trim);
                    
                    {
                        // Re-add intermediate checkpoints
                        // NB: segment i goes from point i to i+1
                        size_type i=first+1;
                        for ( ; i < last; i++) {
                            clipped.AddPoint(
                                self_.GetTime(i),
                                self_.GetKineticEnergy(i),
                                self_.GetPos(i)
                            );
                        }
                        if (last_trim > 0) {
                            // The last segment extends past the end; clip it.
                            clipped.AddPoint(
                                self_.GetTime(i) + last_trim,
                                self_.GetKineticEnergy(i),
                                self_.GetPos(i) + (last_trim*I3Constants::c*self_.GetBeta(last))*self_.GetDir(i)
                            );
                        } else {
                            // The last segment ends inside.
                            clipped.AddPoint(
                                self_.GetTime(i+1),
                                self_.GetKineticEnergy(i+1),
                                self_.GetPos(i+1)
                            );
                        }
                    }
                }

                return result;
            }
            template <typename T>
            result_type operator()(const T &v) const {
                result_type result;
                if (interval_.Contains(self_.GetTime())) {
                    result = self_;
                }
                return result;
            }
        };
        struct FindCrossingTimes : public boost::static_visitor<boost::optional<I3TimeWindow>> {
            const I3MCTrajectory &self_;
            const I3Surfaces::Surface &surface_;
            FindCrossingTimes(const I3MCTrajectory &self, const I3Surfaces::Surface &surface) : self_(self), surface_(surface) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                result_type result;
                size_type first(0), last(v.size());
                double first_trim(0), last_trim(0);
                // find the first segment that is not outside the surface
                // outside means: first > length || (first < 0 && second < 0)
                for ( ; first < last; first++) {
                    auto i = surface_.GetIntersection(self_.GetPos(first), self_.GetDir(first));
                    if (!((i.first > self_.GetLength(first)) || (i.first < 0 && i.second < 0))) {
                        first_trim = std::max(0., i.first);
                        break;
                    }
                }
                // find the last segment that is not outside the surface
                for (; last >= first; last--) {
                    auto i = surface_.GetIntersection(self_.GetPos(last), self_.GetDir(last));
                    if (!((i.first > self_.GetLength(last)) || (i.first < 0 && i.second < 0))) {
                        if (i.second < self_.GetLength(last))
                            last_trim = i.second;
                        break;
                    }
                }
                if (first != v.size()) {
                    result = I3TimeWindow(
                        self_.GetTime(first) + first_trim/(I3Constants::c*self_.GetBeta(first)),
                        last_trim > 0 ?
                            self_.GetTime(last) + last_trim/(I3Constants::c*self_.GetBeta(last))
                            : self_.GetTime(last+1)
                    );
                }

                return result;
            }
            result_type operator()(const InitialState &v) const {
                result_type result;
                auto intersections = surface_.GetIntersection(self_.GetVertexPosition(), self_.GetDir());
                if (intersections.first <= 0 && intersections.second >= 0) {
                    result = I3TimeWindow(self_.time_, std::numeric_limits<double>::infinity());
                }
                return result;
            }
            result_type operator()(const FinalState &v) const {
                result_type result;
                auto intersections = surface_.GetIntersection(self_.GetVertexPosition(), self_.GetDir());
                if (intersections.first <= 0 && intersections.second >= 0) {
                    result = I3TimeWindow(self_.time_, self_.time_);
                }
                return result;
            }
        };
        template <class BinaryPredicate>
        struct Simplify : public boost::static_visitor<I3MCTrajectory> {
            const I3MCTrajectory &self_;
            BinaryPredicate pred_;
            Simplify(const I3MCTrajectory &self, BinaryPredicate&& pred) : self_(self), pred_(pred) {};
            result_type operator()(const std::vector<Checkpoint> &v) const {
                result_type simplified;
                simplified.majorID_ = self_.majorID_;
                simplified.minorID_ = self_.minorID_;
                simplified.pdgEncoding_ = self_.pdgEncoding_;
                simplified.position_ = self_.position_;
                simplified.time_ = self_.time_;
                simplified.energy_ = self_.energy_;
                // Add intermediate points if they pass the predicate
                TrajectoryPoint left(self_.GetTime(0), self_.GetKineticEnergy(0), self_.GetPos(0));
                for (size_type i=1; i < v.size(); i++) {
                    TrajectoryPoint right(self_.GetTime(i), self_.GetKineticEnergy(i), self_.GetPos(i));
                    if (pred_(std::cref(left), std::cref(right))) {
                        simplified.AddPoint(
                            right.GetTime(),
                            right.GetKineticEnergy(),
                            right.GetPos()
                        );
                        left = right;
                    }
                }
                // Unconditionally add last point
                simplified.AddPoint(
                    self_.GetTime(v.size()+1),
                    self_.GetKineticEnergy(v.size()+1),
                    self_.GetPos(v.size()+1)
                );
                return simplified;
            }
            template <typename T>
            result_type operator()(const T &v) const {
                return result_type(self_);
            }
        };
        template <class BinaryPredicate>
        static Simplify<BinaryPredicate> make_Simplify(const I3MCTrajectory &self, BinaryPredicate&& pred) {
            return Simplify<BinaryPredicate>(self, pred);
        }
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
                    std::ostringstream oss;
                    oss << "Checkpoint "<<v.size()+1<<" has energy "<<energy_
                        <<" > checkoint "<<v.size()<<" ("<<v.back().energy<<")";
                    throw std::domain_error(oss.str());
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
            template <typename T>
            result_type operator()(T &v) const
            {
                if (!(reltime_ > 0)) {
                    throw std::domain_error("Checkpoint times must be strictly increasing");
                }
                if (!(energy_ <= self_.energy_)) {
                    throw std::domain_error("Checkpoint energy must be <= initial energy");
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
