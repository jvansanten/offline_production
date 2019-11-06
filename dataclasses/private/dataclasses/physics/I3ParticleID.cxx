
#include "dataclasses/physics/I3ParticleID.h"

#if BOOST_VERSION >= 105300
#define USE_ATOMICS
#endif

#ifdef USE_ATOMICS
#include <boost/thread/lock_guard.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/atomic.hpp>
#endif

#include <unistd.h>
#include <cstdlib>

#include <boost/functional/hash/hash.hpp>

#ifdef USE_ATOMICS
namespace{
	boost::mutex global_id_lock;
	
	boost::atomic<int32_t> global_last_pid_(0);
	boost::atomic<int32_t> global_minor_id_(0);
	boost::atomic<uint64_t> global_major_id_(0);
}
#else
static int32_t global_last_pid_ = 0;
static int32_t global_minor_id_ = 0;
static uint64_t global_major_id_ = 0;
#endif

I3ParticleID I3ParticleID::create(){
  int this_pid = getpid();
  I3ParticleID ID;
  
#ifdef USE_ATOMICS //thread-safe version
  int32_t last_pid=global_last_pid_.load(boost::memory_order_relaxed);
  boost::atomic_thread_fence(boost::memory_order_acquire); //keep memory ops from wandering
  if (this_pid != last_pid) {
    boost::lock_guard<boost::mutex> lg(global_id_lock); //acquire the lock
    //check whether another thread already updated this
    last_pid=global_last_pid_.load(boost::memory_order_relaxed);
    if(this_pid != last_pid){
      log_debug("PID has changed from %i to %i. regenerating I3Particle::majorID.", last_pid, this_pid);
      boost::atomic_thread_fence(boost::memory_order_release);
      global_last_pid_.store(this_pid, boost::memory_order_relaxed);
      global_major_id_.store(0, boost::memory_order_relaxed); // this will cause a new major ID to be generated
      global_minor_id_.store(0, boost::memory_order_relaxed); // reset the minor ID, too
    }
  }
  
  uint64_t old_major_id=global_major_id_.load(boost::memory_order_relaxed);
  boost::atomic_thread_fence(boost::memory_order_acquire); //keep memory ops from wandering
  if(old_major_id==0){
    boost::lock_guard<boost::mutex> lg(global_id_lock); //acquire the lock
    //check whether another thread already updated this
    old_major_id=global_major_id_.load(boost::memory_order_relaxed);
    if(old_major_id==0){
      boost::hash<std::string> string_hash;
      std::stringstream s;
      s<<time(0)<<this_pid<<gethostid();
      boost::atomic_thread_fence(boost::memory_order_release);
      global_major_id_.store(string_hash(s.str()), boost::memory_order_relaxed);
    }
  }
  ID.majorID = global_major_id_.load();
  ID.minorID = global_minor_id_.fetch_add(1);
#else //unsafe version if atomics aren't available
  if (this_pid != global_last_pid_) {
    log_debug("PID has changed from %i to %i. regenerating I3Particle::majorID.", global_last_pid_, this_pid);
    global_last_pid_ = this_pid;
    global_major_id_ = 0; // this will cause a new major ID to be generated
    global_minor_id_ = 0; // reset the minor ID, too
  }
  if(global_major_id_ ==0){
    boost::hash<std::string> string_hash;
    std::stringstream s;
    s<<time(0)<<this_pid<<gethostid();
    global_major_id_ = string_hash(s.str());
  }
  ID.majorID = global_major_id_;
  ID.minorID = global_minor_id_++;
#endif

  return ID;
}
