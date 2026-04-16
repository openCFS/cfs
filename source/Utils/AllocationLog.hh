#ifndef ALLOCATION_LOG_HH
#define ALLOCATION_LOG_HH

<<<<<<< HEAD
#include <boost/atomic/atomic.hpp>
#include <unordered_map>
#include <mutex>
=======
#include <mutex>
#include <boost/unordered/unordered_flat_map.hpp>
>>>>>>> origin/master
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{ 
/** Central log (for static usage) for thread save counting of allocations.
 * Meant for Matrix.hh with the assumption that there are not too many different sizes.
 * All is thread-save! */     
struct AllocationLog
{
<<<<<<< HEAD
  public:
  /** this comes not for free and shall be guarded by progOpts->DoDetailedInfo() */
  inline void AddAllocation(unsigned int size)
  {
    // we use the mutex only if the key needs to be created
    auto it = map_.find(size);
    if(it == map_.end())
    {
      std::lock_guard<std::mutex> lock(mutex_);
      map_[size] = 1;
    }
    else 
    {
      it->second.fetch_add(1);        
    }  
=======
  AllocationLog() { map_.reserve(100); } 

  /** this comes not for free and shall be guarded by progOpts->DoDetailedInfo() */
  inline void AddAllocation(unsigned int size)
  {
    std::scoped_lock lock(mutex_);
    map_[size]++; // the first access to map creates 0
>>>>>>> origin/master
  }

  void ToInfo(const PtrParamNode& in)
  {
<<<<<<< HEAD
=======
    std::scoped_lock lock(mutex_);
>>>>>>> origin/master
    if(!map_.empty())
    {
      // practically it is sufficient to call this only for double
      PtrParamNode pn = in->Get("variants");
<<<<<<< HEAD
      // for performance we have an unordered map but we output ordered
=======
      // for performance reasons we have an unordered map but we output ordered
>>>>>>> origin/master
      std::map<unsigned int, unsigned int> sorted(map_.begin(), map_.end());
      unsigned total = 0;
      for(auto it : sorted) {
       total += it.second;
       pn->Get("size_" + std::to_string(it.first))->SetValue(it.second);
      }
      in->Get("allocations/total")->SetValue(total);
    }
    else
    {
      in->Get("info")->SetValue("no allocations logged, run cfs with -d (detailed info)");
    }
  }

  private:
  /** we write in .info.xml how often we resize which size to have a base for optimization
   * Note that the static stuff is for each template type but usually we want only double and complex */
<<<<<<< HEAD
  std::unordered_map<unsigned int, boost::atomic<unsigned int>> map_;
=======
  boost::unordered_flat_map<unsigned int, unsigned int> map_;
>>>>>>> origin/master
  /** guard for allocation_map */
  std::mutex mutex_;
}; // end of AllocationLog

extern AllocationLog matrixLog;


} //end of namespace
#endif // ALLOCATION_LOG_HH