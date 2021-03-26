#ifndef CC_GUARD
#define CC_GUARD

#include "walker/walker.hpp"

#include <SeqLib/BamRecord.h>

#include <string>

namespace CC {

class cc_walker : public walker::walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);

   cc_walker(const std::string& bam_in, const std::string& interval_list) : walker(bam_in), interval_list_path(interval_list) {}

   protected:
   std::string interval_list_path;
   
   uint32_t curstart = 0;
   uint32_t curend = 0;
   uint16_t curchr = 0;
};

}

#endif
