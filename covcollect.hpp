#ifndef CC_GUARD
#define CC_GUARD

#include "walker/walker.hpp"

#include <SeqLib/BamRecord.h>
#include <SeqLib/GenomicRegionCollection.h>

#include <string>
#include <unordered_map>

namespace CC {

typedef struct r {
   uint32_t start;
   uint32_t end;
   uint32_t ref_consumed;
} read_boundary;

class cc_walker : public walker::walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);
   void load_intervals(uint32_t pad);
   void walk_all();

   cc_walker(const std::string& bam_in, const std::string& interval_list) : walker(bam_in), interval_list_path(interval_list) {}

   protected:
   std::string interval_list_path;
   SeqLib::GenomicRegionCollection<> intervals;

   std::unordered_map<std::string, read_boundary> read_cache;

   uint32_t curstart = 0;
   uint32_t curend = 0;
   uint16_t curchr = 0;
};

}

#endif
