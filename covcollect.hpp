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
} read_boundary_t;

typedef struct t {
   uint32_t n_corrected;
   uint32_t n_uncorrected;
} target_counts_t;

class cc_walker : public walker::walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);
   void load_intervals(uint32_t pad);
   void walk_all();
   uint32_t n_overlap(const SeqLib::GenomicRegion& region, uint32_t start, uint32_t end);

   cc_walker(const std::string& bam_in, const std::string& interval_list) : walker(bam_in), interval_list_path(interval_list) {}

   protected:
   std::string interval_list_path;
   SeqLib::GenomicRegionCollection<> intervals;
   target_counts_t target_coverage = {0, 0};
   size_t cur_region_idx = 0;
   SeqLib::GenomicRegion cur_region;

   std::unordered_map<std::string, read_boundary_t> read_cache;

   uint32_t curstart = 0;
   uint32_t curend = 0;
   uint16_t curchr = 0;
};

}

#endif
