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
   uint32_t pad;

   std::unordered_map<std::string, read_boundary_t> read_cache;

   uint32_t curstart = 0;
   uint32_t curend = 0;
   uint16_t curchr = 0;
};



class cc_bin_walker : public walker::walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);
   void load_intervals(uint32_t pad);
   void set_binwidth(uint32_t binwidth);
   void walk_all();
   uint32_t n_overlap(const uint32_t binstart, uint32_t binend, uint32_t start, uint32_t end);

   cc_bin_walker(const std::string& bam_in, const uint32_t binwidth) : walker(bam_in), binwidth(binwidth) {}

   protected:
   target_counts_t target_coverage = {0, 0};

   std::unordered_map<std::string, read_boundary_t> read_cache;
   std::unordered_map<uint64_t, target_counts_t> active_bins;

   uint32_t curstart = 0;
   uint32_t curend = 0;
   int32_t curchr = 0;

   uint32_t binwidth = 0;
   uint64_t binmin = 0;
   uint64_t binmax = 0;
   std::string curchrname = "";
};

}

#endif
