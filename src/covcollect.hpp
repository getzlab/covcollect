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
   float mean_fraglen;
   float var_fraglen;
   uint32_t n_frags; // total number of analyzed fragments in this bin; used for computing iterative mean/variance
   uint32_t n_tot_reads; // total number of reads in this bin; used for filtering bins with bad reads
   uint32_t n_fail_reads; // number of reads failing filters
} target_counts_t;

typedef struct ea {
   std::string chr;
   uint32_t start;
   uint32_t end;
} extra_args_t;

class cc_walker : public walker::walker {
   public:
   bool filter_read(const SeqLib::BamRecord& record) override;
   bool walk_apply(const SeqLib::BamRecord& record);
   void load_intervals(uint32_t pad, std::string chr, uint32_t start, uint32_t end);
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

class cc_bin_walker : public cc_walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);
   void walk_all(std::string chr, uint32_t start, uint32_t end);
   uint32_t n_overlap(const uint32_t binstart, uint32_t binend, uint32_t start, uint32_t end);

   cc_bin_walker(const std::string& bam_in, const uint32_t binwidth) : cc_walker(bam_in, ""), binwidth(binwidth) {}

   protected:
   target_counts_t target_coverage = {0, 0};

   std::unordered_map<std::string, read_boundary_t> read_cache;
   std::unordered_map<uint64_t, target_counts_t> active_bins;

   uint32_t curstart = 0;
   uint32_t curend = 0;
   int32_t curchr = 0;

   uint32_t binwidth = 0;
   uint64_t binmax = 0;
   std::string curchrname = "";
};

class cc_bin_walker_se : public cc_bin_walker {
   public:
   bool walk_apply(const SeqLib::BamRecord& record);

   cc_bin_walker_se(const std::string& bam_in, const uint32_t binwidth) : cc_bin_walker(bam_in, binwidth) {}
};

}

#endif
