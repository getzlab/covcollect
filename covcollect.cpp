#include "walker/argparse.hpp"
#include "covcollect.hpp"

#include <SeqLib/GenomicRegionCollection.h>

#define MAX(x, y) (((y) > (x)) ? (y) : (x))
#define MIN(x, y) (((y) < (x)) ? (y) : (x))

using namespace std;

namespace CC {

void cc_walker::load_intervals(uint32_t pad) {
   intervals.ReadBED(interval_list_path, header);
   for(auto& region : intervals) region.Pad(pad);
}

void cc_walker::walk_all() {
   walker::walk(intervals);
}

uint32_t cc_walker::n_overlap(const SeqLib::GenomicRegion& region, uint32_t start, uint32_t end) {
   if(start > region.pos2 || region.pos1 > end) return 0;
   return MIN(region.pos2, end) - MAX(start, region.pos1);
}

bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
   // if we've switched regions, flush the cache and write current counts
   size_t region_idx = reader.GetRegionIdx();
   if(region_idx != cur_region_idx) {
      cur_region = intervals[region_idx];
      cur_region_idx = region_idx;
   }

   // this is the first read in the pair; push to cache
   std::string read_name = record.Qname();
   if(read_cache.find(read_name) == read_cache.end()) {
      read_cache.emplace(
        read_name,
        (read_boundary_t) {
	  (uint32_t) record.Position(),
	  (uint32_t) record.PositionEnd(),
	  (uint32_t) record.GetCigar().NumQueryConsumed()
        }
      );

   // this is the second read in the pair; lookup first read in cache to see
   // if it overlaps
   } else {
      // if the reads overlap, consider them as a single entity
      if(MAX(0, (int32_t) read_cache[read_name].end - (int32_t) record.Position()) > 0) {
	 target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, record.PositionEnd());

      // otherwise, consider their contributions separately
      } else {
	 target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end) + 
           n_overlap(cur_region, record.Position(), record.PositionEnd());
      }

      // save coverage that doesn't account for overlaps for demo purposes
      target_coverage.n_uncorrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end) + 
	n_overlap(cur_region, record.Position(), record.PositionEnd());
   }

   return 1;
}

}

int main(int argc, char** argv) { 
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
   w.load_intervals(151);
   w.walk_all();
   return 0;
}
