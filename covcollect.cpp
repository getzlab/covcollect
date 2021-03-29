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
   std::string read_name = record.Qname();

   // if we've switched regions, flush the cache and write current counts
   size_t region_idx = reader.GetRegionIdx();
   if(region_idx != cur_region_idx) {
      cur_region = intervals[region_idx];
      cur_region_idx = region_idx;

      // any reads left in the cache will all be singletons
      for(const auto& read : read_cache) {
	 target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end);
	 target_coverage.n_uncorrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end);
      }
      read_cache.clear();

      fprintf(outfile, "%s\t%d\t%d\t%d\t%d\n",
        header.IDtoName(cur_region.chr).c_str(),
        cur_region.pos1,
        cur_region.pos2,
        target_coverage.n_corrected,
        target_coverage.n_uncorrected
      );
      target_coverage = {0, 0};
   }

   // this is the first read in the pair; push to cache
   if(read_cache.find(read_name) == read_cache.end()) {
      read_cache.emplace(
        read_name,
        (read_boundary_t) {
	  (uint32_t) record.Position(),
	  (uint32_t) record.PositionEnd()
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

      // remove from cache
      read_cache.erase(read_name);
   }

   return 1;
}

}

int main(int argc, char** argv) { 
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
   w.load_intervals(151);
   if(!w.set_output_file(args.output_file)) exit(1);

   w.walk_all();

   return 0;
}
