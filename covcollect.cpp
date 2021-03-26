#include "walker/argparse.hpp"
#include "covcollect.hpp"

#include <SeqLib/GenomicRegionCollection.h>

using namespace std;

namespace CC {

void cc_walker::load_intervals(uint32_t pad) {
   intervals.ReadBED(interval_list_path, header);
   for(auto& region : intervals) region.Pad(pad);
}

void cc_walker::walk_all() {
   walker::walk(intervals);
}

void cc_walker::walk(const SeqLib::GenomicRegion& region) {
   printf("%d:%d-%d\n", region.chr, region.pos1, region.pos2);
   walker::walk(region);
}

bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
   // this is the first read in the pair
   std::string read_name = record.Qname();
   if(read_cache.find(read_name) == read_cache.end()) {
      read_cache.emplace(
        read_name,
        (read_boundary) {
	  (uint32_t) record.Position(),
	  (uint32_t) record.PositionEnd(),
	  (uint32_t) record.GetCigar().NumQueryConsumed()
        }
      );
   } else {
      printf("%d\n", read_cache[read_name].start);
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
