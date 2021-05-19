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
   this->pad = pad;
}

void cc_walker::set_binwidth(uint32_t binwidth) {
	this->binwidth = binwidth;
}

void cc_walker::walk_all() {
   walker::walk();
}

uint32_t cc_walker::n_overlap(const uint32_t binstart, uint32_t binend, uint32_t start, uint32_t end) {
   if(start > binstart - this->pad || binstart + this->pad > end) return 0;
   return MIN(binend - this->pad, end) - MAX(start, binstart + this->pad);
}


bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
	std::cout << "\nNext read start pos:" << record.Position() << "\n";
	std::string read_name = record.Qname();

	// Write and delete
	std::cout << "\nwrite and delete: ";

	for (auto bin = active_bins.begin(); bin != active_bins.end();) {
	   if(bin->first + binwidth < record.Position()) {
		   fprintf(outfile, "%s\t%lu\t%lu\t%d\t%d\n",
				 curchrname.c_str(),
				 bin->first,
				 bin->first + binwidth - 1,
				 bin->second.n_corrected,
				 bin->second.n_uncorrected
				);
		   std::cout << bin->first << ", ";
		   bin = active_bins.erase(bin);
	   }
	   else {
		  if (bin->first <= record.Position()) {
			  binmin = bin->first;
		  }
		  bin++;
	   }
	}

	for(const auto read = read_cache.begin(); read != read_cache.end();) {
		if (read.second.end < binmin) {
			read_cache.erase(read);
		}
		read++;
	}

	// TODO: check if missing regions

	// Add bins
	std::cout << "\nAdd bins: ";
	uint32_t start_new_bin = binmax == 0 ? record.Position() : binmax + binwidth;
	for(uint64_t i = start_new_bin; i < record.PositionEnd() + binwidth; i = i + binwidth) {
		active_bins.emplace(i, (target_counts_t){0, 0});
		std::cout << i << ", ";
		binmax = i;
	}
	std::cout << "new min: " << binmin << "\n";
	std::cout << "new max: " << binmax << "\n";

	std::cout << "active_bins.size() is " << active_bins.size() << '\n';

    // this is the first read in the pair; push to cache
    if(read_cache.find(read_name) == read_cache.end()) {
	   read_cache.emplace(
		 read_name,
		 (read_boundary_t) {
	   (uint32_t) record.Position(),
	   (uint32_t) record.PositionEnd()
		 }
	   );

	   for (auto bin = active_bins.begin(); bin != active_bins.end();) {
		   bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
	   }

    } else {
    	std::cout << "Pair!\n";
    	uint32_t ovlpstart = MAX((int32_t) read_cache[read_name].end, (int32_t) record.Position());

		std::cout << "read_cache[read_name].start: " << read_cache[read_name].start << "\n";
		std::cout << "read_cache[read_name].end: " << read_cache[read_name].end << "\n";
		std::cout << "record.Position(): " << record.Position() << "\n";
		std::cout << "record.PositionEnd(): " << record.PositionEnd() << "\n";
		std::cout << "Looping through bins\n";
		for (auto bin = active_bins.begin(); bin != active_bins.end();) {
			std::cout << "bin->first: " << bin->first << "\n";
			bin->second.n_uncorrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
			bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, ovlpstart, record.PositionEnd());
			bin++;
		}

		// remove from cache
		read_cache.erase(read_name);
    }

	return 1;

}


//bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
//   std::string read_name = record.Qname();
//
//   // if we've switched regions, flush the cache and write current counts
//   size_t region_idx = reader.GetRegionIdx();
//   if(region_idx != cur_region_idx) {
//      // any reads left in the cache will all be singletons
//      for(const auto& read : read_cache) {
//	 target_coverage.n_corrected += n_overlap(cur_region, read.second.start, read.second.end);
//	 target_coverage.n_uncorrected += n_overlap(cur_region, read.second.start, read.second.end);
//      }
//      read_cache.clear();
//
//      fprintf(outfile, "%s\t%d\t%d\t%d\t%d\n",
//        header.IDtoName(cur_region.chr).c_str(),
//        cur_region.pos1 + this->pad,
//        cur_region.pos2 - this->pad,
//        target_coverage.n_corrected,
//        target_coverage.n_uncorrected
//      );
//      target_coverage = {0, 0};
//
//      // we may have skipped over multiple empty regions
//      for(size_t r = cur_region_idx + 1; r < region_idx; r++) {
//	 SeqLib::GenomicRegion gr = intervals[r];
//	 fprintf(outfile, "%s\t%d\t%d\t%d\t%d\n",
//	   header.IDtoName(gr.chr).c_str(),
//	   gr.pos1 + this->pad,
//	   gr.pos2 - this->pad,
//	   0,
//	   0
//	 );
//      }
//
//      // switch to next region
//      cur_region = intervals[region_idx];
//      cur_region_idx = region_idx;
//   }
//
//   // this is the first read in the pair; push to cache
//   if(read_cache.find(read_name) == read_cache.end()) {
//      read_cache.emplace(
//        read_name,
//        (read_boundary_t) {
//	  (uint32_t) record.Position(),
//	  (uint32_t) record.PositionEnd()
//        }
//      );
//
//   // this is the second read in the pair; lookup first read in cache to see
//   // if it overlaps
//   } else {
//      // if the reads overlap, consider them as a single entity
//      if(MAX(0, (int32_t) read_cache[read_name].end - (int32_t) record.Position()) > 0) {
//	 target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, record.PositionEnd());
//
//      // otherwise, consider their contributions separately
//      } else {
//	 target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end) +
//           n_overlap(cur_region, record.Position(), record.PositionEnd());
//      }
//
//      // save coverage that doesn't account for overlaps for demo purposes
//      target_coverage.n_uncorrected += n_overlap(cur_region, read_cache[read_name].start, read_cache[read_name].end) +
//	n_overlap(cur_region, record.Position(), record.PositionEnd());
//
//      // remove from cache
//      read_cache.erase(read_name);
//   }
//
//   return 1;
//}

}

int main(int argc, char** argv) { 
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
   w.set_binwidth(50); // TODO: allow this to be specifiable

   if(!w.set_output_file(args.output_file)) exit(1);

   w.walk_all();

   return 0;
}
