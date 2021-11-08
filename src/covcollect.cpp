#include "walker/argparse.hpp"
#include "covcollect.hpp"

#include <SeqLib/GenomicRegionCollection.h>

#define MAX(x, y) (((y) > (x)) ? (y) : (x))
#define MIN(x, y) (((y) < (x)) ? (y) : (x))

using namespace std;

namespace CC {

void cc_walker::load_intervals(uint32_t pad, int32_t chr_idx, uint32_t start, uint32_t end) {
   SeqLib::GenomicRegionCollection<> all_intervals;
   all_intervals.ReadBED(interval_list_path, header);
   // trim intervals outside of (start, end)
   uint32_t idx = 0;
   for(auto& region : all_intervals) {
      // interval does not lie within trimming region
      if(chr_idx != -1 && (region.chr != chr_idx || region.pos1 < start || region.pos2 > end)) {
	 continue;
      }
      region.Pad(pad);
      intervals.add(region);
   }
   this->pad = pad;
}

void cc_walker::walk_all() {
   cur_region = intervals[0];
   walker::walk(intervals);
}

uint32_t cc_walker::n_overlap(const SeqLib::GenomicRegion& region, uint32_t start, uint32_t end) {
   if(start > region.pos2 - this->pad || region.pos1 + this->pad > end) return 0;
   return MIN(region.pos2 - this->pad, end) - MAX(start, region.pos1 + this->pad);
}

bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
   std::string read_name = record.Qname();

   // if we've switched regions, flush the cache and write current counts
   size_t region_idx = reader.GetRegionIdx();
   if(region_idx != cur_region_idx) {
      // any reads left in the cache will all be singletons
      for(const auto& read : read_cache) {
	 /*if(intervals[read.second.reg_idx].pos1 != cur_region.pos1) {
	    printf("%d %d\n", region_idx, read.second.reg_idx);
	    exit(1);
	 }*/
	 target_coverage.n_corrected += n_overlap(cur_region, read.second.start, read.second.end);
	 target_coverage.n_uncorrected += n_overlap(cur_region, read.second.start, read.second.end);
      }
      read_cache.clear();

      fprintf(outfile, "%s\t%d\t%d\t%d\t%d\n",
        header.IDtoName(cur_region.chr).c_str(),
        cur_region.pos1 + this->pad,
        cur_region.pos2 - this->pad,
        target_coverage.n_corrected,
        target_coverage.n_uncorrected
      );
      target_coverage = {0, 0};

      // we may have skipped over multiple empty regions
      for(size_t r = cur_region_idx + 1; r < region_idx; r++) {
	 SeqLib::GenomicRegion gr = intervals[r];
	 fprintf(outfile, "%s\t%d\t%d\t%d\t%d\n",
	   header.IDtoName(gr.chr).c_str(),
	   gr.pos1 + this->pad,
	   gr.pos2 - this->pad,
	   0,
	   0
	 );
      }

      // switch to next region
      cur_region = intervals[region_idx];
      cur_region_idx = region_idx;
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

void cc_bin_walker::walk_all() {
   walker::walk();
}

uint32_t cc_bin_walker::n_overlap(const uint32_t binstart, uint32_t binend, uint32_t start, uint32_t end) {
   if(start > binstart + this->binwidth || binstart > end) return 0;
   return MIN(binend, end) - MAX(start, binstart);
}


bool cc_bin_walker::walk_apply(const SeqLib::BamRecord &record) {
    std::string read_name = record.Qname();
    int32_t record_chr = record.ChrID();

    if (record_chr != curchr) {
        // print and erase everything
        std::map<uint64_t, target_counts_t> ordered_active_bins(active_bins.begin(), active_bins.end());
        for (auto bin = ordered_active_bins.begin(); bin != ordered_active_bins.end(); ++bin) {
            fprintf(outfile, "%s\t%lu\t%lu\t%d\t%d\n",
              header.IDtoName(curchr).c_str(),
              bin->first,
              bin->first + binwidth - 1,
              bin->second.n_corrected,
              bin->second.n_uncorrected
            );
        }
        active_bins.clear();
        read_cache.clear();
        curchr = record_chr;
        binmin = 0;
        binmax = 0;
    } else {
        std::map<uint64_t, target_counts_t> ordered_active_bins(
                active_bins.begin(), active_bins.end());
        for (auto bin = ordered_active_bins.begin();
                bin->first + binwidth < record.Position() && bin != ordered_active_bins.end();
                    ++bin) {

            fprintf(
              outfile, "%s\t%lu\t%lu\t%d\t%d\n",
              header.IDtoName(curchr).c_str(),
              bin->first,
              bin->first + binwidth - 1,
              bin->second.n_corrected,
              bin->second.n_uncorrected
            );
            active_bins.erase(bin->first);
        }

        binmin = record.Position() - (record.Position() % binwidth);
    }

    // Print gaps
    for (uint64_t i = binmax; i + binwidth < record.Position(); i = i + binwidth) {
        fprintf(
	  outfile, "%s\t%lu\t%lu\t%d\t%d\n",
          header.IDtoName(curchr).c_str(),
          i,
          i + binwidth - 1,
          0,
          0
        );
    }

    binmax = MAX(binmax, (record.Position() / binwidth) * binwidth);

    // Add bins
    for (uint64_t i = binmax; i < record.PositionEnd(); i = i + binwidth) {
        active_bins.emplace(i, (target_counts_t ) { 0, 0 });
    }

    binmax = MAX(binmax, ((record.PositionEnd() / binwidth) + 1) * binwidth);

    // this is the first read in the pair; push to cache
    if (read_cache.find(read_name) == read_cache.end()) {
        read_cache.emplace(
          read_name,
          (read_boundary_t ) {
              (uint32_t) record.Position(),
              (uint32_t) record.PositionEnd()
          }
        );

        for (auto bin = active_bins.begin(); bin != active_bins.end(); bin++) {
            bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
            bin->second.n_uncorrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
        }

    } else {
        uint32_t ovlpstart = MAX((int32_t ) read_cache[read_name].end, (int32_t ) record.Position());

        for (auto bin = active_bins.begin(); bin != active_bins.end(); bin++) {
            bin->second.n_uncorrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
            bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, ovlpstart, record.PositionEnd());
        }

        // remove from cache
        read_cache.erase(read_name);
    }

    return 1;

}

}

int main(int argc, char** argv) { 
   // parse common args
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   // parse args for trimming region
   CC::extra_args_t extra_args = {
     .chr_idx = -1,
     .start = 0,
     .end = 0xFFFFFFFF
   };

   char arg;
   optind = 1;
   while((arg = getopt(argc, argv, "c:s:e:")) != -1) {
      // TODO: parse -L chr:start-end syntax
      switch(arg) {
	 case 'c' :
	    extra_args.chr_idx = atoi(optarg);
	    break;
	 case 's' : // position to start in BAM
	    extra_args.start = atoi(optarg);
	    break;
	 case 'e' : // position to end in BAM
	    extra_args.end = atoi(optarg);
	    break;
      }
   }

   bool use_bins = true;

   for (int i=0; i < args.input_file.length(); i++) {
      if (isdigit(args.input_file[i]) == false) {
         use_bins = false;
         break;
      }
   }

   std::ifstream is_file(args.input_file);

   if (use_bins) {
      CC::cc_bin_walker  w = CC::cc_bin_walker(args.bam_in, stoi(args.input_file));   
      if(!w.set_output_file(args.output_file)) exit(1);
      w.walk_all();
   } else if (is_file) {
      CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
      w.load_intervals(151, extra_args.chr_idx, extra_args.start, extra_args.end);
      if(!w.set_output_file(args.output_file)) exit(1);
      w.walk_all();
   } else {
      throw runtime_error("Invalid input to -i. Provide a positive interer for bin mode, or a valid file path to an interval file for interval mode.");
   }

   return 0;
}
