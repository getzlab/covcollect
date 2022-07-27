#include "walker/argparse.hpp"
#include "covcollect.hpp"

#include <SeqLib/GenomicRegionCollection.h>

#include <math.h>

#define MAX(x, y) (((y) > (x)) ? (y) : (x))
#define MIN(x, y) (((y) < (x)) ? (y) : (x))

using namespace std;

namespace CC {

// override base filter to noop; read filtering occurs within walk_apply, since
// we are interested in tabulating the number of filtered reads
bool cc_walker::filter_read(const SeqLib::BamRecord& record) {
   return false;
}

void cc_walker::load_intervals(uint32_t pad, string chr, uint32_t start, uint32_t end) {
   SeqLib::GenomicRegionCollection<> all_intervals;
   all_intervals.ReadBED(interval_list_path, header);
   // trim intervals outside of (start, end)
   uint32_t idx = 0;
   int32_t chr_idx = chr != "-1" ? header.Name2ID(chr) : -1;
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
   if(intervals.IsEmpty()) return;

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
      // any reads left in the cache will all be in the previous region (which is still called "cur_region")
      for(const auto& read : read_cache) {
	 /*if(intervals[read.second.reg_idx].pos1 != cur_region.pos1) {
	    printf("%d %d\n", region_idx, read.second.reg_idx);
	    exit(1);
	 }*/
	 target_coverage.n_corrected += n_overlap(cur_region, read.second.start, read.second.end);

	 // we do not use these reads for calculating fragment length.
      }
      read_cache.clear();

      fprintf(outfile, "%s\t%d\t%d\t%d\t%0.0f\t%0.0f\t%d\n",
        header.IDtoName(cur_region.chr).c_str(),
        cur_region.pos1 + this->pad,
        cur_region.pos2 - this->pad,
        target_coverage.n_corrected,
        target_coverage.mean_fraglen,
        sqrt(target_coverage.var_fraglen/(target_coverage.n_frags)),
        target_coverage.n_frags
      );
      target_coverage = {0, 0, 0, 0};

      // we may have skipped over multiple empty regions
      for(size_t r = cur_region_idx + 1; r < region_idx; r++) {
	 SeqLib::GenomicRegion gr = intervals[r];
	 fprintf(outfile, "%s\t%d\t%d\t%d\t%0.0f\t%0.0f\t%d\n",
	   header.IDtoName(gr.chr).c_str(),
	   gr.pos1 + this->pad,
	   gr.pos2 - this->pad,
	   0,
	   0.0,
	   0.0,
	   0
	 );
      }

      // switch to next region
      cur_region = intervals[region_idx];
      cur_region_idx = region_idx;
   }

   // apply read filtering; record whether read was filtered
   // TODO


   // this is the first read in the pair; push to cache
   if(read_cache.find(read_name) == read_cache.end()) {
      read_cache.emplace(
        read_name,
        (read_boundary_t) {
	  (uint32_t) record.Position(),
	  (uint32_t) record.PositionEnd()
        }
      );

   // this is the second read in the pair; write coverage to target
   } else {
      target_coverage.n_corrected += n_overlap(cur_region, read_cache[read_name].start, record.PositionEnd());

      // update mean fragment length for this region
      uint32_t fraglen = record.PositionEnd() - read_cache[read_name].start;
      if(target_coverage.n_frags == 0) {
          target_coverage.mean_fraglen = fraglen;
      } else {
          float fld = fraglen - target_coverage.mean_fraglen;
          target_coverage.mean_fraglen += fld/(target_coverage.n_frags + 1);
          target_coverage.var_fraglen += fld*(fraglen - target_coverage.mean_fraglen);
      }
      target_coverage.n_frags++;

      // remove from cache
      read_cache.erase(read_name);
   }

   return 1;
}

void cc_bin_walker::walk_all(string chr, uint32_t start, uint32_t end) {
   int32_t chr_idx = chr != "-1" ? header.Name2ID(chr) : -1;
   if(chr_idx == -1) walker::walk();
   else {
      binmax = start - (start % binwidth);
      curchr = chr_idx;
      walker::walk(SeqLib::GenomicRegion(chr_idx, start < 1000 ? 0 : (start - 1000), end));
      // NOTE: we add a 1kb buffer to the start to capture any readpairs upstream of start that 
      //       nonetheless might overlap it
   }
}

uint32_t cc_bin_walker::n_overlap(const uint32_t binstart, uint32_t binend, uint32_t start, uint32_t end) {
   if(start > binstart + this->binwidth || binstart > end) return 0;
   return MIN(binend, end) - MAX(start, binstart);
}


bool cc_bin_walker::walk_apply(const SeqLib::BamRecord &record) {
   std::string read_name = record.Qname();
   int32_t record_chr = record.ChrID();

   // we've switched chromosomes; flush cache
   if (record_chr != curchr) {
       std::map<uint64_t, target_counts_t> ordered_active_bins(active_bins.begin(), active_bins.end());
       for (auto bin = ordered_active_bins.begin(); bin != ordered_active_bins.end(); ++bin) {
           fprintf(outfile, "%s\t%lu\t%lu\t%d\t%0.0f\t%0.0f\t%d\n",
             header.IDtoName(curchr).c_str(),
             bin->first,
             bin->first + binwidth - 1,
             bin->second.n_corrected,
             bin->second.mean_fraglen,
             sqrt(bin->second.var_fraglen/(bin->second.n_frags)),
             bin->second.n_frags
           );
       }
       active_bins.clear();
       read_cache.clear();
       curchr = record_chr;
       binmax = 0;
   } else {
       std::map<uint64_t, target_counts_t> ordered_active_bins(
               active_bins.begin(), active_bins.end());
       for (auto bin = ordered_active_bins.begin();
               bin->first + binwidth < record.Position() && bin != ordered_active_bins.end();
                   ++bin) {

           fprintf(
             outfile, "%s\t%lu\t%lu\t%d\t%0.0f\t%0.0f\t%d\t%d\t%d\n",
             header.IDtoName(curchr).c_str(),
             bin->first,
             bin->first + binwidth - 1,
             bin->second.n_corrected,
             bin->second.mean_fraglen,
             sqrt(bin->second.var_fraglen/(bin->second.n_frags)),
             bin->second.n_frags,
             bin->second.n_tot_reads,
             bin->second.n_fail_reads
           );
           active_bins.erase(bin->first);
       }
   }

   // Print gaps
   for (uint64_t i = binmax; i + binwidth < record.Position(); i = i + binwidth) {
       fprintf(
         outfile, "%s\t%lu\t%lu\t%d\t%0.0f\t%0.0f\t%d\t%d\t%d\n",
         header.IDtoName(curchr).c_str(),
         i,
         i + binwidth - 1,
         0,
         0.0,
         0.0,
         0,
         0,
         0
       );
   }

   binmax = MAX(binmax, (record.Position() / binwidth) * binwidth);

   // Add bins
   for (uint64_t i = binmax; i <= record.PositionEnd(); i = i + binwidth) {
       active_bins.emplace(i, (target_counts_t ) { 0, 0, 0, 0, 0, 0 });
   }

   binmax = MAX(binmax, ((record.PositionEnd() / binwidth) + 1) * binwidth);

   // apply read filtering; record whether read was filtered
   bool fail = walker::filter_read(record) || !record.ProperPair();
   for (auto bin = active_bins.begin(); bin != active_bins.end(); bin++) {
       if(record.Position() >= bin->first && record.Position() < bin->first + binwidth) {
           bin->second.n_fail_reads += (int) fail;
           bin->second.n_tot_reads++;
           break;
       }
   }
   if(fail) return 1;

   // this is the first read in the pair; push to cache
   if (read_cache.find(read_name) == read_cache.end()) {
       read_cache.emplace(
         read_name,
         (read_boundary_t ) {
             (uint32_t) record.Position(),
             (uint32_t) record.PositionEnd()
         }
       );

   // this is the second read in the pair; write coverage to bins
   } else {
       uint32_t midpoint = (read_cache[read_name].start + record.PositionEnd())/2;

       for (auto bin = active_bins.begin(); bin != active_bins.end(); bin++) {
           bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, read_cache[read_name].start, record.PositionEnd());

           // update mean fragment length for this bin, if read's midpoint is in it
           if(midpoint >= bin->first && midpoint < bin->first + binwidth) {
               uint32_t fraglen = record.PositionEnd() - read_cache[read_name].start;
               if(bin->second.n_frags == 0) {
                   bin->second.mean_fraglen = fraglen;
               } else {
                   float fld = fraglen - bin->second.mean_fraglen;
                   bin->second.mean_fraglen += fld/(bin->second.n_frags + 1);
                   bin->second.var_fraglen += fld*(fraglen - bin->second.mean_fraglen);
               }
               bin->second.n_frags++;
           }
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
     .chr = "-1",
     .start = 0,
     .end = 0xFFFFFFFF
   };

   char arg;
   optind = 1;
   while((arg = getopt(argc, argv, "c:s:e:")) != -1) {
      // TODO: parse -L chr:start-end syntax
      switch(arg) {
	 case 'c' :
	    extra_args.chr = string(optarg);
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
      CC::cc_bin_walker w = CC::cc_bin_walker(args.bam_in, stoi(args.input_file));   
      if(!w.set_output_file(args.output_file)) exit(1);
      w.walk_all(extra_args.chr, extra_args.start, extra_args.end);
   } else if (is_file) {
      CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
      w.load_intervals(151, extra_args.chr, extra_args.start, extra_args.end);
      if(!w.set_output_file(args.output_file)) exit(1);
      w.walk_all();
   } else {
      throw runtime_error("Invalid input to -i. Provide a positive interer for bin mode, or a valid file path to an interval file for interval mode.");
   }

   return 0;
}
