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
	std::string read_name = record.Qname();
	std::cout << "Record chromosome: " << record.ChrID() << "\n";
	int32_t record_chr = record.ChrID();

//	std::map<uint64_t, target_counts_t> ordered_active_bins(active_bins.begin(), active_bins.end());
	for (auto bin = active_bins.begin(); bin != active_bins.end();) {
		// TODO: check chromosome changed
		// TODO: print in order (binmin?
	   if(bin->first + binwidth < record.Position()) {
		   fprintf(outfile, "%s\t%lu\t%lu\t%d\t%d\n",
				 curchrname.c_str(),
				 bin->first,
				 bin->first + binwidth - 1,
				 bin->second.n_corrected,
				 bin->second.n_uncorrected
				);
		   bin = active_bins.erase(bin);
	   }
	   else {
		  if (bin->first <= record.Position()) {
			  binmin = bin->first;
		  }
		  bin++;
	   }
//		std::cout << bin->first << "\t" << bin->second.n_uncorrected << "\n";
	}

//	if (active_bins.size() > 10) {
//		throw runtime_error("here");
//	}

	// TODO: print missing bins

	// TODO: Clear if move to new chromosome

	for(auto read = read_cache.begin(); read != read_cache.end();) {
		if (read->second.end < binmin) {
			read = read_cache.erase(read);
		} else {
			read++;
		}
	}

	// Add bins
	uint32_t start_new_bin = binmax == 0 ? record.Position() : binmax + binwidth;
	for(uint64_t i = start_new_bin; i < record.PositionEnd() + binwidth; i = i + binwidth) {
		active_bins.emplace(i, (target_counts_t){0, 0});
		binmax = i;
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

	   for (auto bin = active_bins.begin(); bin != active_bins.end(); bin++) {
		   bin->second.n_corrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
		   bin->second.n_uncorrected += n_overlap(bin->first, bin->first + binwidth, record.Position(), record.PositionEnd());
	   }

    } else {
    	uint32_t ovlpstart = MAX((int32_t) read_cache[read_name].end, (int32_t) record.Position());

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
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
   w.set_binwidth(50); // TODO: allow this to be specifiable

   if(!w.set_output_file(args.output_file)) exit(1);

   w.walk_all();

   return 0;
}
