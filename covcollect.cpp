#include "walker/argparse.hpp"
#include "covcollect.hpp"

#include <SeqLib/GenomicRegionCollection.h>

using namespace std;

namespace CC {

void cc_walker::load_intervals() {
   intervals.ReadBED(interval_list_path, header);
}

void cc_walker::walk() {
   walker::walk(intervals);
}

bool cc_walker::walk_apply(const SeqLib::BamRecord& record) {
   return 0;
}

}

int main(int argc, char** argv) { 
   walker::basic_arg_t args = {};
   if(!walker::basic_argparse(argc, argv, &args)) exit(1);

   CC::cc_walker w = CC::cc_walker(args.bam_in, args.input_file);
   w.load_intervals();
   w.walk();
   return 0;
}
