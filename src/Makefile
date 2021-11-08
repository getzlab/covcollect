CC = g++
CFLAGS = -std=c++11
LDFLAGS = -L/home/cchu/SeqLib/src
LDLIBS = -lhts -lcurl -lcrypto -lz -lpthread -llzma -lbz2 -lseqlib
INC = -I/home/cchu/SeqLib/

GAC: covcollect.cpp walker/walker.cpp covcollect.hpp walker/walker.hpp
	$(CC) $(EXTRA) $(CFLAGS) walker/argparse.cpp walker/walker.cpp covcollect.cpp -o covcollect $(LDFLAGS) $(LDLIBS) $(INC)
