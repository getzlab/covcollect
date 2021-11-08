CC = g++
CFLAGS = -std=c++11
LDFLAGS = -L/mnt/j/proj/misc/20181119_seqlib/src
LDLIBS = -lhts -lcurl -lcrypto -lz -lpthread -llzma -lbz2 -lseqlib
INC = -I/mnt/j/proj/misc/20181119_seqlib/

GAC: covcollect.cpp walker/walker.cpp covcollect.hpp walker/walker.hpp
	$(CC) $(EXTRA) $(CFLAGS) walker/argparse.cpp walker/walker.cpp covcollect.cpp -o covcollect $(LDFLAGS) $(LDLIBS) $(INC)
