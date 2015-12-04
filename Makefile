CC=gcc
CXX=g++

SHARED_DIR = ../shared
SHARED_LIB = gen_func.cpp read_data.cpp threadpool.c

CFLAGS = -I$(SHARED_DIR)
#DFLAGS = -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
DFLAGS = -O3 -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LIB = -lgsl -lgslcblas -lz -lpthread



all: $(SHARED_LIB) parse_args ngsLD
	$(CXX) $(DFLAGS) *.o $(LIB) -o ngsLD



$(SHARED_LIB):
	$(CXX) $(CFLAGS) $(DFLAGS) -c $(SHARED_DIR)/$@

parse_args:
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

ngsLD:
	$(CXX) $(CFLAGS) $(DFLAGS) -c ngsLD.cpp $(LIB)

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -f *~ *.o ngsLD examples/testLD* examples/test.log
