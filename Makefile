CC=gcc
CXX=g++

CFLAGS = -g -Wall -I./shared
#CFLAGS = -O3 -Wall -I./shared
DFLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LIB = -lgsl -lgslcblas -lz -lpthread

SHARED_LIB = gen_func.cpp read_data.cpp threadpool.c



all: $(SHARED_LIB) parse_args ngsLD
	$(CXX) $(DFLAGS) *.o $(LIB) -o ngsLD



$(SHARED_LIB):
	$(CXX) $(CFLAGS) $(DFLAGS) -c shared/$@

parse_args:
	$(CXX) $(CFLAGS) $(DFLAGS) -c parse_args.cpp

ngsLD:
	$(CXX) $(CFLAGS) $(DFLAGS) -c ngsLD.cpp $(LIB)

test:
	@cd examples/; bash test.sh 2> test.log; cd ../

clean:
	@rm -f *~ *.o ngsLD examples/testA_* examples/test.log
