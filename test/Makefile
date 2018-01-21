# get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo)
-include ../configMake

IDIR := ../intervals
GDIR := ../grids
SDIR := ../src_itvl
RDIR := ../

INCS += -I$(RDIR)


TARGETS = testInterval testIvec testGrid testPaver


ISRCS = $(wildcard $(IDIR)/*.cpp)
GSRCS = $(wildcard $(GDIR)/*.cpp)
SSRCS = $(wildcard $(SDIR)/*.cpp)

IOBJS = $(patsubst %.cpp, %.o, $(ISRCS))
GOBJS = $(patsubst %.cpp, %.o, $(ISRCS) $(GSRCS))
SOBJS = $(patsubst %.cpp, %.o, $(ISRCS) $(GSRCS) $(SSRCS))
DEPS = $(patsubst %.cpp, %.h, $(SRCS))


.PHONY : clean

all: $(TARGETS)

testInterval: testIval.cpp $(IOBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testIvec: testIvec.cpp $(IOBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testGrid: testGrid.cpp $(GOBJS) $(IOBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testPaver: testPaver.cpp $(GOBJS) $(IOBJS) $(SOBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^


%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@


DELS := $(addsuffix .dSYM, $(TARGETS))

clean:
	rm -rf $(GOBJS) $(IOBJS) $(SOBJS) $(DELS) $(TARGETS)
