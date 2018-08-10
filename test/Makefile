# get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo)
-include ../cflags

# define variables
# better to use relative directory in case there are spaces in absolute directory 
RDIR := ../
SDIR := $(RDIR)src
ODIR := $(RDIR)build

INCS += -I$(RDIR)


SRCS := $(wildcard $(SDIR)/*.cpp)
DEPS := $(patsubst %.cpp, %.h, $(SRCS))
OBJS := $(addprefix $(ODIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

ESRCS := $(wildcard *.cpp)
EDIRS := $(sort $(dir $(ESRCS)))  # use "sort" to remove the duplicate
EINCS := $(wildcard *.hpp)

TARGETS := $(basename $(ESRCS))

.PHONY : clean


# rules
all: $(TARGETS)

testAD: ad.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testGrid: grids.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testIval: ival.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testIvec: ivec.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testPaver: paver.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testPrecision: precision.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testDSolver: dsolver.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^


display:
	@echo $(ODIR)

clean:
	rm -rf *.o *.dSYM $(TARGETS)



$(OBJS): | $(ODIR)
$(ODIR):
	mkdir $(ODIR)
$(ODIR)/%.o: $(SDIR)/%.cpp $(SDIR)/%.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

ad.o: testAD.cpp $(SDIR)/flow_taylor.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

grids.o: testGrid.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

ival.o: testIval.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

ivec.o: testIvec.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

paver.o: testPaver.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

precision.o: testPrecision.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

dsolver.o: testDSolver.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
