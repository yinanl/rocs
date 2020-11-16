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

testAD: testAD.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testGrid: testGrid.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testIval: testIval.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testIvec: testIvec.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testPaver: testPaver.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testPrecision: testPrecision.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testDSolver: testDsolver.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testHDF5: testHDF5.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

testAbst: testAbst.o $(OBJS)
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

%.o: %.cpp $(SDIR)/flow_taylor.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
