# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := ../../src
ODIR := ../../build

INCS += -I$(RDIR)

ifdef MATDIR
SRCS := $(wildcard $(SDIR)/*.cpp)
else
SRCSALL := $(wildcard $(SDIR)/*.cpp)
SRCS := $(filter-out $(SDIR)/matlabio.cpp,$(SRCSALL))
endif

DEPS := $(patsubst %.cpp, %.h, $(SRCS))
OBJS := $(addprefix $(ODIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

ESRCS := $(wildcard *.cpp)
EOBJS := $(patsubst %.cpp, %.o, $(ESRCS))
EDIR := $(dir $(ESRCS))

TARGETS := $(basename $(ESRCS))


.PHONY : clean


# rules
all: $(TARGETS)

roa: roa.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)
reach: reach.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)
real: real.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

clean:
	rm -rf *.o *.dSYM $(TARGETS)

cleanall:
	rm -rf *.o *.dSYM $(TARGETS) $(ODIR)/*.o

display:
	@echo $(OBJS)



$(OBJS): | $(ODIR)
$(ODIR):
	mkdir $(ODIR)
$(ODIR)/%.o: $(SDIR)/%.cpp $(SDIR)/%.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

roa.o: roa.cpp $(SDIR)/system.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

real.o: roa_real.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

reach.o: reach.cpp $(SDIR)/system.hpp $(SDIR)/timer.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
