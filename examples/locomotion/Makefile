# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := ../../src
ODIR := ../../build

INCS += -I$(RDIR)

SRCS := $(wildcard $(SDIR)/*.cpp)
DEPS := $(patsubst %.cpp, %.h, $(SRCS))
OBJS := $(addprefix $(ODIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

ESRCS := $(wildcard *.cpp)
EOBJS := $(patsubst %.cpp, %.o, $(ESRCS))
EDIR := $(dir $(ESRCS))

TARGETS := $(basename $(ESRCS))

.PHONY : clean


# rules
all: $(TARGETS)

# $(EDIR)/%: $(EDIR)/%.o $(OBJS)
# 	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
main_pipm_unistep: pipm_unistep.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
main_pipm_semisteps: pipm_semisteps.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
test: test.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
walk_pipm_ppm: pipm2ppm.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
walk_pipms: pipms.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
integrated_pipm_ppm: int_pipm2ppm.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
integrated_pipms: int_pipms.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
integrated_ppm_pipm: int_ppm2pipm.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
integrated_pipm_mcm: int_pipm2mcm.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
integrated_mcm_pipm: int_mcm2pipm.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

clean:
	rm -rf *.o *.dSYM $(TARGETS)

cleanall:
	rm -rf *.o *.dSYM $(TARGETS) $(ODIR)/*.o


$(OBJS): | $(ODIR)
$(ODIR):
	mkdir $(ODIR)
$(ODIR)/%.o: $(SDIR)/%.cpp $(SDIR)/%.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

pipm_unistep.o: main_pipm_unistep.cpp pipm2.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
pipm_semisteps.o: main_pipm_semisteps.cpp pipm.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
test.o: test.cpp pipm.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
pipm2ppm.o: walk_pipm_ppm.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
pipms.o: walk_pipms.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
int_pipm2ppm.o: integrated_pipm_ppm.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
int_pipms.o: integrated_pipms.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
int_ppm2pipm.o: integrated_ppm_pipm.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
int_pipm2mcm.o: integrated_pipm_mcm.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
int_mcm2pipm.o: integrated_mcm_pipm.cpp walkstep.hpp $(SDIR)/abstraction.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

# $(EDIR)/%.o: $(EDIR)/%.cpp $(SDIR)/abstraction.hpp
# 	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
