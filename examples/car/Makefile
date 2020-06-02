# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := $(RDIR)src
ODIR := $(RDIR)build

INCS += -I$(RDIR)
CXXFLAGS += -DVERBOSE

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

carTraj: carReachTraj.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carAbst: carAbst.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carReach1: carReach1.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carReach2: carReach2.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carBuchi: carBuchi.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carParking: carParking.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carFullmap: carFullmap.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

testConstraints: testConstraints.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

showtargets:
	@echo $(TARGETS)

clean:
	rm -rf *.o *.dSYM $(TARGETS)

cleanall:
	rm -rf *.o *.dSYM $(TARGETS) $(ODIR)/*.o




$(OBJS): | $(ODIR)
$(ODIR):
	mkdir $(ODIR)
$(ODIR)/%.o: $(SDIR)/%.cpp $(SDIR)/%.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

%.o: %.cpp $(SDIR)/system.hpp car.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
