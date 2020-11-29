# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := $(RDIR)src
ODIR := $(RDIR)build

INCS += -I$(RDIR)
CXXFLAGS += -DVERBOSE -DNDEBUG

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

carAbst: carAbst.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carAbst4: carAbst_dba4.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carAbst_ref: carAbst_ref.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carDBA: carDBA.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carDBA4: carDBA_dba4.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carDBA0: carDBA_0.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

carDBAp4: carDBA_p4.o $(OBJS)
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

sim_carAbst: sim_carAbst.o $(OBJS)
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
$(ODIR)/%.o: $(SDIR)/%.c $(SDIR)/%.h
	gcc -Wall -Wextra -o3 -c $< -o $@

%.o: %.cpp $(SDIR)/system.hpp car.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
