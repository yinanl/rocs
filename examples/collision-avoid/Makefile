# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := ../../src
ODIR := ../../build

INCS += -I$(RDIR)
CXXFLAGS += -DVERBOSE #-DLOGGING

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

safeset: safeset.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
safesetAbst: safesetAbst.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
planning: planning.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
simulation: simulation.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

testReachset: testReachset.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
testWinGraph: testWinGraph.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
ref: carAbst_ref.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^
sim_ref: sim_carAbst.o $(OBJS)
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
$(ODIR)/%.o: $(SDIR)/%.c $(SDIR)/%.h
	gcc -Wall -Wextra -o3 -c $< -o $@

%.o: %.cpp $(SDIR)/system.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
