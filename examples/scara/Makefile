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

TARGETS := reach dba abst
TESTS := testReachset testInitialize

.PHONY : clean


# rules
all: $(TARGETS)

reach: dbintReach.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

dba: dbintDBA.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

abst: dbintAbst.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

testReachset: testReachset.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

testInitialize: testInitialize.o $(OBJS)
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

%.o: %.cpp $(SDIR)/system.hpp scara.hpp scara_decom.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
