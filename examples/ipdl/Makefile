# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := ../../src
ODIR := ../../build

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

ipdl: ipdl.o $(OBJS)
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

%.o: %.cpp $(SDIR)/system.hpp
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

# $(EDIR)/%.o: $(EDIR)/%.cpp $(SDIR)/system.hpp
#	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@
