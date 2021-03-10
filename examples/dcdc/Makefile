# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../cflags

RDIR := ../../
SDIR := ../../src
ODIR := ../../build

INCS += -I$(RDIR)
# CXXFLAGS += -DVERBOSE

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

dcdcInv: dcdcInv.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

dcdcRS: dcdcRS.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

dcdcAbstInv: dcdcAbstInv.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

dcdcAbstDBA: dcdcAbstDBA.o $(OBJS) $(ODIR)/buchi.o
	$(CXX) $(CXXFLAGS) $(INCS) $^ -o $@ $(LDFLAGS) $(LIBS)

display:
	@echo $(OBJS)\
	@echo $(EDIR)

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
