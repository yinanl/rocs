# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../../configMake

IDIR := ../../intervals
GDIR := ../../grids
SDIR := ../../src_itvl
RDIR := ../../

INCS += -I$(RDIR)

# define the source files
SRCIVEC = $(IDIR)/interval.cpp $(IDIR)/interval_vector.cpp
SRCSP = $(SDIR)/interval_paver.cpp $(SRCIVEC)
SRCCP = $(SDIR)/csolver.cpp $(SDIR)/problem.cpp $(SRCSP)
SRCGRID = $(GDIR)/grid.cpp

SRC = ipdl.cpp $(SRCGRID) $(SRCCP)


# define the executable file
EXE = ipdl

.PHONY : clean


# rules
$(EXE): $(SRC) ipdl.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $(EXE) $(SRC)

clean:
	rm -rf *.o *.dSYM $(EXE)
