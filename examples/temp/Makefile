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

SRC = temp.cpp $(SRCCP)


# define the executable file
EXECUTABLE = temp

.PHONY : clean


# rules

all: $(EXECUTABLE)

$(EXECUTABLE):$(SRC) temp.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $(EXECUTABLE) $(SRC)

clean:
	rm -rf *.o *.dSYM $(EXECUTABLE)
