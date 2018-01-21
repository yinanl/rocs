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

SRC01 = carReach1.cpp $(SRCGRID) $(SRCCP)
SRC02 = carReach2.cpp $(SRCGRID) $(SRCCP)
SRCBUCHI = carBuchi.cpp $(SRCGRID) $(SRCCP)


# define the executable file
EXE01 = carReach1
EXE02 = carReach2
EXEBUCHI = carBuchi
EXECUTABLE = carReach1 carReach2 carBuchi

.PHONY : clean


# rules
all: $(EXECUTABLE)

$(EXE01): $(SRC01) car.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $(EXE01) $(SRC01)

$(EXE02): $(SRC02) car.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $(EXE02) $(SRC02)

$(EXEBUCHI): $(SRCBUCHI) car.hpp
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $(EXEBUCHI) $(SRCBUCHI)

clean:
	rm -rf *.o *.dSYM $(EXECUTABLE)
