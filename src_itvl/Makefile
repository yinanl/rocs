

############################################################################

# the simplest makefile
#
#mcis:mcis.cpp
#	g++ -Wall -g mcis.cpp -o mcis



############################################################################
# REMEMBER TO UPDATE THE ENVIRONMENT VARIABLE DYLD_LIBRARY_PATH in terminal
# export DYLD_LIBRARY_PATH='/Applications/MATLAB_R2016a.app/bin/maci64:/Applications/MATLAB_R2016a.app/sys/os/maci64':$DYLD_LIBRARY_PATH


# ROOTs
BOOSTROOT = /usr/local/Cellar/boost/1.62.0
MATLABROOT = /Applications/MATLAB_R2016a.app
LINALGROOT = /usr/local/Cellar/armadillo/7.500.2

INTERVALROOT = ../intervals
GRIDROOT = ../grids
SOLVERROOT = ../../src_itvl

# for armadillo
#DEBUG = -DARMA_EXTRA_DEBUG
OPT = -O3
#OPT = -O0
#EXTRA_OPT = -fwhole-program
FINAL = -DARMA_NO_DEBUG


# define the c++ compiler to use
CXX = clang++
#CXX = gcc

# define any compile-time flags
CXXFLAGS = -Wall -g $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT) -std=c++11 -stdlib=libc++


# define any directories containing header files other than /usr/include
#
INCLUDES = -I$(BOOSTROOT)/include -I$(LINALGROOT)/include -I$(MATLABROOT)/extern/include #-I/usr/local/include


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LDFLAGS = -L$(BOOSTROOT)/lib -L$(LINALGROOT)/lib -L$(MATLABROOT)/bin/maci64 -rpath $(MATLABROOT)/bin/maci64


# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -larmadillo -lboost_unit_test_framework-mt -lmat -leng -lmx


# define the source files
SRCIVEC = $(INTERVALROOT)/interval.cpp $(INTERVALROOT)/interval_vector.cpp
SRCGRID = $(GRIDROOT)/grid.cpp
SRCSP = $(SOLVERROOT)/interval_paver.cpp $(SRCIVEC)
SRCCP = $(SOLVERROOT)/csolver.cpp $(SOLVERROOT)/problem.cpp $(SRCSP)

SRCPAVER = testPaver.cpp $(SRCSP)
SRCVF = testVF.cpp $(SRCIVEC) $(SRCGRID)


# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#


# define the executable file
EXEPAVER = testPaver
EXEVF = testVF
EXECUTABLE = testVF testPaver


#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#
.PHONY : clean

all: $(EXECUTABLE)

$(EXEPAVER):$(SRCPAVER)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) -o $(EXEPAVER) $(SRCPAVER)

$(EXEVF):$(SRCVF)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) -o $(EXEVF) $(SRCVF)

clean:
	rm -rf $(EXECUTABLE)
