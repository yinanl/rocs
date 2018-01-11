

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

# for armadillo
#DEBUG = -DARMA_EXTRA_DEBUG
OPT = -O2
#EXTRA_OPT = -fwhole-program
#FINAL = -DARMA_NO_DEBUG


# define the c++ compiler to use
CXX = clang++

# define any compile-time flags
CXXFLAGS = -Wall -g $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT) -std=c++11 -stdlib=libc++


# define any directories containing header files other than /usr/include
#
INCLUDES = -I$(BOOSTROOT)/include -I$(LINALGROOT)/include -I$(MATLABROOT)/extern/include #-I/usr/local/include


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LDFLAGS = -L$(BOOSTROOT)/lib -L$(LINALGROOT)/lib -L$(MATLABROOT)/bin/maci64 -rpath /Applications/MATLAB_R2016a.app/bin/maci64


# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -larmadillo -lboost_unit_test_framework-mt -lmat -leng -lmx


# define the source files
SRCIVEC = interval.cpp interval_vector.cpp


SRC01 = testIval.cpp interval.cpp
SRC02 = testIvec.cpp $(SRCIVEC)


# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#


# define the executable file
EXE01 = testIval
EXE02 = testIvec



#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#



all: prog01 prog02 prog03 prog04 prog05 prog06

prog01: $(SRC01)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) -o $(EXE01) $(SRC01)

prog02: $(SRC02)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) -o $(EXE02) $(SRC02)


clean:
	rm -rf $(EXECUTABLE)
