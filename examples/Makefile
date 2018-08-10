# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../cflags


# define variables
SDIR := ../src
ODIR := ../build

INCS += -I../

SRCS := $(wildcard $(SDIR)/*.cpp)
DEPS := $(patsubst %.cpp, %.h, $(SRCS))
OBJS := $(addprefix $(ODIR)/, $(patsubst %.cpp, %.o, $(notdir $(SRCS))))

ESRCS := $(wildcard */*.cpp)
EDIRS := $(sort $(dir $(ESRCS)))  # use "sort" to remove the duplicate
EINCS := $(wildcard */*.hpp)

TARGETS := $(basename $(ESRCS))

.PHONY : clean


$(ODIR)/%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

# rules
examples:
	for dir in $(EDIRS); do \
		$(MAKE) -C $$dir;\
	done

car:
	cd ./car && $(MAKE)
dcdc:
	cd ./dcdc && $(MAKE)
ipdl:
	cd ./ipdl && $(MAKE)
temp:
	cd ./temp && $(MAKE)
vdp:
	cd ./vdp && $(MAKE)

showtargets:
	@echo $(TARGETS)


DELS := $(addsuffix .o, $(TARGETS))


clean:
	rm -rf $(OBJS) $(DELS) $(TARGETS)
