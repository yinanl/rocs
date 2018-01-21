# 1) get directories of external packages: 
# BSTDIR (boost), MATDIR (matlab), LADIR (armadillo).
# 2) get compiling and linking flags:
# CXX, CXXFLAGS, INCS, LIBS, LIBTEST, LDFLAGS.
-include ../configMake


# define variables
IDIR := ../intervals
GDIR := ../grids
SDIR := ../src_itvl
RDIR := ../


INCS += -I$(RDIR)


ESRCS := $(wildcard */*.cpp)
EDIRS := $(sort $(dir $(ESRCS)))
EINCS := $(wildcard */*.hpp)

TARGETS := $(basename $(ESRCS))


SRCS := $(wildcard $(IDIR)/*.cpp) \
	$(wildcard $(GDIR)/*.cpp) \
	$(wildcard $(SDIR)/*.cpp)

OBJS := $(patsubst %.cpp, %.o, $(SRCS))
DEPS := $(patsubst %.cpp, %.h, $(SRCS))


.PHONY : clean

define rules =
 $(1): $$($(1).cpp) $$($(OBJS))
	$$($(CXX)) $$($(CXXFLAGS)) $$($(INCS)) $$($(LDFLAGS)) $$($(LIBS)) -o $$($(1)) $$($(1).cpp) $$($(OBJS))
endef

$(foreach case, $(TARGETS), $(eval $(call rules, $(case))))

# rules
all: $(TARGETS)

#$(TARGETS):
#	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) -o $@ $^

examples:
	for dir in $(EDIRS); do \
		$(MAKE) -C $$dir;\
	done

testExamples: testVF.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCS) $(LDFLAGS) $(LIBS) $(LIBTEST) -o $@ $^

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@


disp:
	@echo $(TARGETS)


DELS := $(addsuffix .dSYM, $(TARGETS) testExamples)

clean:
	rm -rf $(OBJS) $(DELS) $(TARGETS) testExamples
