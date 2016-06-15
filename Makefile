PROGNAME:=hdbdd

buildtype:=release
APP:=app/
SRC:=src/
CSRCS:=$(SRC)bitvector.c $(SRC)bp.c $(SRC)bpcore.c $(SRC)darray.c $(SRC)mmap.c
CSRCS_ORG:=$(CSRCS) $(SRC)bddc.c
CPPSRCS:=$(SRC)hdbdd.cpp $(SRC)bitio.cpp
CPPSRCS_ORG:=
CCSRCS:=$(SRC)BDD.cc $(SRC)ZBDD.cc $(SRC)SOP.cc $(SRC)CtoI.cc
CC:=gcc
CXX:=g++
CFLAGS:=-I./src/ -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2
CXXFLAGS+=-I./src/ -std=c++0x -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64
DEFS:=-DB_64
LDFLAGS+=
LIBS+=
ifeq ($(buildtype),release)
	CFLAGS+=-O3
	CXXFLAGS+=-O3
else ifeq ($(buildtype),debug)
	CFLAGS+=-O0 -g
	CXXFLAGS+=-O0 -g
else
	$(error buildtype must be release, debug, profile or coverage)
endif

OUTDIR:=./build/$(buildtype)/

PROG:=$(OUTDIR)$(PROGNAME)
COBJS:=$(CSRCS:%.c=$(OUTDIR)%.o)
COBJS_ORG:=$(CSRCS_ORG:%.c=$(OUTDIR)%.o)
CPPOBJS:=$(CPPSRCS:%.cpp=$(OUTDIR)%.o)
CPPOBJS_ORG:=$(CPPSRCS_ORG:%.cpp=$(OUTDIR)%.o)
CCOBJS:=$(CCSRCS:%.cc=$(OUTDIR)%.o)
CCOBJS_ORG:=$(CCSRCS:%.cc=$(OUTDIR)%_org.o)
CDEPS:=$(CSRCS:%.c=$(OUTDIR)%.d)
CDEPS_ORG:=$(CSRCS_ORG:%.c=$(OUTDIR)%.d)
CPPDEPS:=$(CPPSRCS:%.cpp=$(OUTDIR)%.d)
CPPDEPS_ORG:=$(CPPSRCS_ORG:%.cpp=$(OUTDIR)%_org.d)
CCDEPS:=$(CCSRCS:%.cc=$(OUTDIR)%.d)
CCDEPS_ORG:=$(CCSRCS:%.cc=$(OUTDIR)%_org.d)
DEPS:=$(CDEPS) $(CDEPS_ORG) $(CPPDEPS) $(CPPDEPS_ORG) $(CCDEPS) $(CCDEPS_ORG)

.PHONY: clean distclean run

all: n_queen n_queen_org powerset_union powerset_union_org min_logic min_logic_org lcm lcm_org

-include $(DEPS)

n_queen: $(COBJS) $(CPPOBJS) ./app/n_queen.cpp
	$(CXX) -DUSE_DENSE $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)/n_queen $^ $(LIBS)

n_queen_org: $(COBJS_ORG) $(CPPOBJS_ORG) ./app/n_queen.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)/n_queen_org $^ $(LIBS)

powerset_union: $(COBJS) $(CPPOBJS) ./app/powerset_union.cpp
	$(CXX) -DUSE_DENSE $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)/powerset_union $^ $(LIBS)

powerset_union_org: $(COBJS_ORG) $(CPPOBJS_ORG) ./app/powerset_union.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)/powerset_union_org $^ $(LIBS)

min_logic: $(COBJS) $(CPPOBJS) $(CCOBJS) ./app/min_logic.cpp
	$(CXX) -DUSE_DENSE $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)min_logic $^ $(LIBS)

min_logic_org: $(COBJS_ORG) $(CPPOBJS_ORG) $(CCOBJS_ORG) ./app/min_logic.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)min_logic_org $^ $(LIBS)

lcm: ./app/lcm.cc $(OUTDIR)$(SRC)lcm-vsop.o $(COBJS) $(CPPOBJS) $(CCOBJS)
	$(CXX) -DUSE_DENSE $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)lcm $^ $(LIBS)

lcm_org: $(COBJS_ORG) $(CPPOBJS_ORG) $(CCOBJS_ORG) $(OUTDIR)$(SRC)lcm-vsop_org.o ./app/lcm.cc
	$(CXX) $(CXXFLAGS) $(DEFS) $(LDFLAGS) -o $(OUTDIR)lcm_org $^ $(LIBS)

$(OUTDIR)$(SRC)lcm-vsop.o: $(SRC)lcm-vsop.cc
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) -O3 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DB_64 -c $< -o $@ -D_NO_MAIN_ -DUSE_DENSE

$(OUTDIR)$(SRC)lcm-vsop_org.o: $(SRC)lcm-vsop.cc
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) -O3 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DB_64 -c $< -o $@ -D_NO_MAIN_

$(OUTDIR)%.o: %.c
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) $(DEFS) -o $@ -c -MMD -MP -MF $(@:%.o=%.d) $<

$(OUTDIR)%.o: %.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) $(DEFS) -o $@ -c -MMD -MP -MF $(@:%.o=%.d) $<

$(OUTDIR)%.o: %.cc
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) $(DEFS) -DUSE_DENSE -o $@ -c -MMD -MP -MF $(@:%.o=%.d) $<

$(OUTDIR)%_org.o: %.cc
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) $(DEFS) -o $@ -c -MMD -MP -MF $(@:%.o=%.d) $<
	
clean:
	rm -rf $(OUTDIR)

distclean:
	rm -rf ./build/
