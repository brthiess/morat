.PHONY: clean fresh run gendeps

#CXX=clang++
CPPFLAGS += -std=c++0x -Wall
LDFLAGS  += -pthread

SYS := $(shell $(CXX) -dumpmachine)
ifneq (, $(findstring linux, $(SYS)))
	ALARM = lib/alarm.o
else
	ALARM = lib/alarm-timer.o
endif

ifdef DEBUG
	CPPFLAGS += -g3
else
	CPPFLAGS += -O3 -funroll-loops

	ifneq (, $(findstring darwin, $(SYS)))
		CPPFLAGS += -m64
		LDFLAGS += -m64
	else
		CPPFLAGS += -march=native
	endif
endif


all: castro nhex chex moy trex pentagod

test: \
		lib/test.o \
		lib/fileio.o \
		lib/outcome.o \
		lib/outcome_test.o \
		lib/sgf_test.o \
		lib/string.o \
		lib/string_test.o \
		lib/zobrist.o \
		havannah/agentmcts.o \
		havannah/agentmctsthread.o \
		havannah/agentmcts_test.o \
		havannah/agentpns.o \
		havannah/agentpns_test.o \
		havannah/board_test.o \
		hex/agentmcts.o \
		hex/agentmctsthread.o \
		hex/agentmcts_test.o \
		hex/agentpns.o \
		hex/agentpns_test.o \
		hex/board_test.o \
		pentago/agentmcts.o \
		pentago/agentmctsthread.o \
		pentago/agentmcts_test.o \
		pentago/agentpns.o \
		pentago/agentpns_test.o \
		pentago/board.o \
		rex/agentmcts.o \
		rex/agentmctsthread.o \
		rex/agentmcts_test.o \
		rex/agentpns.o \
		rex/agentpns_test.o \
		cylindrical_hex/agentmcts.o \
		cylindrical_hex/agentmctsthread.o \
		cylindrical_hex/agentmcts_test.o \
		cylindrical_hex/agentpns.o \
		cylindrical_hex/agentpns_test.o \
		cylindrical_hex/board_test.o \
		y/agentmcts.o \
		y/agentmctsthread.o \
		y/agentmcts_test.o \
		y/agentpns.o \
		y/agentpns_test.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)
	./test

castro: \
		havannah/main.o \
		havannah/agentmcts.o \
		havannah/agentmctsthread.o \
		havannah/agentpns.o \
		havannah/gtpgeneral.o \
		havannah/gtpagent.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		lib/zobrist.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

pentagod: \
		pentago/main.o \
		pentago/agentab.o \
		pentago/agentmcts.o \
		pentago/agentmctsthread.o \
		pentago/agentpns.o \
		pentago/board.o \
		pentago/gtpgeneral.o \
		pentago/gtpagent.o \
		pentago/moveiterator.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

moy: \
		y/main.o \
		y/agentmcts.o \
		y/agentmctsthread.o \
		y/agentpns.o \
		y/gtpagent.o \
		y/gtpgeneral.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		lib/zobrist.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

nhex: \
		hex/main.o \
		hex/agentmcts.o \
		hex/agentmctsthread.o \
		hex/agentpns.o \
		hex/gtpagent.o \
		hex/gtpgeneral.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		lib/zobrist.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)
	
chex: \
		cylindrical_hex/main.o \
		cylindrical_hex/agentmcts.o \
		cylindrical_hex/agentmctsthread.o \
		cylindrical_hex/agentpns.o \
		cylindrical_hex/gtpagent.o \
		cylindrical_hex/gtpgeneral.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		lib/zobrist.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

trex: \
		rex/main.o \
		rex/agentmcts.o \
		rex/agentmctsthread.o \
		rex/agentpns.o \
		rex/gtpagent.o \
		rex/gtpgeneral.o \
		lib/fileio.o \
		lib/gtpcommon.o \
		lib/outcome.o \
		lib/string.o \
		lib/zobrist.o \
		$(ALARM)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

clean:
	rm -f */*.o test castro moy pentagod nhex chex trex .Makefile

fresh: clean all

profile:
	valgrind --tool=callgrind

gendeps: .Makefile

.Makefile: # contains the actual dependencies for all the .o files above
	./gendeps.sh > .Makefile

include .Makefile
