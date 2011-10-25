CHARM_PATH = #
INCPATH = $(CHARM_PATH)/include
STRUCTURES_PATH = #

VPATH = $(STRUCTURES_PATH)

APP_FLAGS = #
OPTS = -O3 -g $(APP_FLAGS)
CXXFLAGS += $(OPTS) -I$(INCPATH) -I$(STRUCTURES_PATH)
LDFLAGS += $(OPTS) -L. -language charm++ -module RandCentLB -module RotateLB -module GreedyLB -memory os #-tracemode projections

CHARMC = $(CHARM_PATH)/bin/charmc

CXX = $(CHARMC)
CC = $(CXX)
AR = ar q 
CXX_DEPEND = $(CXX) -M -MM -MG -Wall $(APP_FLAGS) 
CFLAGS = $(OPTS) $(DEFINE_FLAGS) -g 

OBJECTS = Main.o DataManager.o TreePiece.o util.o \
	  Reduction.o Worker.o Request.o
SRC = Main.cc DataManager.cc TreePiece.cc \
      util.cc Reduction.cc Worker.cc Request.cc \
      gen_util.cc plummer.cc

TARGET = barnes 
all: $(TARGET) plummer

$(TARGET): $(OBJECTS) Makefile.dep 
	$(CHARMC) -o $(TARGET) $(LDFLAGS) $(OBJECTS)

plummer.o: plummer.cc 
	g++ -I$(STRUCTURES_PATH) -c plummer.cc 

gen_util.o: gen_util.cc 
	g++ -I$(STRUCTURES_PATH) -c gen_util.cc 

plummer: plummer.o gen_util.o 
	g++ -I$(STRUCTURES_PATH) -o plummer plummer.o gen_util.o

%.decl.h %.def.h : %.ci
	$(CHARMC) $(APP_FLAGS) -E $<

%.o: Makefile

clean:
	rm -f core* *.a $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host gen plummer

depends:
	$(CXX_DEPEND) $(SRC) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep 

.PHONY: all docs depends 

include Makefile.dep
