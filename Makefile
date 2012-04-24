CHARM_PATH = #
CHARM_LIBS = $(CHARM_PATH)/lib
INCPATH = $(CHARM_PATH)/include
STRUCTURES_PATH = utility/structures

VPATH = $(STRUCTURES_PATH)

APP_FLAGS = -DCHECK_INTER -DPHASE_BARRIERS -DNODE_LEVEL_MERGE #-DVERBOSE_TRAVERSAL -DVERBOSE_TRAVERSAL_INTERACTION
OPTS = -O3 -g $(APP_FLAGS)
CPPFLAGS += -I$(INCPATH) -I$(STRUCTURES_PATH) 
CXXFLAGS += $(OPTS) $(CPPFLAGS)
LDFLAGS += $(OPTS) -L$(STRUCTURES_PATH) -lTipsy -L. -language charm++ -module RandCentLB -module RotateLB -module GreedyLB -module Orb3dLB_notopo -module NDMeshStreamer -module completion -memory os #-tracemode projections

CHARMC = $(CHARM_PATH)/bin/charmc

CXX = $(CHARMC)
CC = $(CXX)
AR = ar q 
CXX_DEPEND = $(CXX) -M -MM -MG -Wall $(APP_FLAGS) 
CFLAGS = $(OPTS) $(DEFINE_FLAGS) -g 

OBJECTS = Main.o DataManager.o TreePiece.o util.o \
	  Orb3dLB_notopo.o  TreeMerger.o \
	  Reduction.o Worker.o Request.o State.o 
SRC = Main.cpp DataManager.cpp TreePiece.cpp \
      util.cpp Reduction.cpp Worker.cpp Request.cpp \
      Orb3dLB_notopo.cpp TreeMerger.cpp \
      State.cpp gen_util.cpp plummer.cpp tipsyPlummer.cpp  

TARGET = barnes
all: $(STRUCTURES_PATH)/libTipsy.a $(TARGET) 

$(TARGET): $(OBJECTS) Makefile.dep libmoduleOrb3dLB_notopo.a $(STRUCTURES_PATH)/libTipsy.a 
	$(CHARMC) -o $(TARGET) $(LDFLAGS) $(OBJECTS)


libmoduleOrb3dLB_notopo.a: Orb3dLB_notopo.o
	$(CHARMC) -o libmoduleOrb3dLB_notopo.a Orb3dLB_notopo.o  



$(STRUCTURES_PATH)/libTipsy.a: $(STRUCTURES_PATH)/Makefile
	cd $(STRUCTURES_PATH); $(MAKE) libTipsy.a

$(STRUCTURES_PATH)/Makefile:  
	cd $(STRUCTURES_PATH); ./configure

plummer.o: plummer.cpp 
	g++ $(CPPFLAGS) -c plummer.cpp

tipsyPlummer.o: tipsyPlummer.cpp 
	g++ $(CPPFLAGS) -g -O0 -c tipsyPlummer.cpp

gen_util.o: gen_util.cpp 
	g++ $(CPPFLAGS) -c gen_util.cpp

plummer: plummer.o gen_util.o 
	g++ $(CPPFLAGS) -o plummer plummer.o gen_util.o

tipsyPlummer: tipsyPlummer.o gen_util.o 
	g++ $(CPPFLAGS) -g -O0 -o tipsyPlummer tipsyPlummer.o gen_util.o -lTipsy -L$(STRUCTURES_PATH)

%.decl.h %.def.h : %.ci
	$(CHARMC) $(APP_FLAGS) -E $(CPPFLAGS) $<

%.o: Makefile

clean:
	rm -f core* *.a *.o $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host gen plummer tipsyPlummer

depends:
	$(CXX_DEPEND) $(SRC) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep 

.PHONY: all docs depends 

include Makefile.dep
