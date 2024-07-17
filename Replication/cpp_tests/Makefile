GUROBI_VER = 110
ifdef GUROBI_HOME
	USE_GUROBI=1
	INCLUDE_GUROBI = -I$(GUROBI_HOME)/include
	LIB_GUROBI = -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi$(GUROBI_VER)
else
	USE_GUROBI=0
endif
INCLUDE_BOOST = /usr/include/boost
INCLUDE_LASPATED = -I../../Model_Calibration/Cpp

CARGS   	= -m64 -std=c++14 -Wfatal-errors -Wno-ignored-attributes $(CXX_FLAGS) -lboost_program_options
# CONFIG 		= -g -O0
CONFIG 		= -O3
CPP			= g++
EXE			= laspated

all: main.cpp tests.h
	$(CPP) -o $(EXE) main.cpp $(CARGS) -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) $(INCLUDE_LASPATED) $(LIB_GUROBI) -I. $(CONFIG)

clean:
	$(RM) $(EXE)
