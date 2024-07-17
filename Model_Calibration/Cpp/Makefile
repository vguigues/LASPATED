GUROBI_VER = 110
INCLUDE_BOOST = /usr/include/boost

ifdef GUROBI_HOME
	USE_GUROBI=1
	INCLUDE_GUROBI = -I$(GUROBI_HOME)/include
	LIB_GUROBI = -L$(GUROBI_HOME)/lib -lgurobi_c++ $(GUROBI_LIB) -lgurobi$(GUROBI_VER)
else
	USE_GUROBI=0
endif

CARGS   	= -m64 -std=c++14 -Wfatal-errors -Wno-ignored-attributes $(CXX_FLAGS) -lboost_program_options
# CONFIG 		= -g -O0
CONFIG 		= -O3
CPP			= g++
EXE			= laspated
TEST_EXE	= test_problems

all: laspated.cpp
	$(CPP) -o $(EXE) laspated.cpp $(CARGS) -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) $(INCLUDE_LASPATED) -I. $(LIB_GUROBI)  $(CONFIG)


test: $(TEST_EXE)

$(TEST_EXE): test_problems.cpp
	$(CPP) -o $(TEST_EXE) test_problems.cpp -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) -I. $(LIB_GUROBI) -lboost_program_options $(CONFIG)

clean:
	$(RM) $(EXE) $(TEST_EXE)
