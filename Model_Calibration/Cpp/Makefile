GUROBI_VER = 110
INCLUDE_BOOST = /usr/include/boost

ifdef GUROBI_HOME
	USE_GUROBI=1
	INCLUDE_GUROBI = -I$(GUROBI_HOME)/include
	LIB_GUROBI = -L$(GUROBI_HOME)/lib -lgurobi_c++ $(GUROBI_LIB) -lgurobi$(GUROBI_VER)
else
	USE_GUROBI=0
endif

# ifndef XTL_INCLUDE
# $(error undefined XTL_INCLUDE. Please provide xtl location)
# else 
# INCLUDE_XTL = -I$(XTL_INCLUDE)
# endif
# ifndef XTENSOR_INCLUDE
# $(error undefined XTENSOR_INCLUDE. Please provide xtensor location)
# else
# INCLUDE_XTENSOR = -I$(XTENSOR_INCLUDE)
# endif

CARGS   	= -m64 -std=c++14 -Wfatal-errors -Wno-ignored-attributes $(CXX_FLAGS) -lboost_program_options
# CONFIG 		= -g -O0
CONFIG 		= -O3
CPP			= g++
EXE			= laspated

all: laspated.cpp
	$(CPP) -o $(EXE) laspated.cpp $(CARGS) -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) $(INCLUDE_LASPATED) $(LIB_GUROBI) $(CONFIG)

clean:
	$(RM) $(EXE)
