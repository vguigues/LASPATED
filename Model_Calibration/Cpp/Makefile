GUROBI_LIB_VER = -lgurobi110
INCLUDE_BOOST = /usr/include/boost

ifdef GUROBI_HOME
	USE_GUROBI=1
	INCLUDE_GUROBI = -I$(GUROBI_HOME)/include
	LIB_GUROBI = -L$(GUROBI_HOME)/lib -lgurobi_c++ $(GUROBI_LIB_VER)
else
	ifdef GUROBI_INCLUDE
		ifdef GUROBI_LIB
			USE_GUROBI=1
			INCLUDE_GUROBI = -I$(GUROBI_INCLUDE)
			LIB_GUROBI = -L$(GUROBI_LIB)
		else
			$(error Missing GUROBI_INCLUDE or GUROBI_LIB. Please provide both or none of them.)
		endif
	else
		ifdef GUROBI_LIB
			$(error Missing GUROBI_INCLUDE or GUROBI_LIB. Please provide both or none of them.)	
		else
			USE_GUROBI=0
		endif
	endif
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
