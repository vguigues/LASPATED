
CXX_FLAGS =  -DBOOST_TEST_DYN_LINK -DBOOST_SPIRIT_USE_PHOENIX_V3 \
			-DBOOST_RESULT_OF_USE_DECLTYPE -DBOOST_FILESYSTEM_NO_DEPRECATED -rdynamic


BOOST_INCLUDE 	= -I/usr/local/include/boost
XTL_INCLUDE 	= -I/usr/local/include/xtl
XTENSOR_INCLUDE = -I/usr/local/include/xtensor
# FMT_INCLUDE     = -I/usr/local/include/fmt
GUROBI_INCLUDE 	= -I$(GUROBI_HOME)/include

PLATFORM	= linux64
# INC     	= -I/usr/local/include/boost -I$(GUROBI_HOME)/include -I/usr/include/lua5.2 \
# 			-I/usr/local/include -I/usr/local/include/osrm -I/usr/local/include/xtensor \
# 			-I/usr/local/include/xtl -I/usr/local/include/fmt
INC     	:= $(BOOST_INCLUDE) $(FMT_INCLUDE) \
			$(XTL_INCLUDE) $(XTENSOR_INCLUDE) $(GUROBI_INCLUDE) \

CPP     	= /usr/bin/g++-9
CARGS   	= -m64 -std=c++14 -Wfatal-errors -Wno-ignored-attributes $(CXX_FLAGS)
CPPLIB  	= -L$(GUROBI_HOME)/lib -lboost_program_options -lboost_system -pthread -lpthread -ldl \
			-lgurobi_c++ -lgurobi100 -L/usr/local/lib -lrt -lz -fuse-ld=gold \
			-Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-O1 -Wl,--hash-style=gnu \
			-Wl,--sort-common -lboost_date_time -lboost_chrono -lboost_filesystem \
			-lboost_iostreams -lboost_thread -lboost_system -ltbb -ltbbmalloc -lfmt
INCLUDE_DIR	= include
SRC_DIR		= src
OBJ_DIR		= obj
SRC 		= $(wildcard $(SRC_DIR)/*.cpp)
OBJ 		= $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
EXE 		= laspated
# CONFIG		= -O1 -g
# CONFIG		= -O0 -g
CONFIG 		= -O3

all: $(EXE)

$(EXE): $(OBJ)
	$(CPP) $^ -lm -o $@ $(CPPLIB) $(CONFIG)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDE_DIR)/%.h
	$(CPP) $(CARGS) -c $< -o $@ $(INC) $(CPPLIB) $(CONFIG)

clean:
	$(RM) $(OBJ)
