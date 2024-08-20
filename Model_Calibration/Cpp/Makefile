GUROBI_VER = 110
INCLUDE_BOOST = /usr/include/boost

ifeq ($(USE_GUROBI), 1)
	INCLUDE_GUROBI = -I$(GUROBI_HOME)/include
	LIB_GUROBI = -L$(GUROBI_HOME)/lib -lgurobi_c++ $(GUROBI_LIB) -lgurobi$(GUROBI_VER)
else
	USE_GUROBI=0
endif

# Detectar o SO
UNAME_S := $(shell uname -s)

# Variáveis padrão
OS := Unknown
ARCH := Unknown
CFLAGS :=

# Detectar o SO e a arquitetura
ifeq ($(UNAME_S),Linux)
    OS := Linux
    ARCH := $(shell uname -m)
else ifeq ($(UNAME_S),Darwin)
    OS := macOS
    ARCH := $(shell uname -m)
else ifeq ($(OS),Windows_NT)
    OS := Windows
    ARCH := $(shell wmic os get osarchitecture | findstr /C:"64" > nul && echo x86_64 || echo x86)
endif

# Configurar CFLAGS com base na arquitetura detectada
ifneq ($(filter x86_64,$(ARCH)),)
    CFLAGS += -m64
endif

# Exibir as variáveis detectadas (para depuração)
$(info OS detected: $(OS))
$(info Architecture detected: $(ARCH))
$(info CFLAGS: $(CFLAGS))

CARGS   	= -std=c++14 -Wfatal-errors -Wno-ignored-attributes $(CXX_FLAGS) -lboost_program_options
# CONFIG 		= -g -O0
CONFIG 		= -O3
CPP			= g++
EXE			= laspated
TEST_EXE	= test_problems

all: laspated.cpp
	$(CPP) -o $(EXE) laspated.cpp $(CARGS) $(CFLAGS) -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) $(INCLUDE_LASPATED) -I. $(LIB_GUROBI)  $(CONFIG)


test: $(TEST_EXE)

$(TEST_EXE): test_problems.cpp
	$(CPP) -o $(TEST_EXE) test_problems.cpp -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) -I. $(LIB_GUROBI) -lboost_program_options $(CONFIG)

clean:
	$(RM) $(EXE) $(TEST_EXE)
