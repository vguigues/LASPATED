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

all: main.cpp tests.h
	$(CPP) -o $(EXE) main.cpp $(CARGS) $(CFLAGS) -DUSE_GUROBI=$(USE_GUROBI) $(INCLUDE_GUROBI) $(INCLUDE_LASPATED) $(LIB_GUROBI) -I. $(CONFIG)

clean:
	$(RM) $(EXE)
