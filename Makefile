#CC = clang++
CXXFLAGS = -O3 -mcx16 -march=native -std=c++20 -Wall -Wextra
CC = clang++
INCLUDE_PATH = -Iparlaylib/include/
all: build_pch query
.PHONY: gen run clean

ifdef CILKPLUS
CC = clang++
CXXFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CXXFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CXXFLAGS += -DPARLAY_SEQUENTIAL
else
CXXFLAGS += -pthread
endif

ifdef DEBUG
CC = g++
CXXFLAGS += -DDEBUG -Og
else ifdef PERF
CC = g++
CXXFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CXXFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL
else
CXXFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CXXFLAGS += -DPARLAY_USE_STD_ALLOC
endif

build_pch:	build_pch.cpp build_pch.hpp graph.hpp dijkstra.hpp
	$(CC) $(CXXFLAGS) $(INCLUDE_PATH) build_pch.cpp -o build_pch

query:	query_test.cpp graph.hpp query.hpp
	$(CC) $(CXXFLAGS) $(INCLUDE_PATH) query_test.cpp -o query

clean:
	rm build_pch
	rm query
