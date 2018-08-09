CXX=g++
RM=rm
LDLIBS=-larmadillo

all: base dp

base: base_main.cpp
	$(CXX) base_main.cpp -O2 -std=c++11 -o base_sim $(LDLIBS)

dp: dp.cpp
	$(CXX) base_model.cpp dp_main.cpp -O2 -std=c++11 -o dp_sim $(LDLIBS)
