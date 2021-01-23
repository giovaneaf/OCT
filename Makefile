
DIR=opt/gurobi911/linux64
FLAGS=-std=c++11 -O2 -Wall -Wextra
INCLUDE=-I/$(DIR)/include/
LIBRARY=-L/$(DIR)/lib/

all: evolutionary.cpp
	g++ $(FLAGS) evolutionary.cpp -o evolutionary $(INCLUDE) $(LIBRARY) -lgurobi_c++ -lgurobi91

heuristic: MRCTheuristic.cpp
	g++ $(FLAGS) MRCTheuristic.cpp -o heuristic

grb: mipgrb.cpp
	g++ $(FLAGS) mipgrb.cpp -o mipgrb $(INCLUDE) $(LIBRARY) -lgurobi_c++ -lgurobi91

clean:
	rm evolutionary