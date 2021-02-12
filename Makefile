
DIR=opt/gurobi911/linux64
FLAGS=-std=c++11 -O2 -Wall -Wextra
INCLUDE=-I/$(DIR)/include/
LIBRARY=-L/$(DIR)/lib/

all: evolutionary.cpp
	g++ $(FLAGS) evolutionary.cpp -o evolutionary $(INCLUDE) $(LIBRARY) -lgurobi_c++ -lgurobi91

heuristic: MRCTheuristic.cpp
	g++ $(FLAGS) MRCTheuristic.cpp -o heuristic

gls: gls.cpp
	g++ $(FLAGS) gls.cpp -o gls

simpleEvo: simpleEvo.cpp
	g++ $(FLAGS) simpleEvo.cpp -o simpleEvo

clean:
	rm evolutionary heuristic gls simpleEvo