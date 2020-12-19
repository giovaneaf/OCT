
DIR=opt/gurobi911/linux64

all: evoalgorithm.cpp
	g++ -std=c++11 -O2 evoalgorithm.cpp -o evoalgorithm -I/$(DIR)/include/ -L/$(DIR)/lib/ -lgurobi_c++ -lgurobi91
clean: evoalgorithm
	rm evoalgorithm