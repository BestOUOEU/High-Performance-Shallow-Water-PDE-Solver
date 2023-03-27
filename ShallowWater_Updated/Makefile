CC       = g++
CXX      = g++
CXXFLAGS = -g -std=c++11 -fopenmp -O3
CCFLAGS  = -g -std=c++11 -fopenmp -O3
LDFLAGS  = -g -std=c++11 -fopenmp -O3
LDLIBS   = -llapack -lblas -lboost_program_options
TARGET   = ShallowWaterMain
TARGET2  = ShallowWater
HDRS     = ShallowWater.h

default: $(TARGET)

$(TARGET): $(TARGET).o $(TARGET2).o

$(TARGET).o: $(TARGET).cpp 

$(TARGET2).o: $(TARGET2).cpp $(HDRS)

.PHONY: test1 test2 test3 test4 clean git plot

test1: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 1 --T 80 --dt 0.1 --case1 0

test2: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 2 --T 80 --dt 0.1 --case1 0

test3: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 3 --T 80 --dt 0.1 --case1 0

test4: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 4 --T 80 --dt 0.1 --case1 0

test5: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 1 --T 80 --dt 0.1 --case1 1

test6: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 2 --T 80 --dt 0.1 --case1 1

test7: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 3 --T 80 --dt 0.1 --case1 1

test8: $(TARGET)
	time ./$(TARGET) --Nx 100 --Ny 100 --ic 4 --T 80 --dt 0.1 --case1 1

clean:
	rm -f $(TARGET) *.o
	rm -f $(TARGET2) *.o
	rm -f $(TARGET)
	rm -f repository.log

git:
	git log --name-status > repository.log
