CC_MONO=g++ -Wall -mcmodel=medium

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: sltraj

vectorXYZ.o: vectorXYZ.cpp vectorXYZ.h
	$(CC) -c vectorXYZ.cpp 

ioutil.o: ioutil.cpp ioutil.h
	$(CC) -c ioutil.cpp 

velocity.o: velocity.cpp velocity.h
	$(CC) -c velocity.cpp -lnetcdf_c++ -fopenmp

constants.o: constants.cpp constants.h
	$(CC) -c constants.cpp

lagrangian_engine.o: lagrangian_engine.cpp ioutil.o
	$(CC) -c lagrangian_engine.cpp

sltraj.o: sltraj.cpp 
	$(CC) -c sltraj.cpp

sltraj: sltraj.o vectorXYZ.o velocity.o lagrangian_engine.o ioutil.o constants.o
	$(CC)  sltraj.o ioutil.o vectorXYZ.o velocity.o lagrangian_engine.o constants.o -o sltraj -lnetcdf_c++  -fopenmp

clean:
	$(RM) *.o 
