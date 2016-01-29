CC_MONO=g++ -Wall -mcmodel=medium

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: sltraj mkitracer idisp

vectorXYZ.o: vectorXYZ.cpp vectorXYZ.h
	$(CC) -c vectorXYZ.cpp 

mkitracer.o: mkitracer.cpp 
	$(CC) -c mkitracer.cpp 

ioutil.o: ioutil.cpp ioutil.h
	$(CC) -c ioutil.cpp 

velocity.o: velocity.cpp velocity.h
	$(CC) -c velocity.cpp -lnetcdf_c++ -fopenmp

constants.o: constants.cpp constants.h
	$(CC) -c constants.cpp

random.o: random.cpp random.h
	$(CC) -c random.cpp

lagrangian_engine.o: lagrangian_engine.cpp ioutil.o random.o
	$(CC) -c lagrangian_engine.cpp

sltraj.o: sltraj.cpp 
	$(CC) -c sltraj.cpp

idisp.o: idisp.cpp 
	$(CC) -c idisp.cpp

sltraj: sltraj.o vectorXYZ.o velocity.o lagrangian_engine.o ioutil.o constants.o random.o
	$(CC)  sltraj.o ioutil.o vectorXYZ.o velocity.o lagrangian_engine.o constants.o random.o -o sltraj -lnetcdf_c++  -fopenmp

idisp: idisp.o vectorXYZ.o ioutil.o constants.o
	$(CC)  idisp.o ioutil.o vectorXYZ.o constants.o -o idisp

mkitracer: mkitracer.o  vectorXYZ.o constants.o
	$(CC)  mkitracer.o  vectorXYZ.o constants.o -o mkitracer

clean:
	$(RM) *.o 
