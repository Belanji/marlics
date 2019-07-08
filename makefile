########### Gnu:
#COMPILER = g++
#FLAGS= -static -Ofast -fipa-pta -fopenmp
#LIB = -lm -lgsl -lgslcblas -lgomp -pthread

############ Intel:

#COMPILER = icpc
#FLAGS=     -ipo -O3  -no-prec-div -xHost -simd -qopenmp -fp-model fast=2  -g
#LIB = -mkl -lgsl 

############ Intel 2 :

COMPILER = icpc
FLAGS= -ipo -O3  -no-prec-div -xAVX -simd -qopenmp -fp-model fast=2 -static
LIB = -mkl -lgsl 

############# Pgi:
#FLAGS =  -fast -Mipa=fast,inline -O3
CPPS := $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJS  := $(patsubst src/%.cpp,build/%.o,${CPPS})

tensorQ: ${OBJS}
	@${COMPILER}  ${FLAGS}  ${OBJS} ${LIB} -o tensorQ
${OBJS}: build/%.o: src/%.cpp | build
	@${COMPILER}  ${FLAGS} -c $< -o $@
${OBJS}: ${HEADER}
build:
	@mkdir build

clean:
	@rm -f tensorQ
	@rm -fr build
	@rm -f *.csv
	@rm -f *.oo