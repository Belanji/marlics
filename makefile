########### Gnu:
#COMPILER = g++
#FLAGS= -static -Ofast -fipa-pta -fopenmp
#LIB = -lm -lgsl -lgslcblas -lgomp -pthread

############ Intel:

COMPILER = icpc
FLAGS=     -ipo -O3  -no-prec-div -xHost -simd -qopenmp -fp-model fast=2  -g
LIB = -mkl -lgsl 

############ Intel 2 :

#COMPILER = icpc
#FLAGS= -ipo -O3  -no-prec-div -xAVX -simd -qopenmp -fp-model fast=2 -static
#LIB = -mkl -lgsl 


############ Intel Debug :

#COMPILER = icpc
#FLAGS= -O0 -g -static
#LIB = -mkl -lgsl 


CPPS := $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJS  := $(patsubst src/%.cpp,build/%.o,${CPPS})

marlics: ${OBJS}
	@${COMPILER}  ${FLAGS}  ${OBJS} ${LIB} -o marlics
${OBJS}: build/%.o: src/%.cpp | build
	@${COMPILER}  ${FLAGS} -c $< -o $@
${OBJS}: ${HEADER}
build:
	@mkdir build

clean:
	@rm -f marlics
	@rm -fr build
	@rm -f *.csv
	@rm -f *.oo
