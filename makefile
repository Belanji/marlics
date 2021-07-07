########### Gnu:
 COMPILER :=g++
 FLAGS    :=-Ofast -fipa-pta -fopenmp -static -march=native
 LIB      :=-lm -lgsl -lgslcblas -lgomp -pthread

############ Intel:
# COMPILER := icpc
# FLAGS    := -static -ipo -fast -no-prec-div -xHost -simd -qopenmp -fp-model fast=2
# LIB      := -mkl -lgsl 


CPPS   := $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJS   := $(patsubst src/%.cpp,build/%.o,${CPPS})

marlics: ${OBJS}
	@${COMPILER} ${FLAGS} ${OBJS} ${LIB} -o marlics
######### Build all cpp's into objects
${OBJS}: build/%.o: src/%.cpp | build
	@${COMPILER} ${FLAGS} -c $< -o $@
######### 	
${OBJS}: ${HEADER}
build:
	@mkdir build
	
######### Compile a debug version	
debug:	$(patsubst %.o,%dbg.o,${OBJS})
	@${COMPILER}  -O0 -g $(filter-out -O% -fast -static, ${FLAGS})  $(patsubst %.o,%dbg.o,${OBJS}) ${LIB} -o marlics_debug
	@rm build/*dbg.o
$(patsubst %.o,%dbg.o,${OBJS}): build/%dbg.o: src/%.cpp | build
	@${COMPILER}  -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@
$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

########
clean:
	@rm -f marlics marlics_debug
	@rm -fr build
	@rm -f *.csv
	@rm -f *.o
