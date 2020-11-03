
ifeq (${CC},cc)
 ifeq (, $(shell which icpc))
  CC=g++
  ifeq (, $(shell which g++))
   $(error "No g++ or icpc in $$PATH.)
  endif
 else
  CC=icpc
 endif
endif

ifeq (${CC},g++)
########### Gnu:
 COMPILER := g++
 FLAGS    :=  -Ofast -fipa-pta -fopenmp -Wno-unused-result -Wno-format -static -march=znver1 -mprefer-vector-width=256
 LIB := -lm -lgsl -lgslcblas -lgomp -pthread
endif

ifeq ($(CC),icpc)  
 ifeq (, $(shell which icpc))
  $(error "No icpc in $$PATH.)
 endif
############ Intel:
 COMPILER := icpc
 FLAGS    := -ipo -fast -no-prec-div -xHost -simd -qopenmp -fp-model fast=2
 LIB := -mkl -lgsl 
endif
ifdef static
FLAGS:=${FLAGS} -static
endif

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
debug:	$(patsubst %.o,%dbg.o,${OBJS})
	@${COMPILER}  -O0 -g $(filter-out -O% -fast -static, ${FLAGS})  $(patsubst %.o,%dbg.o,${OBJS}) ${LIB} -o marlics_debug
	@rm build/*dbg.o
$(patsubst %.o,%dbg.o,${OBJS}): build/%dbg.o: src/%.cpp | build
	@${COMPILER}  -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@
$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

clean:
	@rm -f marlics marlics_debug
	@rm -fr build
	@rm -f *.csv
	@rm -f *.o
