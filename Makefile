IDIR=include
GCC=g++
LIBS=-lgsl -lgslcblas -lgmp
CPPFLAGS= -I$(IDIR)   -Wall  -std=c++17 -O3 -g -rdynamic -march=native 
SOURCES=$(wildcard src/*.cpp)
OBJ_1=$(patsubst %.cpp,%.o,$(SOURCES))
OBJ=$(patsubst src/%, obj/%, $(OBJ_1))
DEP1= $(patsubst %.cpp,%.d,$(SOURCES))
DEP=$(patsubst src/%, dep/%, $(DEP1))
$(shell mkdir -p obj)
$(shell mkdir -p dep)
dep/%.d: src/%.cpp
	@set -e; rm -f $@; \
	$(GCC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

obj/%.o: src/%.cpp
	$(GCC) -c  $(CPPFLAGS) -o $@ $<


Klnll: $(OBJ)
	$(GCC)  -o $@ $^ $(CPPFLAGS) $(LIBS)



.PHONY: clean

clean :
	rm -f Klnll $(OBJ) $(DEP)
