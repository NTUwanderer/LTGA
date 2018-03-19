
CXX = g++
#CXXFLAGS = -O0 -g -std=c++11
CXXFLAGS = -O2 -Wall -march=native -std=c++11
INCLUDE = 
TLIB = -lm

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

.SUFFIXES: .cpp

.o :
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) 

#-----File Dependencies----------------------


SRC = $(SRC1)
# SRC = $(SRC1) $(SRC2) $(SRC3)

SRC1 = LTGA.cpp spin.cpp

SRC2 = sweep.cpp spin.cpp

# SRC3 = genZobrist.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

OBJ1 = $(addsuffix .o, $(basename $(SRC1)))
OBJ2 = $(addsuffix .o, $(basename $(SRC2)))
# OBJ3 = $(addsuffix .o, $(basename $(SRC3)))

all: LTGA sweep


LTGA: $(OBJ1)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ1)

sweep: $(OBJ2)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ2)

# genZobrist: $(OBJ3)
# 	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ3)

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ)

# DO NOT DELETE

spin.o: spin.h
