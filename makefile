#############################################################################
# Makefile for building: Faina
#############################################################################

SRC          = ./Src
INCLUDE_DIRS = -I. -I$(SRC)
VPATH        = ./:$(SRC)/Math:$(SRC)

####### Compiler, tools and options

CC       = gcc
#CFLAGS   = -c -g -Wundef
#CFLAGS   = -c -O3 -Wundef
LDFLAGS  = -lm
CXX           = g++
#CFLAGS        = -pipe -g -Wall -W -fPIC
#CXXFLAGS      = -pipe -fopenmp -g -std=gnu++11 -Wall -W -fPIC
CFLAGS        = -pipe -O3 -Wall -W -fPIC
CXXFLAGS      = -pipe -fopenmp -O3 -std=gnu++11 -Wall -W -fPIC
INCPATH       = -I.
LINK          = g++
LFLAGS        = 
LIBS          = $(SUBLIBS) -fopenmp
DEL_FILE      = rm -f   

####### Output directory

OBJECTS_DIR   = ./

####### Files

HEADERS       = ./Src/constants.h \
		./Src/bremsstrahlung.h \
		./Src/coordinateTransform.h \
		./Src/examples.h \
		./Src/inverseCompton.h \
		./Src/KPIevaluator.h \
		./Src/massiveParticleDistribution.h \
		./Src/optimization.h \
		./Src/particleDistribution.h \
		./Src/photonDistribution.h \
		./Src/pionDecay.h \
		./Src/radiation.h \
		./Src/radiationSource.h \
		./Src/synchrotron.h \
		./Src/util.h \
		./Src/diffusion.h \
		./Src/Math/largeVectorBasis.h \
		./Src/Math/matrixElement.h \
		./Src/Math/specialmath.h

SOURCES       = Src/KPIevaluator.cpp \
		Src/bremsstrahlung.cpp \
		Src/coordinateTransform.cpp \
		Src/examples.cpp \
		Src/inverseCompton.cpp \
		Src/massiveParticleDistribution.cpp \
		Src/optimization.cpp \
		Src/photonDistribution.cpp \
		Src/pionDecay.cpp \
		Src/radiation.cpp \
		Src/radiationSource.cpp \
		Src/radiationSourceFactory.cpp \
		Src/synchrotron.cpp \
		Src/util.cpp \
		Src/diffusion.cpp \
		Src/Math/largeVectorBasis.cpp \
		Src/Math/matrixElement.cpp \
		Src/Math/specialmath.cpp \
		main.cpp 
OBJS      = KPIevaluator.o \
		bremsstrahlung.o \
		coordinateTransform.o \
		examples.o \
		inverseCompton.o \
		massiveParticleDistribution.o \
		optimization.o \
		photonDistribution.o \
		pionDecay.o \
		radiation.o \
		radiationSource.o \
		radiationSourceFactory.o \
		synchrotron.o \
		util.o \
		diffusion.o \
		largeVectorBasis.o \
		matrixElement.o \
		specialmath.o \
		main.o

TARGET        = Faina


first: Faina
####### Build rules

Faina:  $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDPATHS) $(LDFLAGS) $(LIBS)

%.o: vpath %.cpp
	$(CXX) $(CXX_FLAGS) -c -o $@ $< $(LDPATHS) $(LDFLAGS) $(LIBS)

clean:
	@rm -f	*.o
	@echo make clean: done

# ---------------------------------------------------------
#          Dependencies for object files
# ---------------------------------------------------------

$(OBJ):  $(HEADERS)

