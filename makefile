CC=mpic++
CFLAGS=-I -O2 -fopenmp
DEPS = ./Src/bremsstrahlung.h ./Src/constants.h ./Src/coordinateTransform.h ./Src/examples.h ./Src/inverseCompton.h ./Src/massiveParticleDistribution.h ./Src/optimization.h ./Src/particleDistribution.h ./Src/photonDistribution.h ./Src/pionDecay.h ./Src/radiation.h ./Src/radiationSource.h ./Src/synchrotron.h ./Src/util.h
OBJ = bremsstrahlung.o coordinateTransform.o examples.o inverseCompton.o massiveParticleDistribution.o optimization.o photonDistribution.o pionDecay.o radiationSource.o radiation.o synchrotron.o util.o main.o

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS) -o $@ $<

# the build target executable:
TARGET = a.out

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

default:
	$(TARGET)

all:
	$(TARGET)

clean:
	rm -f *.o
	rm -f *.err
	rm -f *.out