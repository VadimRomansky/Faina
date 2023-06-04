CC=mpic++
CFLAGS=-I -O2 -fopenmp
DEPS = bremsstrahlung.h constants.h coordinateTransform.h examples.h inverseCompton.h massiveParticleDistribution.h optimization.h particleDistribution.h photonDistribution.h pionDecay.h radiation.h radiationSource.h synchrotron.h util.h
OBJ = bremsstrahlung.o coordinateTransform.o examples.o inverseCompton.o massiveParticleDistribution.o optimization.o photonDistribution.o pionDecay.o radiationSource.o synchrotron.o util.o main.o

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