TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

SOURCES += \
        Src/KPIevaluator.cpp \
        Src/bremsstrahlung.cpp \
        Src/coordinateTransform.cpp \
        Src/examples.cpp \
        Src/inverseCompton.cpp \
        Src/main.cpp \
        Src/massiveParticleDistribution.cpp \
        Src/optimization.cpp \
        Src/photonDistribution.cpp \
        Src/pionDecay.cpp \
        Src/radiation.cpp \
        Src/radiationSource.cpp \
        Src/radiationSourceFactory.cpp \
        Src/synchrotron.cpp \
        Src/util.cpp

HEADERS += \
    Src/KPIevaluator.h \
    Src/bremsstrahlung.h \
    Src/constants.h \
    Src/coordinateTransform.h \
    Src/examples.h \
    Src/inverseCompton.h \
    Src/massiveParticleDistribution.h \
    Src/optimization.h \
    Src/particleDistribution.h \
    Src/photonDistribution.h \
    Src/pionDecay.h \
    Src/radiation.h \
    Src/radiationSource.h \
    Src/radiationSourceFactory.h \
    Src/synchrotron.h \
    Src/util.h
