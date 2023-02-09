TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        coordinateTransform.cpp \
        inverseCompton.cpp \
        main.cpp \
        massiveParticleDistribution.cpp \
        optimization.cpp \
        photonDistribution.cpp \
        pionDecay.cpp \
        radiationSource.cpp \
        synchrotron.cpp \
        util.cpp

HEADERS += \
    constants.h \
    coordinateTransform.h \
    inverseCompton.h \
    massiveParticleDistribution.h \
    optimization.h \
    particleDistribution.h \
    photonDistribution.h \
    pionDecay.h \
    radiationSource.h \
    synchrotron.h \
    util.h
