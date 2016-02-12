QT += core
QT += gui

TARGET = MRVM
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXX = clang-omp++
QMAKE_CC = clang-omp
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
QMAKE_CXXFLAGS_RELEASE += $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$QMAKE_CXXFLAGS

QMAKE_LFLAGS += -F/Library/Frameworks/

INCLUDEPATH += /usr/local/include \
               /usr/local/include/libiomp

DEPENDPATH += /usr/local/include

QMAKE_LIBS += -liomp5
LIBS += -L/usr/local/lib/ -lafopencl.3.2.2
LIBS += -framework GDAL
LIBS += -L/usr/local/lib/ -liomp5

HEADERS += ./include/mrvm.h \
           ./include/mrvm_global.h \
           ./include/roots.h

PRECOMPILED_HEADER += ./headers/stdafx.h

SOURCES += ./src/stdafx.cpp \
           ./src/main.cpp \
           ./src/kernel.cpp \
           ./src/roots.cpp \
           ./src/mrvm.cpp \
           ./src/mrvmitem.cpp \
           ./src/realmrvmitem.cpp \
           ./src/realarraymrvmitem.cpp \
           ./src/rasteritem.cpp \
           ./src/realrastermrvmitem.cpp \
           ./src/categoricalmrvmitem.cpp \
           ./src/categoricalrastermrvmitem.cpp \
           ./src/rasterbootstrap.cpp


