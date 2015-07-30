#-------------------------------------------------
#
# Project created by QtCreator 2015-07-29T20:54:44
#
#-------------------------------------------------

QT       -= gui

TARGET = MLTools
TEMPLATE = lib

DEFINES += MLTOOLS_LIBRARY

SOURCES += mltools.cpp

HEADERS += mltools.h\
        mltools_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
