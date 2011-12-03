#-------------------------------------------------
#
# Project created by QtCreator 2011-12-03T12:40:37
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = vp-tree
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -fopenmp -g


SOURCES += main.cpp

HEADERS += \
    vp-tree.h

LIBS += -lgomp

