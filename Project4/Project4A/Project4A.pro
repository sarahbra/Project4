QT += core
QT -= gui

CONFIG += c++11

TARGET = Project4A
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

INCLUDEPATH += Lib
INCLUDEPATH += /etc/alternatives/mpi
SOURCES += main.cpp
SOURCES += lib.cpp
