QT += core
QT -= gui

CONFIG += c++11

TARGET = Project4d
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app
INCLUDEPATH += Lib
SOURCES += main.cpp
SOURCES += Lib/lib.cpp
