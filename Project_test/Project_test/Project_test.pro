QT += core
QT -= gui

CONFIG += c++11

TARGET = Project_test
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

INCLUDEPATH += Lib
SOURCES += main.cpp
SOURCES += lib.cpp
