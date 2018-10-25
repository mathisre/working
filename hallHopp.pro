TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    CGanalysis.cpp \
    fileutils.cpp \
    main.cpp \
    mersenne.cpp \
    MKcurrent.cpp \
    paramfile.cpp \
    stringutils.cpp \
    systemclass.cpp \
    treeutils.cpp \
    params.cpp \
    run.cpp

HEADERS += \
    CGanalysis.h \
    fileutils.h \
    MKcurrent.h \
    paramfile.h \
    randomc.h \
    stringutils.h \
    systemclass.h \
    treeutils.h \
    params.h \
    run.h \
    mersenne.h
