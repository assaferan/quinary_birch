bin_PROGRAMS = birch birch_dbg
lib_LTLIBRARIES = libbirch.la

AM_CXXFLAGS = -O3 -DNDEBUG -Wall -Werror -std=c++11
DBG_AM_CXXFLAGS = -g -pg -Wall -Werror -std=c++11 -DDEBUG # -DDEBUG_LEVEL_FULL
AM_CXXFLAGS += -fvar-tracking-assignments-toggle
AM_CXXFLAGS += -fomit-frame-pointer
AM_CXXFLAGS += -funroll-all-loops
DBG_AM_LDFLAGS = -lm -lgmp -lgmpxx

SOURCES = birch.h
SOURCES += birch_util.cpp
SOURCES += birch_util.h
SOURCES += EuclideanDomain.h
SOURCES += EuclideanDomain.inl
SOURCES += Integer.h
SOURCES += Integer.inl
SOURCES += Ring.h
SOURCES += Ring.inl
SOURCES += testInteger.h
SOURCES += testInteger.inl

birch_dbg_SOURCES = birch.cpp $(SOURCES)
birch_dbg_CXXFLAGS = $(DBG_AM_CXXFLAGS)

birch_SOURCES = birch.cpp $(SOURCES)
birch_CXXFLAGS = $(AM_CXXFLAGS)

libbirch_la_SOURCES = $(SOURCES)
libbirch_la_CXXFLAGS = -shared -fPIC $(AM_CXXFLAGS)

birch_dbg_LDFLAGS = ${DBG_AM_LDFLAGS}

AM_LDFLAGS = ${DBG_AM_LDFLAGS}
AM_LDFLAGS += -O3 -fPIC

birch_LDFLAGS = ${AM_LDFLAGS}
libbirch_la_LDFLAGS = ${AM_LDFLAGS}


