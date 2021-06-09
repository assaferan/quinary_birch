bin_PROGRAMS = birch birch_dbg
lib_LTLIBRARIES = libbirch.la

OPT_CXXFLAGS = -O3 -DNDEBUG -Wall -Werror -std=c++11
DBG_CXXFLAGS = -g -pg -Wall -Werror -std=c++11 -DDEBUG # -DDEBUG_LEVEL_FULL
OPT_CXXFLAGS += -fvar-tracking-assignments-toggle
OPT_CXXFLAGS += -fomit-frame-pointer
OPT_CXXFLAGS += -funroll-all-loops
DBG_LDFLAGS = -lm -lgmp -lgmpxx

SOURCES = birch.h
SOURCES += birch_util.cpp
SOURCES += birch_util.h
SOURCES += EuclideanDomain.h
SOURCES += EuclideanDomain.inl
SOURCES += Integer.h
SOURCES += Integer.inl
SOURCES += Rational.h
SOURCES += Rational.inl
SOURCES += Ring.h
SOURCES += Ring.inl
SOURCES += testInteger.h
SOURCES += testInteger.inl

birch_dbg_SOURCES = birch.cpp $(SOURCES)
birch_dbg_CXXFLAGS = $(DBG_CXXFLAGS)

birch_SOURCES = birch.cpp $(SOURCES)
birch_CXXFLAGS = $(OPT_CXXFLAGS)

libbirch_la_SOURCES = $(SOURCES)
libbirch_la_CXXFLAGS = -shared -fPIC $(OPT_CXXFLAGS)

birch_dbg_LDFLAGS = ${DBG_LDFLAGS}

OPT_LDFLAGS = ${DBG_LDFLAGS}
OPT_LDFLAGS += -O3 -fPIC

birch_LDFLAGS = ${OPT_LDFLAGS}
libbirch_la_LDFLAGS = ${OPT_LDFLAGS}


