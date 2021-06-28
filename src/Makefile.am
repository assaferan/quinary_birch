bin_PROGRAMS = birch birch_dbg
lib_LTLIBRARIES = libbirch.la

OPT_CXXFLAGS = -O3 -DNDEBUG -Wall -Werror -std=c++11 -msse3 -mavx
DBG_CXXFLAGS = -g -pg -Wall -Werror -std=c++11 -DDEBUG # -DDEBUG_LEVEL_FULL
OPT_CXXFLAGS += -fvar-tracking-assignments-toggle
OPT_CXXFLAGS += -fomit-frame-pointer
OPT_CXXFLAGS += -funroll-all-loops
DBG_LDFLAGS = -lm -lgmp -lgmpxx -L../libs -lantic -lflint -lbliss -lnauty

SOURCES = birch.h
SOURCES += birch_util.cpp
SOURCES += birch_util.h
SOURCES += birch_util.inl
SOURCES += BirchExample.h
SOURCES += BirchExample.inl
SOURCES += Eigenvector.h
SOURCES += Eigenvector.inl
SOURCES += EuclideanDomainElement.h
SOURCES += EuclideanDomainElement.inl
SOURCES += FieldElement.h
SOURCES += FieldElement.inl
SOURCES += Fp.cpp
SOURCES += Fp.h
SOURCES += Fp.inl
SOURCES += FpElement.h
SOURCES += FpElement.inl
SOURCES += Integer.cpp
SOURCES += Integer.h
SOURCES += Integer.inl
SOURCES += IntegerRing.h
SOURCES += Isometry.h
SOURCES += Isometry.inl
SOURCES += Matrix.h
SOURCES += Matrix.inl
SOURCES += NumberField.h
SOURCES += NumberField.inl
SOURCES += NumberFieldElement.h
SOURCES += NumberFieldElement.inl
SOURCES += ParseNipp.cpp
SOURCES += ParseNipp.h
SOURCES += Polynomial.h
SOURCES += Polynomial.inl
SOURCES += QuadForm.h
SOURCES += QuadForm.inl
SOURCES += QuadFormFp.h
SOURCES += QuadFormFp.inl
SOURCES += QuadFormInt.cpp
SOURCES += QuadFormInt.h
SOURCES += QuadFormInt.inl
SOURCES += Rational.h
SOURCES += Rational.inl
SOURCES += RationalField.h
SOURCES += Ring.h
SOURCES += RingElement.h
SOURCES += RingElement.inl
SOURCES += SetCover.cpp
SOURCES += SetCover.h
SOURCES += SquareMatrix.h
SOURCES += SquareMatrix.inl
SOURCES += SquareMatrixInt.cpp
SOURCES += SquareMatrixInt.h
SOURCES += SquareMatrixInt.inl
SOURCES += TestBirch.h
SOURCES += TestBirch.inl
SOURCES += testInteger.h
SOURCES += testInteger.inl
SOURCES += testRational.h
SOURCES += testRational.inl
SOURCES += UnivariatePoly.h
SOURCES += UnivariatePoly.inl
SOURCES += UnivariatePolyFp.h
SOURCES += UnivariatePolyFp.inl
SOURCES += UnivariatePolyInt.cpp
SOURCES += UnivariatePolyInt.h
SOURCES += UnivariatePolyInt.inl
SOURCES += Vector.h
SOURCES += Vector.inl
SOURCES += VectorInt.h
SOURCES += VectorInt.inl
SOURCES += polyhedral_common/src_basic/Basic_string.cpp
SOURCES += polyhedral_common/src_basic/ExceptionEnding.cpp
SOURCES += polyhedral_common/src_basic/hash_functions.cpp
SOURCES += polyhedral_common/src_basic/Temp_common.cpp
SOURCES += polyhedral_common/src_basic/TypeConversion.cpp
SOURCES += polyhedral_common/src_number/NumberTheory.cpp

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


