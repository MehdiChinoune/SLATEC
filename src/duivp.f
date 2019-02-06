*DECK DUIVP
      SUBROUTINE DUIVP (X, Y, YP)
C***BEGIN PROLOGUE  DUIVP
C***PURPOSE  Dummy routine for DBVSUP quick check.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (UIVP-S, DUIVP-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This routine is never called;  it is here to prevent loaders from
C   complaining about undefined externals while testing DBVSUP.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920401  Variables declaration and TYPE sections added.  (WRB)
C***END PROLOGUE  DUIVP
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     .. Array Arguments ..
      DOUBLE PRECISION Y(*), YP(*)
C***FIRST EXECUTABLE STATEMENT  DUIVP
      STOP
      END
