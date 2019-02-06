*DECK UVEC
      SUBROUTINE UVEC (X, Y, YP)
C***BEGIN PROLOGUE  UVEC
C***PURPOSE  Dummy routine for BVSUP quick check.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (UVEC-S, DUVEC-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   This routine is never called;  it is here to prevent loaders from
C   complaining about undefined externals while testing BVSUP.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920401  Variables declaration and TYPE sections added.  (WRB)
C***END PROLOGUE  UVEC
C     .. Scalar Arguments ..
      REAL X
C     .. Array Arguments ..
      REAL Y(*), YP(*)
C***FIRST EXECUTABLE STATEMENT  UVEC
      STOP
      END
