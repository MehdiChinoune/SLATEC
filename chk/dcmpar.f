*DECK DCMPAR
      SUBROUTINE DCMPAR (ICNT, ITEST)
C***BEGIN PROLOGUE  DCMPAR
C***PURPOSE  Compare values in COMMON block DCHECK for quick check
C            routine DPFITT.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (CMPARE-S, DCMPAR-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DCHECK
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890921  Realigned order of variables in the COMMON block.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920214  Minor improvements to code for readability.  (WRB)
C***END PROLOGUE  DCMPAR
C     .. Scalar Arguments ..
      INTEGER ICNT
C     .. Array Arguments ..
      INTEGER ITEST(9)
C     .. Scalars in Common ..
      DOUBLE PRECISION EPS, RP, SVEPS, TOL
      INTEGER IERP, IERR, NORD, NORDP
C     .. Arrays in Common ..
      DOUBLE PRECISION R(11)
C     .. Local Scalars ..
      DOUBLE PRECISION RPP, SS
      INTEGER IERPP, NRDP
C     .. Local Arrays ..
      INTEGER ITEMP(4)
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     .. Common blocks ..
      COMMON /DCHECK/ EPS, R, RP, SVEPS, TOL, NORDP, NORD, IERP, IERR
C***FIRST EXECUTABLE STATEMENT  DCMPAR
      ICNT = ICNT + 1
      ITEMP(1) = 0
      ITEMP(2) = 0
      ITEMP(3) = 0
      ITEMP(4) = 0
      SS = SVEPS - EPS
      NRDP = NORDP - NORD
      RPP = RP - R(11)
      IERPP = IERP - IERR
      IF (ABS(SS).LE.TOL .OR. ICNT.LE.2 .OR. ICNT.GE.6) ITEMP(1) = 1
      IF (ABS(NRDP) .EQ. 0) ITEMP(2) = 1
      IF (ICNT .EQ. 2) ITEMP(2) = 1
      IF (ABS(RPP) .LE. TOL) ITEMP(3) = 1
      IF (ABS(IERPP) .EQ. 0) ITEMP(4) = 1
C
C     Check to see if all four tests were good.
C     If so, set the test number equal to 1.
C
      ITEST(ICNT) = ITEMP(1)*ITEMP(2)*ITEMP(3)*ITEMP(4)
      RETURN
      END
