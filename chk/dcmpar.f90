!*==DCMPAR.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DCMPAR
      SUBROUTINE DCMPAR(Icnt,Itest)
      IMPLICIT NONE
!*--DCMPAR5
!***BEGIN PROLOGUE  DCMPAR
!***PURPOSE  Compare values in COMMON block DCHECK for quick check
!            routine DPFITT.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (CMPARE-S, DCMPAR-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DCHECK
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   890921  Realigned order of variables in the COMMON block.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920214  Minor improvements to code for readability.  (WRB)
!***END PROLOGUE  DCMPAR
!     .. Scalar Arguments ..
      INTEGER Icnt
!     .. Array Arguments ..
      INTEGER Itest(9)
!     .. Scalars in Common ..
      DOUBLE PRECISION EPS , RP , SVEps , TOL
      INTEGER IERp , IERr , NORd , NORdp
!     .. Arrays in Common ..
      DOUBLE PRECISION R(11)
!     .. Local Scalars ..
      DOUBLE PRECISION rpp , ss
      INTEGER ierpp , nrdp
!     .. Local Arrays ..
      INTEGER itemp(4)
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     .. Common blocks ..
      COMMON /DCHECK/ EPS , R , RP , SVEps , TOL , NORdp , NORd , IERp , IERr
!***FIRST EXECUTABLE STATEMENT  DCMPAR
      Icnt = Icnt + 1
      itemp(1) = 0
      itemp(2) = 0
      itemp(3) = 0
      itemp(4) = 0
      ss = SVEps - EPS
      nrdp = NORdp - NORd
      rpp = RP - R(11)
      ierpp = IERp - IERr
      IF ( ABS(ss)<=TOL.OR.Icnt<=2.OR.Icnt>=6 ) itemp(1) = 1
      IF ( ABS(nrdp)==0 ) itemp(2) = 1
      IF ( Icnt==2 ) itemp(2) = 1
      IF ( ABS(rpp)<=TOL ) itemp(3) = 1
      IF ( ABS(ierpp)==0 ) itemp(4) = 1
!
!     Check to see if all four tests were good.
!     If so, set the test number equal to 1.
!
      Itest(Icnt) = itemp(1)*itemp(2)*itemp(3)*itemp(4)
      END SUBROUTINE DCMPAR
