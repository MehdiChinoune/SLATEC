!*==DUIVP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DUIVP
      SUBROUTINE DUIVP(X,Y,Yp)
      IMPLICIT NONE
!*--DUIVP5
!***BEGIN PROLOGUE  DUIVP
!***PURPOSE  Dummy routine for DBVSUP quick check.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (UIVP-S, DUIVP-D)
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   This routine is never called;  it is here to prevent loaders from
!   complaining about undefined externals while testing DBVSUP.
!
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920401  Variables declaration and TYPE sections added.  (WRB)
!***END PROLOGUE  DUIVP
!     .. Scalar Arguments ..
      DOUBLE PRECISION X
!     .. Array Arguments ..
      DOUBLE PRECISION Y(*) , Yp(*)
!***FIRST EXECUTABLE STATEMENT  DUIVP
      STOP
      END SUBROUTINE DUIVP
