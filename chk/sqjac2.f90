!*==SQJAC2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SQJAC2
      SUBROUTINE SQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
      IMPLICIT NONE
!*--SQJAC25
!*** Start of declarations inserted by SPAG
      REAL Fjac , Fvec , X
      INTEGER Iflag , Ldfjac , N
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SQJAC2
!***PURPOSE  Evaluate full Jacobian for SNSQE test.
!***LIBRARY   SLATEC
!***KEYWORDS  QUICK CHECK
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     SUBROUTINE TO EVALUATE THE FULL JACOBIAN FOR TEST PROBLEM USED
!     IN QUICK CHECK OF SNSQE.
!
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  SQJAC2
      DIMENSION X(*) , Fvec(*) , Fjac(Ldfjac,*)
!***FIRST EXECUTABLE STATEMENT  SQJAC2
      Fjac(1,1) = -1.E0
      Fjac(1,2) = 0.E0
      Fjac(2,1) = -2.E1*X(1)
      Fjac(2,2) = 1.E1
      END SUBROUTINE SQJAC2
