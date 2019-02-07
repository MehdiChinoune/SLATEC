!*==DQJAC2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DQJAC2
SUBROUTINE DQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
  IMPLICIT NONE
  !*--DQJAC25
  !***BEGIN PROLOGUE  DQJAC2
  !***PURPOSE
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     SUBROUTINE TO EVALUATE THE FULL JACOBIAN FOR TEST PROBLEM USED
  !     IN QUICK CHECK OF DNSQE.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DQJAC2
  INTEGER Iflag , Ldfjac , N
  REAL(8) :: Fjac(Ldfjac,*) , Fvec(*) , X(*)
  !***FIRST EXECUTABLE STATEMENT  DQJAC2
  Fjac(1,1) = -1.0D0
  Fjac(1,2) = 0.0D0
  Fjac(2,1) = -2.0D1*X(1)
  Fjac(2,2) = 1.0D1
END SUBROUTINE DQJAC2
