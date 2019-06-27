!** D9CHU
REAL(DP) ELEMENTAL FUNCTION D9CHU(A,B,Z)
  !> Evaluate for large Z  Z**A * U(A,B,Z) where U is the logarithmic confluent
  !  hypergeometric function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C11
  !***
  ! **Type:**      DOUBLE PRECISION (R9CHU-S, D9CHU-D)
  !***
  ! **Keywords:**  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate for large Z  Z**A * U(A,B,Z)  where U is the logarithmic
  ! confluent hypergeometric function.  A rational approximation due to Y.
  ! L. Luke is used.  When U is not in the asymptotic region, i.e., when A
  ! or B is large compared with Z, considerable significance loss occurs.
  ! A warning is provided when the computed result is less than half
  ! precision.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: A, B, Z
  INTEGER :: i, j
  REAL(DP) :: aa(4), bb(4), ab, anbn, bp, ct1, ct2, ct3, c2, d1z, g1, g2, g3, sab, x2i1
  REAL(DP), PARAMETER :: eps = 4._DP*D1MACH(4), sqeps = SQRT(D1MACH(4))
  !* FIRST EXECUTABLE STATEMENT  D9CHU
  bp = 1._DP + A - B
  ab = A*bp
  ct2 = 2._DP*(Z-ab)
  sab = A + bp
  !
  bb(1) = 1._DP
  aa(1) = 1._DP
  !
  ct3 = sab + 1._DP + ab
  bb(2) = 1._DP + 2._DP*Z/ct3
  aa(2) = 1._DP + ct2/ct3
  !
  anbn = ct3 + sab + 3._DP
  ct1 = 1._DP + 2._DP*Z/anbn
  bb(3) = 1._DP + 6._DP*ct1*Z/ct3
  aa(3) = 1._DP + 6._DP*ab/anbn + 3._DP*ct1*ct2/ct3
  !
  DO i = 4, 300
    x2i1 = 2*i - 3
    ct1 = x2i1/(x2i1-2._DP)
    anbn = anbn + x2i1 + sab
    ct2 = (x2i1-1._DP)/anbn
    c2 = x2i1*ct2 - 1._DP
    d1z = x2i1*2._DP*Z/anbn
    !
    ct3 = sab*ct2
    g1 = d1z + ct1*(c2+ct3)
    g2 = d1z - c2
    g3 = ct1*(1._DP-ct3-2._DP*ct2)
    !
    bb(4) = g1*bb(3) + g2*bb(2) + g3*bb(1)
    aa(4) = g1*aa(3) + g2*aa(2) + g3*aa(1)
    IF( ABS(aa(4)*bb(1)-aa(1)*bb(4))<eps*ABS(bb(4)*bb(1)) ) GOTO 100
    !
    ! IF OVERFLOWS OR UNDERFLOWS PROVE TO BE A PROBLEM, THE STATEMENTS
    ! BELOW COULD BE ALTERED TO INCORPORATE A DYNAMICALLY ADJUSTED SCALE
    ! FACTOR.
    !
    DO j = 1, 3
      aa(j) = aa(j+1)
      bb(j) = bb(j+1)
    END DO
  END DO
  ERROR STOP 'D9CHU : NO CONVERGENCE IN 300 TERMS'
  !
  100  D9CHU = aa(4)/bb(4)
  !
  ! IF( D9CHU<sqeps .OR. D9CHU>1._DP/sqeps ) CALL XERMSG('D9CHU','ANSWER LT HALF PRECISION',2,1)
  !
END FUNCTION D9CHU