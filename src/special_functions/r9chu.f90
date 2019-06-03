!** R9CHU
REAL FUNCTION R9CHU(A,B,Z)
  !>
  !  Evaluate for large Z  Z**A * U(A,B,Z) where U is the
  !            logarithmic confluent hypergeometric function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C11
  !***
  ! **Type:**      SINGLE PRECISION (R9CHU-S, D9CHU-D)
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
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL :: A, B, Z
  INTEGER:: i, j
  REAL :: aa(4), ab, anbn, bb(4), bp, c2, ct1, ct2, ct3, d1z, g1, g2, g3, sab, x2i1
  REAL, PARAMETER :: eps = 4.0*R1MACH(4), sqeps = SQRT(R1MACH(4))
  !* FIRST EXECUTABLE STATEMENT  R9CHU
  !
  bp = 1.0 + A - B
  ab = A*bp
  ct2 = 2.0*(Z-ab)
  sab = A + bp
  !
  bb(1) = 1.0
  aa(1) = 1.0
  !
  ct3 = sab + 1.0 + ab
  bb(2) = 1.0 + 2.0*Z/ct3
  aa(2) = 1.0 + ct2/ct3
  !
  anbn = ct3 + sab + 3.0
  ct1 = 1.0 + 2.0*Z/anbn
  bb(3) = 1.0 + 6.0*ct1*Z/ct3
  aa(3) = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3
  !
  DO i = 4, 300
    x2i1 = 2*i - 3
    ct1 = x2i1/(x2i1-2.0)
    anbn = anbn + x2i1 + sab
    ct2 = (x2i1-1.0)/anbn
    c2 = x2i1*ct2 - 1.0
    d1z = x2i1*2.0*Z/anbn
    !
    ct3 = sab*ct2
    g1 = d1z + ct1*(c2+ct3)
    g2 = d1z - c2
    g3 = ct1*(1.0-ct3-2.0*ct2)
    !
    bb(4) = g1*bb(3) + g2*bb(2) + g3*bb(1)
    aa(4) = g1*aa(3) + g2*aa(2) + g3*aa(1)
    IF ( ABS(aa(4)*bb(1)-aa(1)*bb(4))<eps*ABS(bb(4)*bb(1)) ) GOTO 100
    !
    ! IF OVERFLOWS OR UNDERFLOWS PROVE TO BE A PROBLEM, THE STATEMENTS
    ! BELOW COULD BE ALTERED TO INCORPORATE A DYNAMICALLY ADJUSTED SCALE
    ! FACTOR.
    !
    DO j = 1, 3
      bb(j) = bb(j+1)
      aa(j) = aa(j+1)
    END DO
  END DO
  CALL XERMSG('R9CHU','NO CONVERGENCE IN 300 TERMS',1,2)
  !
  100  R9CHU = aa(4)/bb(4)
  !
  IF ( R9CHU<sqeps.OR.R9CHU>1.0/sqeps )&
    CALL XERMSG('R9CHU','ANSWER LESS THAN HALF PRECISION',2,1)
  !
END FUNCTION R9CHU
