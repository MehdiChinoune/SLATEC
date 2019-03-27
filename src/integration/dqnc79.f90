!** DQNC79
SUBROUTINE DQNC79(FUN,A,B,Err,Ans,Ierr,K)
  IMPLICIT NONE
  !>
  !***
  !  Integrate a function using a 7-point adaptive Newton-Cotes
  !            quadrature rule.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A1A1
  !***
  ! **Type:**      DOUBLE PRECISION (QNC79-S, DQNC79-D)
  !***
  ! **Keywords:**  ADAPTIVE QUADRATURE, INTEGRATION, NEWTON-COTES
  !***
  ! **Author:**  Kahaner, D. K., (NBS)
  !           Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract  *** a DOUBLE PRECISION routine ***
  !       DQNC79 is a general purpose program for evaluation of
  !       one dimensional integrals of user defined functions.
  !       DQNC79 will pick its own points for evaluation of the
  !       integrand and these will vary from problem to problem.
  !       Thus, DQNC79 is not designed to integrate over data sets.
  !       Moderately smooth integrands will be integrated efficiently
  !       and reliably.  For problems with strong singularities,
  !       oscillations etc., the user may wish to use more sophis-
  !       ticated routines such as those in QUADPACK.  One measure
  !       of the reliability of DQNC79 is the output parameter K,
  !       giving the number of integrand evaluations that were needed.
  !
  !     Description of Arguments
  !
  !     --Input--* FUN, A, B, ERR are DOUBLE PRECISION *
  !       FUN  - name of external function to be integrated.  This name
  !              must be in an EXTERNAL statement in your calling
  !              program.  You must write a Fortran function to evaluate
  !              FUN.  This should be of the form
  !                    DOUBLE PRECISION FUNCTION FUN (X)
  !              C
  !              C     X can vary from A to B
  !              C     FUN(X) should be finite for all X on interval.
  !              C
  !                    FUN = ...
  !                    RETURN
  !                    END
  !       A    - lower limit of integration
  !       B    - upper limit of integration (may be less than A)
  !       ERR  - is a requested error tolerance.  Normally, pick a value
  !              0 .LT. ERR .LT. 1.0D-8.
  !
  !     --Output--
  !       ANS  - computed value of the integral.  Hopefully, ANS is
  !              accurate to within ERR * integral of ABS(FUN(X)).
  !       IERR - a status code
  !            - Normal codes
  !               1  ANS most likely meets requested error tolerance.
  !              -1  A equals B, or A and B are too nearly equal to
  !                  allow normal integration.  ANS is set to zero.
  !            - Abnormal code
  !               2  ANS probably does not meet requested error tolerance.
  !       K    - the number of function evaluations actually used to do
  !              the integration.  A value of K .GT. 1000 indicates a
  !              difficult problem; other programs may be more efficient.
  !              DQNC79 will gracefully give up if K exceeds 2000.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, I1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920218  Code redone to parallel QNC79.  (WRB)
  !   930120  Increase array size 80->99, and KMX 2000->5000 for SUN -r8
  !           wordlength.  (RWC)

  !     .. Scalar Arguments ..
  REAL(8) :: A, Ans, B, Err
  INTEGER Ierr, K
  !     .. Function Arguments ..
  REAL(8) :: FUN
  EXTERNAL :: FUN
  !     .. Local Scalars ..
  REAL(8) :: ae, area, bank, blocal, c, ce, ee, ef, eps, &
    q13, q7, q7l, sq2, test, tol, vr, w1, w2, w3, w4
  INTEGER i, kml, kmx, l, lmn, lmx, nbits, nib, nlmn, nlmx
  !     .. Local Arrays ..
  REAL(8) :: aa(99), f(13), f1(99), f2(99), f3(99), f4(99), &
    f5(99), f6(99), f7(99), hh(99), q7r(99), vl(99)
  INTEGER lr(99)
  !     .. External Functions ..
  REAL(8), EXTERNAL :: D1MACH
  INTEGER, EXTERNAL :: I1MACH
  !     .. External Subroutines ..
  EXTERNAL :: XERMSG
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, LOG, MAX, MIN, SIGN, SQRT
  !     .. Save statement ..
  SAVE nbits, nlmx, sq2, w1, w2, w3, w4
  !     .. Data statements ..
  DATA kml/7/, kmx/5000/, nlmn/2/
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DQNC79
  IF ( first ) THEN
    w1 = 41.0D0/140.0D0
    w2 = 216.0D0/140.0D0
    w3 = 27.0D0/140.0D0
    w4 = 272.0D0/140.0D0
    nbits = INT( D1MACH(5)*I1MACH(14)/0.30102000D0 )
    nlmx = MIN(99,(nbits*4)/5)
    sq2 = SQRT(2.0D0)
    first = .FALSE.
  ENDIF
  Ans = 0.0D0
  Ierr = 1
  ce = 0.0D0
  IF ( A==B ) GOTO 400
  lmx = nlmx
  lmn = nlmn
  IF ( B/=0.0D0 ) THEN
    IF ( SIGN(1.0D0,B)*A>0.0D0 ) THEN
      c = ABS(1.0D0-A/B)
      IF ( c<=0.1D0 ) THEN
        IF ( c<=0.0D0 ) GOTO 400
        nib = INT( 0.5D0 - LOG(c)/LOG(2.0D0) )
        lmx = MIN(nlmx,nbits-nib-4)
        IF ( lmx<2 ) GOTO 400
        lmn = MIN(lmn,lmx)
      ENDIF
    ENDIF
  ENDIF
  tol = MAX(ABS(Err),2.0D0**(5-nbits))
  IF ( Err==0.0D0 ) tol = SQRT(D1MACH(4))
  eps = tol
  hh(1) = (B-A)/12.0D0
  aa(1) = A
  lr(1) = 1
  DO i = 1, 11, 2
    f(i) = FUN(A+(i-1)*hh(1))
  ENDDO
  blocal = B
  f(13) = FUN(blocal)
  K = 7
  l = 1
  area = 0.0D0
  q7 = 0.0D0
  ef = 256.0D0/255.0D0
  bank = 0.0D0
  !
  !     Compute refined estimates, estimate the error, etc.
  !
  100 CONTINUE
  DO i = 2, 12, 2
    f(i) = FUN(aa(l)+(i-1)*hh(l))
  ENDDO
  K = K + 6
  !
  !     Compute left and right half estimates
  !
  q7l = hh(l)*((w1*(f(1)+f(7))+w2*(f(2)+f(6)))+(w3*(f(3)+f(5))+w4*f(4)))
  q7r(l) = hh(l)&
    *((w1*(f(7)+f(13))+w2*(f(8)+f(12)))+(w3*(f(9)+f(11))+w4*f(10)))
  !
  !     Update estimate of integral of absolute value
  !
  area = area + (ABS(q7l)+ABS(q7r(l))-ABS(q7))
  !
  !     Do not bother to test convergence before minimum refinement level
  !
  IF ( l>=lmn ) THEN
    !
    !     Estimate the error in new value for whole interval, Q13
    !
    q13 = q7l + q7r(l)
    ee = ABS(q7-q13)*ef
    !
    !     Compute nominal allowed error
    !
    ae = eps*area
    !
    !     Borrow from bank account, but not too much
    !
    test = MIN(ae+0.8D0*bank,10.0D0*ae)
    !
    !     Don't ask for excessive accuracy
    !
    test = MAX(test,tol*ABS(q13),0.00003D0*tol*area)
    !
    !     Now, did this interval pass or not?
    !
    IF ( ee<=test ) THEN
      !
      !     On good intervals accumulate the theoretical estimate
      !
      ce = ce + (q7-q13)/255.0D0
    ELSE
      !
      !     Consider the left half of next deeper level
      !
      IF ( K>kmx ) lmx = MIN(kml,lmx)
      IF ( l<lmx ) GOTO 200
      !
      !     Have hit maximum refinement level -- penalize the cumulative error
      !
      ce = ce + (q7-q13)
    ENDIF
    !
    !     Update the bank account.  Don't go into debt.
    !
    bank = bank + (ae-ee)
    IF ( bank<0.0D0 ) bank = 0.0D0
    !
    !     Did we just finish a left half or a right half?
    !
    IF ( lr(l)<=0 ) THEN
      !
      !     Proceed to right half at this level
      !
      vl(l) = q13
      GOTO 300
    ELSE
      !
      !     Left and right halves are done, so go back up a level
      !
      vr = q13
      DO WHILE ( l>1 )
        IF ( l<=17 ) ef = ef*sq2
        eps = eps*2.0D0
        l = l - 1
        IF ( lr(l)<=0 ) THEN
          vl(l) = vl(l+1) + vr
          GOTO 300
        ELSE
          vr = vl(l+1) + vr
        ENDIF
      ENDDO
      !
      !     Exit
      !
      Ans = vr
      IF ( ABS(ce)>2.0D0*tol*area ) THEN
        Ierr = 2
        CALL XERMSG('SLATEC','DQNC79',&
          'ANS is probably insufficiently accurate.',2,1)
      ENDIF
      RETURN
    ENDIF
  ENDIF
  200  l = l + 1
  eps = eps*0.5D0
  IF ( l<=17 ) ef = ef/sq2
  hh(l) = hh(l-1)*0.5D0
  lr(l) = -1
  aa(l) = aa(l-1)
  q7 = q7l
  f1(l) = f(7)
  f2(l) = f(8)
  f3(l) = f(9)
  f4(l) = f(10)
  f5(l) = f(11)
  f6(l) = f(12)
  f7(l) = f(13)
  f(13) = f(7)
  f(11) = f(6)
  f(9) = f(5)
  f(7) = f(4)
  f(5) = f(3)
  f(3) = f(2)
  GOTO 100
  300  q7 = q7r(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 12.0D0*hh(l)
  f(1) = f1(l)
  f(3) = f2(l)
  f(5) = f3(l)
  f(7) = f4(l)
  f(9) = f5(l)
  f(11) = f6(l)
  f(13) = f7(l)
  GOTO 100
  400  Ierr = -1
  CALL XERMSG('SLATEC','DQNC79',&
    'A and B are too nearly equal to allow normal integration. $$ANS is set to zero and IERR to -1.',-1,-1)
  RETURN
END SUBROUTINE DQNC79
