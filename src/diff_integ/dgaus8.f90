!** DGAUS8
SUBROUTINE DGAUS8(FUN,A,B,Err,Ans,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Integrate a real function of one variable over a finite
  !            interval using an adaptive 8-point Legendre-Gauss
  !            algorithm.  Intended primarily for high accuracy
  !            integration or integration of smooth functions.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A1A1
  !***
  ! **Type:**      DOUBLE PRECISION (GAUS8-S, DGAUS8-D)
  !***
  ! **Keywords:**  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
  !             GAUSS QUADRATURE, NUMERICAL INTEGRATION
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract  *** a DOUBLE PRECISION routine ***
  !        DGAUS8 integrates real functions of one variable over finite
  !        intervals using an adaptive 8-point Legendre-Gauss algorithm.
  !        DGAUS8 is intended primarily for high accuracy integration
  !        or integration of smooth functions.
  !
  !        The maximum number of significant digits obtainable in ANS
  !        is the smaller of 18 and the number of digits carried in
  !        double precision arithmetic.
  !
  !     Description of Arguments
  !
  !        Input--* FUN, A, B, ERR are DOUBLE PRECISION *
  !        FUN - name of external function to be integrated.  This name
  !              must be in an EXTERNAL statement in the calling program.
  !              FUN must be a DOUBLE PRECISION function of one DOUBLE
  !              PRECISION argument.  The value of the argument to FUN
  !              is the variable of integration which ranges from A to B.
  !        A   - lower limit of integration
  !        B   - upper limit of integration (may be less than A)
  !        ERR - is a requested pseudorelative error tolerance.  Normally
  !              pick a value of ABS(ERR) so that DTOL .LT. ABS(ERR) .LE.
  !              1.0D-3 where DTOL is the larger of 1.0D-18 and the
  !              double precision unit roundoff D1MACH(4).  ANS will
  !              normally have no more error than ABS(ERR) times the
  !              integral of the absolute value of FUN(X).  Usually,
  !              smaller values of ERR yield more accuracy and require
  !              more function evaluations.
  !
  !              A negative value for ERR causes an estimate of the
  !              absolute error in ANS to be returned in ERR.  Note that
  !              ERR must be a variable (not a constant) in this case.
  !              Note also that the user must reset the value of ERR
  !              before making any more calls that use the variable ERR.
  !
  !        Output--* ERR,ANS are double precision *
  !        ERR - will be an estimate of the absolute error in ANS if the
  !              input value of ERR was negative.  (ERR is unchanged if
  !              the input value of ERR was non-negative.)  The estimated
  !              error is solely for information to the user and should
  !              not be used as a correction to the computed integral.
  !        ANS - computed value of integral
  !        IERR- a status code
  !            --Normal codes
  !               1 ANS most likely meets requested error tolerance,
  !                 or A=B.
  !              -1 A and B are too nearly equal to allow normal
  !                 integration.  ANS is set to zero.
  !            --Abnormal code
  !               2 ANS probably does not meet requested error tolerance.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, I1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  INTERFACE
    REAL(8) FUNCTION FUN(X)
      REAL(8), INTENT(IN) :: X
    END FUNCTION
  END INTERFACE
  INTEGER Ierr, k, l, lmn, lmx, lr(60), mxl, nbits, nib, nlmx
  INTEGER I1MACH
  REAL(8) :: A, aa(60), ae, anib, Ans, area, B, c, ce, ee, ef, &
    eps, Err, est, gl, glr, gr(60), hh(60), tol, vl(60), vr
  REAL(8) :: D1MACH
  REAL(8), PARAMETER :: x1 = 1.83434642495649805D-01, x2 = 5.25532409916328986D-01, &
    x3 =7.96666477413626740D-01 , x4 = 9.60289856497536232D-01
  REAL(8), PARAMETER :: w1 = 3.62683783378361983D-01, w2 = 3.13706645877887287D-01, &
    w3 = 2.22381034453374471D-01, w4 = 1.01228536290376259D-01
  REAL(8), PARAMETER :: sq2 = 1.41421356D0
  INTEGER, PARAMETER :: nlmn = 1, kmx = 5000, kml = 6
  !* FIRST EXECUTABLE STATEMENT  DGAUS8
  !
  !     Initialize
  !
  k = I1MACH(14)
  anib = D1MACH(5)*k/0.30102000D0
  nbits = INT( anib )
  nlmx = MIN(60,(nbits*5)/8)
  Ans = 0.0D0
  Ierr = 1
  ce = 0.0D0
  IF ( A==B ) THEN
    IF ( Err<0.0D0 ) Err = ce
    RETURN
  ELSE
    lmx = nlmx
    lmn = nlmn
    IF ( B/=0.0D0 ) THEN
      IF ( SIGN(1.0D0,B)*A>0.0D0 ) THEN
        c = ABS(1.0D0-A/B)
        IF ( c<=0.1D0 ) THEN
          IF ( c<=0.0D0 ) THEN
            IF ( Err<0.0D0 ) Err = ce
            RETURN
          ELSE
            anib = 0.5D0 - LOG(c)/0.69314718D0
            nib = INT( anib )
            lmx = MIN(nlmx,nbits-nib-7)
            IF ( lmx<1 ) THEN
              Ierr = -1
              CALL XERMSG('SLATEC','DGAUS8',&
                'A and B are too nearly equal to allow normal integration. $$ANS is set to zero and IERR to -1.',1,-1)
              IF ( Err<0.0D0 ) Err = ce
              RETURN
            ELSE
              lmn = MIN(lmn,lmx)
            END IF
          END IF
        END IF
      END IF
    END IF
    tol = MAX(ABS(Err),2.0D0**(5-nbits))/2.0D0
    IF ( Err==0.0D0 ) tol = SQRT(D1MACH(4))
    eps = tol
    hh(1) = (B-A)/4.0D0
    aa(1) = A
    lr(1) = 1
    l = 1
    est = G8(aa(l)+2.0D0*hh(l),2.0D0*hh(l))
    k = 8
    area = ABS(est)
    ef = 0.5D0
    mxl = 0
  END IF
  100 CONTINUE
  DO
    !
    !     Compute refined estimates, estimate the error, etc.
    !
    gl = G8(aa(l)+hh(l),hh(l))
    gr(l) = G8(aa(l)+3.0D0*hh(l),hh(l))
    k = k + 16
    area = area + (ABS(gl)+ABS(gr(l))-ABS(est))
    !     IF (L .LT .LMN) GO TO 11
    glr = gl + gr(l)
    ee = ABS(est-glr)*ef
    ae = MAX(eps*area,tol*ABS(glr))
    IF ( ee<=ae ) EXIT
    !
    !     Consider the left half of this level
    !
    IF ( k>kmx ) lmx = kml
    IF ( l>=lmx ) THEN
      mxl = 1
      EXIT
    ELSE
      l = l + 1
      eps = eps*0.5D0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5D0
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
    END IF
  END DO
  ce = ce + (est-glr)
  IF ( lr(l)<=0 ) THEN
    !
    !     Proceed to right half at this level
    !
    vl(l) = glr
  ELSE
    !
    !     Return one level
    !
    vr = glr
    DO WHILE ( l>1 )
      l = l - 1
      eps = eps*2.0D0
      ef = ef*sq2
      IF ( lr(l)<=0 ) THEN
        vl(l) = vl(l+1) + vr
        GOTO 200
      ELSE
        vr = vl(l+1) + vr
      END IF
    END DO
    !
    !     Exit
    !
    Ans = vr
    IF ( (mxl/=0).AND.(ABS(ce)>2.0D0*tol*area) ) THEN
      Ierr = 2
      CALL XERMSG('SLATEC','DGAUS8',&
        'ANS is probably insufficiently accurate.',3,1)
    END IF
    IF ( Err<0.0D0 ) Err = ce
    RETURN
  END IF
  200  est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0D0*hh(l)
  GOTO 100
  RETURN
CONTAINS
  REAL(8) FUNCTION G8(x,h)
    REAL(8), INTENT(IN) :: x, h
    G8 = h*((w1*(FUN(x-x1*h)+FUN(x+x1*h))+w2*(FUN(x-x2*h)+FUN(x+x2*h)))&
      +(w3*(FUN(x-x3*h)+FUN(x+x3*h))+w4*(FUN(x-x4*h)+FUN(x+x4*h))))
  END FUNCTION G8
END SUBROUTINE DGAUS8