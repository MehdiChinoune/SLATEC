!** PPGQ8
SUBROUTINE PPGQ8(FUN,Ldc,C,Xi,Lxi,Kk,Id,A,B,Inppv,Err,Ans,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to PFQAD
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PPGQ8-S, DPPGQ8-D)
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        PPGQ8, a modification of GAUS8, integrates the
  !        product of FUN(X) by the ID-th derivative of a spline
  !        PPVAL(LDC,C,XI,LXI,KK,ID,X,INPPV)  between limits A and B.
  !
  !     Description of arguments
  !
  !        INPUT--
  !        FUN - Name of external function of one argument which
  !              multiplies PPVAL.
  !        LDC - Leading dimension of matrix C, LDC.GE.KK
  !        C   - Matrix of Taylor derivatives of dimension at least
  !              (K,LXI)
  !        XI  - Breakpoint vector of length LXI+1
  !        LXI - Number of polynomial pieces
  !        KK  - Order of the spline, KK.GE.1
  !        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
  !        A   - Lower limit of integral
  !        B   - Upper limit of integral (may be less than A)
  !        INPPV- Initialization parameter for PPVAL
  !        ERR - Is a requested pseudorelative error tolerance.  Normally
  !              pick a value of ABS(ERR).LT.1E-3.  ANS will normally
  !              have no more error than ABS(ERR) times the integral of
  !              the absolute value of FUN(X)*PPVAL(LDC,C,XI,LXI,KK,ID,X,
  !              INPPV).
  !
  !        OUTPUT--
  !        ERR - Will be an estimate of the absolute error in ANS if the
  !              input value of ERR was negative.  (ERR is unchanged if
  !              the input value of ERR was nonnegative.)  The estimated
  !              error is solely for information to the user and should
  !              not be used as a correction to the computed integral.
  !        ANS - Computed value of integral
  !        IERR- A status code
  !            --Normal codes
  !               1 ANS most likely meets requested error tolerance,
  !                 or A=B.
  !              -1 A and B ARE too nearly equal to allow normal
  !                 integration.  ANS is set to zero.
  !            --Abnormal code
  !               2 ANS probably does not meet requested error tolerance.
  !
  !***
  ! **See also:**  PFQAD
  !***
  ! **Routines called:**  I1MACH, PPVAL, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)

  !
  INTERFACE
    REAL FUNCTION FUN(X)
      REAL, INTENT(IN) :: X
    END FUNCTION
  END INTERFACE
  INTEGER Id, Ierr, Inppv, k, Kk, l, Ldc, lmn, lmx, lr(30), Lxi, mxl, nbits, nib, nlmx
  INTEGER I1MACH
  REAL A, aa(30), ae, anib, Ans, area, B, be, C(Ldc,*), cc, ee, ef, eps, &
    Err, est, gl, glr, gr(30), hh(30), tol, vl(30), vr, Xi(*)
  REAL R1MACH
  REAL, PARAMETER :: x1 = 1.83434642495649805E-01, x2 = 5.25532409916328986E-01, &
    x3 =7.96666477413626740E-01 , x4 = 9.60289856497536232E-01
  REAL, PARAMETER ::  w1 =3.62683783378361983E-01 , w2 = 3.13706645877887287E-01, &
    w3 = 2.22381034453374471E-01, w4 = 1.01228536290376259E-01
  REAL, PARAMETER :: sq2 = 1.41421356E0
  INTEGER, PARAMETER :: nlmn = 1, kmx = 5000, kml = 6
  !
  !     INITIALIZE
  !
  !* FIRST EXECUTABLE STATEMENT  PPGQ8
  k = I1MACH(11)
  anib = R1MACH(5)*k/0.30102000E0
  nbits = INT(anib)
  nlmx = (nbits*5)/8
  Ans = 0.0E0
  Ierr = 1
  be = 0.0E0
  IF ( A==B ) THEN
    IF ( Err<0.0E0 ) Err = be
    RETURN
  ELSE
    lmx = nlmx
    lmn = nlmn
    IF ( B/=0.0E0 ) THEN
      IF ( SIGN(1.0E0,B)*A>0.0E0 ) THEN
        cc = ABS(1.0E0-A/B)
        IF ( cc<=0.1E0 ) THEN
          IF ( cc<=0.0E0 ) THEN
            IF ( Err<0.0E0 ) Err = be
            RETURN
          ELSE
            anib = 0.5E0 - LOG(cc)/0.69314718E0
            nib = INT(anib)
            lmx = MIN(nlmx,nbits-nib-7)
            IF ( lmx<1 ) THEN
              Ierr = -1
              CALL XERMSG('SLATEC','PPGQ8',&
                'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.',1,-1)
              IF ( Err<0.0E0 ) Err = be
              RETURN
            ELSE
              lmn = MIN(lmn,lmx)
            END IF
          END IF
        END IF
      END IF
    END IF
    tol = MAX(ABS(Err),2.0E0**(5-nbits))/2.0E0
    IF ( Err==0.0E0 ) tol = SQRT(R1MACH(4))
    eps = tol
    hh(1) = (B-A)/4.0E0
    aa(1) = A
    lr(1) = 1
    l = 1
    est = G8(aa(l)+2.0E0*hh(l),2.0E0*hh(l))
    k = 8
    area = ABS(est)
    ef = 0.5E0
    mxl = 0
  END IF
  100 CONTINUE
  DO
    !
    !     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
    !
    gl = G8(aa(l)+hh(l),hh(l))
    gr(l) = G8(aa(l)+3.0E0*hh(l),hh(l))
    k = k + 16
    area = area + (ABS(gl)+ABS(gr(l))-ABS(est))
    glr = gl + gr(l)
    ee = ABS(est-glr)*ef
    ae = MAX(eps*area,tol*ABS(glr))
    IF ( ee<=ae ) EXIT
    !
    !     CONSIDER THE LEFT HALF OF THIS LEVEL
    !
    IF ( k>kmx ) lmx = kml
    IF ( l>=lmx ) THEN
      mxl = 1
      EXIT
    ELSE
      l = l + 1
      eps = eps*0.5E0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5E0
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
    END IF
  END DO
  be = be + (est-glr)
  IF ( lr(l)<=0 ) THEN
    !
    !     PROCEED TO RIGHT HALF AT THIS LEVEL
    !
    vl(l) = glr
  ELSE
    !
    !     RETURN ONE LEVEL
    !
    vr = glr
    DO WHILE ( l>1 )
      l = l - 1
      eps = eps*2.0E0
      ef = ef*sq2
      IF ( lr(l)<=0 ) THEN
        vl(l) = vl(l+1) + vr
        GOTO 200
      ELSE
        vr = vl(l+1) + vr
      END IF
    END DO
    !
    !      EXIT
    !
    Ans = vr
    IF ( (mxl/=0).AND.(ABS(be)>2.0E0*tol*area) ) THEN
      Ierr = 2
      CALL XERMSG('SLATEC','PPGQ8',&
        'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.',3,1)
    END IF
    IF ( Err<0.0E0 ) Err = be
    RETURN
  END IF
  200  est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0E0*hh(l)
  GOTO 100
  RETURN
CONTAINS
  REAL FUNCTION G8(x,h)
    REAL, INTENT(IN) :: x, h
    REAL, EXTERNAL :: PPVAL
    G8 = h*((w1*(FUN(x-x1*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x-x1*h,Inppv)+FUN(&
      x+x1*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x+x1*h,Inppv))&
      +w2*(FUN(x-x2*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x-x2*h,Inppv)&
      +FUN(x+x2*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x+x2*h,Inppv)))&
      +(w3*(FUN(x-x3*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x-x3*h,Inppv)&
      +FUN(x+x3*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x+x3*h,Inppv))&
      +w4*(FUN(x-x4*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x-x4*h,Inppv)&
      +FUN(x+x4*h)*PPVAL(Ldc,C,Xi,Lxi,Kk,Id,x+x4*h,Inppv))))
  END FUNCTION G8
END SUBROUTINE PPGQ8