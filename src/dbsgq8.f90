!*==DBSGQ8.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBSGQ8
SUBROUTINE DBSGQ8(FUN,Xt,Bc,N,Kk,Id,A,B,Inbv,Err,Ans,Ierr,Work)
  IMPLICIT NONE
  !*--DBSGQ85
  !***BEGIN PROLOGUE  DBSGQ8
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBFQAD
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BSGQ8-S, DBSGQ8-D)
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract    **** A DOUBLE PRECISION routine ****
  !
  !        DBSGQ8, a modification of GAUS8, integrates the
  !        product of FUN(X) by the ID-th derivative of a spline
  !        DBVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
  !
  !     Description of Arguments
  !
  !        INPUT-- FUN,XT,BC,A,B,ERR are DOUBLE PRECISION
  !        FUN - Name of external function of one argument which
  !              multiplies DBVALU.
  !        XT  - Knot array for DBVALU
  !        BC  - B-coefficient array for DBVALU
  !        N   - Number of B-coefficients for DBVALU
  !        KK  - Order of the spline, KK.GE.1
  !        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1
  !        A   - Lower limit of integral
  !        B   - Upper limit of integral (may be less than A)
  !        INBV- Initialization parameter for DBVALU
  !        ERR - Is a requested pseudorelative error tolerance.  Normally
  !              pick a value of ABS(ERR).LT.1D-3.  ANS will normally
  !              have no more error than ABS(ERR) times the integral of
  !              the absolute value of FUN(X)*DBVALU(XT,BC,N,KK,X,ID,
  !              INBV,WORK).
  !
  !
  !        OUTPUT-- ERR,ANS,WORK are DOUBLE PRECISION
  !        ERR - Will be an estimate of the absolute error in ANS if the
  !              input value of ERR was negative.  (ERR is unchanged if
  !              the input value of ERR was nonnegative.)  The estimated
  !              error is solely for information to the user and should
  !              not be used as a correction to the computed integral.
  !        ANS - Computed value of integral
  !        IERR- A status code
  !            --Normal Codes
  !               1 ANS most likely meets requested error tolerance,
  !                 or A=B.
  !              -1 A and B are too nearly equal to allow normal
  !                 integration.  ANS is set to zero.
  !            --Abnormal Code
  !               2 ANS probably does not meet requested error tolerance.
  !        WORK- Work vector of length 3*K for DBVALU
  !
  !***SEE ALSO  DBFQAD
  !***ROUTINES CALLED  D1MACH, DBVALU, I1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  !***END PROLOGUE  DBSGQ8
  !
  INTERFACE
    REAL(8) FUNCTION FUN(X)
      REAL(8), INTENT(IN) :: X
    END FUNCTION
  END INTERFACE
  INTEGER Id , Ierr , Inbv , k , Kk , kml , kmx , l , lmn , lmx , lr , mxl ,&
    N , nbits , nib , nlmn , nlmx
  INTEGER I1MACH
  REAL(8) :: A , aa , ae , anib , Ans , area , B , Bc , c , ce , ee ,&
    ef , eps , Err , est , gl , glr , gr , hh , sq2 , tol ,&
    vl , vr , Work , w1 , w2 , w3 , w4 , Xt , x1 , x2 , x3 ,&
    x4 , x , h
  REAL(8) :: D1MACH , DBVALU
  DIMENSION Xt(*) , Bc(*) , Work(*)
  DIMENSION aa(60) , hh(60) , lr(60) , vl(60) , gr(60)
  SAVE x1 , x2 , x3 , x4 , w1 , w2 , w3 , w4 , sq2 , nlmn , kmx , kml
  DATA x1 , x2 , x3 , x4/1.83434642495649805D-01 , 5.25532409916328986D-01 ,&
    7.96666477413626740D-01 , 9.60289856497536232D-01/
  DATA w1 , w2 , w3 , w4/3.62683783378361983D-01 , 3.13706645877887287D-01 ,&
    2.22381034453374471D-01 , 1.01228536290376259D-01/
  DATA sq2/1.41421356D0/
  DATA nlmn/1/ , kmx/5000/ , kml/6/
  !
  !     INITIALIZE
  !
  !***FIRST EXECUTABLE STATEMENT  DBSGQ8
  k = I1MACH(14)
  anib = D1MACH(5)*k/0.30102000D0
  nbits = INT(anib)
  nlmx = MIN((nbits*5)/8,60)
  Ans = 0.0D0
  Ierr = 1
  ce = 0.0D0
  IF ( A==B ) THEN
    IF ( Err<0.0D0 ) Err = ce
    GOTO 99999
  ELSE
    lmx = nlmx
    lmn = nlmn
    IF ( B/=0.0D0 ) THEN
      IF ( SIGN(1.0D0,B)*A>0.0D0 ) THEN
        c = ABS(1.0D0-A/B)
        IF ( c<=0.1D0 ) THEN
          IF ( c<=0.0D0 ) THEN
            IF ( Err<0.0D0 ) Err = ce
            GOTO 99999
          ELSE
            anib = 0.5D0 - LOG(c)/0.69314718D0
            nib = INT(anib)
            lmx = MIN(nlmx,nbits-nib-7)
            IF ( lmx<1 ) THEN
              Ierr = -1
              CALL XERMSG('SLATEC','DBSGQ8',&
                'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. '&
                //' ANS IS SET TO ZERO AND IERR TO -1.',1,-1)
              IF ( Err<0.0D0 ) Err = ce
              GOTO 99999
            ELSE
              lmn = MIN(lmn,lmx)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
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
  ENDIF
  100  DO
  !
  !     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
  !
  gl = G8(aa(l)+hh(l),hh(l))
  gr(l) = G8(aa(l)+3.0D0*hh(l),hh(l))
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
    eps = eps*0.5D0
    ef = ef/sq2
    hh(l) = hh(l-1)*0.5D0
    lr(l) = -1
    aa(l) = aa(l-1)
    est = gl
  ENDIF
ENDDO
ce = ce + (est-glr)
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
    eps = eps*2.0D0
    ef = ef*sq2
    IF ( lr(l)<=0 ) THEN
      vl(l) = vl(l+1) + vr
      GOTO 200
    ELSE
      vr = vl(l+1) + vr
    ENDIF
  ENDDO
  !
  !      EXIT
  !
  Ans = vr
  IF ( (mxl/=0).AND.(ABS(ce)>2.0D0*tol*area) ) THEN
    Ierr = 2
    CALL XERMSG('SLATEC','DBSGQ8',&
      'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.',3,1)
  ENDIF
  IF ( Err<0.0D0 ) Err = ce
  GOTO 99999
ENDIF
200  est = gr(l-1)
lr(l) = 1
aa(l) = aa(l) + 4.0D0*hh(l)
GOTO 100
99999 RETURN
CONTAINS
REAL(8) FUNCTION G8(x,h)
  REAL(8), INTENT(IN) :: x, h
  REAL(8), EXTERNAL :: DBVALU
  G8 = h*((w1*(FUN(x-x1*h)*DBVALU(Xt,Bc,N,Kk,Id,x-x1*h,Inbv,Work)+FUN(x+&
    x1*h)*DBVALU(Xt,Bc,N,Kk,Id,x+x1*h,Inbv,Work))&
    +w2*(FUN(x-x2*h)*DBVALU(Xt,Bc,N,Kk,Id,x-x2*h,Inbv,Work)&
    +FUN(x+x2*h)*DBVALU(Xt,Bc,N,Kk,Id,x+x2*h,Inbv,Work)))&
    +(w3*(FUN(x-x3*h)*DBVALU(Xt,Bc,N,Kk,Id,x-x3*h,Inbv,Work)&
    +FUN(x+x3*h)*DBVALU(Xt,Bc,N,Kk,Id,x+x3*h,Inbv,Work))&
    +w4*(FUN(x-x4*h)*DBVALU(Xt,Bc,N,Kk,Id,x-x4*h,Inbv,Work)&
    +FUN(x+x4*h)*DBVALU(Xt,Bc,N,Kk,Id,x+x4*h,Inbv,Work))))
END FUNCTION G8
END SUBROUTINE DBSGQ8
