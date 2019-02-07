!*==DBSPEV.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBSPEV
SUBROUTINE DBSPEV(T,Ad,N,K,Nderiv,X,Inev,Svalue,Work)
  IMPLICIT NONE
  !*--DBSPEV5
  !***BEGIN PROLOGUE  DBSPEV
  !***PURPOSE  Calculate the value of the spline and its derivatives from
  !            the B-representation.
  !***LIBRARY   SLATEC
  !***CATEGORY  E3, K6
  !***TYPE      DOUBLE PRECISION (BSPEV-S, DBSPEV-D)
  !***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract    **** a double precision routine ****
  !         DBSPEV is the BSPLEV routine of the reference.
  !
  !         DBSPEV calculates the value of the spline and its derivatives
  !         at X from the B-representation (T,A,N,K) and returns them in
  !         SVALUE(I),I=1,NDERIV, T(K) .LE. X .LE. T(N+1).  AD(I) can be
  !         the B-spline coefficients A(I), I=1,N) if NDERIV=1.  Otherwise
  !         AD must be computed before hand by a call to DBSPDR (T,A,N,K,
  !         NDERIV,AD).  If X=T(I),I=K,N), right limiting values are
  !         obtained.
  !
  !         To compute left derivatives or left limiting values at a
  !         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
  !
  !         DBSPEV calls DINTRV, DBSPVN
  !
  !     Description of Arguments
  !
  !         Input      T,AD,X, are double precision
  !          T       - knot vector of length N+K
  !          AD      - vector of length (2*N-NDERIV+1)*NDERIV/2 containing
  !                    the difference table from DBSPDR.
  !          N       - number of B-spline coefficients
  !                    N = sum of knot multiplicities-K
  !          K       - order of the B-spline, K .GE. 1
  !          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K.
  !                    NDERIV=1 gives the zero-th derivative =
  !                    function value
  !          X       - argument, T(K) .LE. X .LE. T(N+1)
  !          INEV    - an initialization parameter which must be set
  !                    to 1 the first time DBSPEV is called.
  !
  !         Output     SVALUE,WORK are double precision
  !          INEV    - INEV contains information for efficient process-
  !                    ing after the initial call and INEV must not
  !                    be changed by the user.  Distinct splines require
  !                    distinct INEV parameters.
  !          SVALUE  - vector of length NDERIV containing the spline
  !                    value in SVALUE(1) and the NDERIV-1 derivatives
  !                    in the remaining components.
  !          WORK    - work vector of length 3*K
  !
  !     Error Conditions
  !         Improper input is a fatal error.
  !
  !***REFERENCES  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***ROUTINES CALLED  DBSPVN, DINTRV, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DBSPEV
  !
  INTEGER i , id , Inev , iwork , jj , K , kp1 , kp1mn , l , left , ll , &
    mflag , N , Nderiv
  DOUBLE PRECISION Ad , Svalue , sum , T , Work , X
  !     DIMENSION T(N+K)
  DIMENSION T(*) , Ad(*) , Svalue(*) , Work(*)
  !***FIRST EXECUTABLE STATEMENT  DBSPEV
  IF ( K<1 ) THEN
    !
    !
    CALL XERMSG('SLATEC','DBSPEV','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( N<K ) THEN
    CALL XERMSG('SLATEC','DBSPEV','N DOES NOT SATISFY N.GE.K',2,1)
    RETURN
  ELSEIF ( Nderiv<1.OR.Nderiv>K ) THEN
    CALL XERMSG('SLATEC','DBSPEV','NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K'&
      ,2,1)
    RETURN
  ELSE
    id = Nderiv
    CALL DINTRV(T,N+1,X,Inev,i,mflag)
    IF ( X>=T(K) ) THEN
      IF ( mflag/=0 ) THEN
        IF ( X>T(i) ) GOTO 100
        DO WHILE ( i/=K )
          i = i - 1
          IF ( X/=T(i) ) GOTO 20
        ENDDO
        CALL XERMSG('SLATEC','DBSPEV',&
          'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)',2,1)
        GOTO 99999
      ENDIF
      !
      ! *I* HAS BEEN FOUND IN (K,N) SO THAT T(I) .LE. X .LT. T(I+1)
      !     (OR .LE. T(I+1), IF T(I) .LT. T(I+1) = T(N+1) ).
      20       kp1mn = K + 1 - id
      kp1 = K + 1
      CALL DBSPVN(T,kp1mn,K,1,X,i,Work(1),Work(kp1),iwork)
      jj = (N+N-id+2)*(id-1)/2
      DO
        !     ADIF(LEFTPL,ID) = AD(LEFTPL-ID+1 + (2*N-ID+2)*(ID-1)/2)
        !     LEFTPL = LEFT + L
        left = i - kp1mn
        sum = 0.0D0
        ll = left + jj + 2 - id
        DO l = 1 , kp1mn
          sum = sum + Work(l)*Ad(ll)
          ll = ll + 1
        ENDDO
        Svalue(id) = sum
        id = id - 1
        IF ( id==0 ) THEN
          !
          RETURN
        ELSE
          jj = jj - (N-id+1)
          kp1mn = kp1mn + 1
          CALL DBSPVN(T,kp1mn,K,2,X,i,Work(1),Work(kp1),iwork)
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  100  CALL XERMSG('SLATEC','DBSPEV','X IS NOT IN T(K).LE.X.LE.T(N+1)',2,1)
  RETURN
  99999 END SUBROUTINE DBSPEV
