!** DBSQAD
SUBROUTINE DBSQAD(T,Bcoef,N,K,X1,X2,Bquad,Work)
  !>
  !***
  !  Compute the integral of a K-th order B-spline using the
  !            B-representation.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A2A1, E3, K6
  !***
  ! **Type:**      DOUBLE PRECISION (BSQAD-S, DBSQAD-D)
  !***
  ! **Keywords:**  INTEGRAL OF B-SPLINES, QUADRATURE
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract    **** a double precision routine ****
  !
  !         DBSQAD computes the integral on (X1,X2) of a K-th order
  !         B-spline using the B-representation (T,BCOEF,N,K).  Orders
  !         K as high as 20 are permitted by applying a 2, 6, or 10
  !         point Gauss formula on subintervals of (X1,X2) which are
  !         formed by included (distinct) knots.
  !
  !         If orders K greater than 20 are needed, use DBFQAD with
  !         F(X) = 1.
  !
  !         The maximum number of significant digits obtainable in
  !         DBSQAD is the smaller of 18 and the number of digits
  !         carried in double precision arithmetic.
  !
  !     Description of Arguments
  !         Input      T,BCOEF,X1,X2 are double precision
  !           T      - knot array of length N+K
  !           BCOEF  - B-spline coefficient array of length N
  !           N      - length of coefficient array
  !           K      - order of B-spline, 1 .LE. K .LE. 20
  !           X1,X2  - end points of quadrature interval in
  !                    T(K) .LE. X .LE. T(N+1)
  !
  !         Output     BQUAD,WORK are double precision
  !           BQUAD  - integral of the B-spline over (X1,X2)
  !           WORK   - work vector of length 3*K
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***
  ! **References:**  D. E. Amos, Quadrature subroutines for splines and
  !                 B-splines, Report SAND79-1825, Sandia Laboratories,
  !                 December 1979.
  !***
  ! **Routines called:**  DBVALU, DINTRV, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, il1, il2, ilo, inbv, jf, K, left, m, mf, mflag, N, npk, np1
  REAL(8) :: a, aa, b, bb, Bcoef(*), bma, bpa, Bquad, c1, gx, q, summ(5), T(*), &
    ta, tb, Work(*), X1, X2, y1, y2
  !
  REAL(8), PARAMETER :: gpts(9) = [ 5.77350269189625764D-01, 2.38619186083196909D-01, &
    6.61209386466264514D-01, 9.32469514203152028D-01, 1.48874338981631211D-01, &
    4.33395394129247191D-01, 6.79409568299024406D-01, 8.65063366688984511D-01, &
    9.73906528517171720D-01 ]
  REAL(8), PARAMETER :: gwts(9) = [ 1.00000000000000000D+00, 4.67913934572691047D-01, &
    3.60761573048138608D-01, 1.71324492379170345D-01, 2.95524224714752870D-01, &
    2.69266719309996355D-01, 2.19086362515982044D-01, 1.49451349150580593D-01, &
    6.66713443086881376D-02 ]
  !
  !* FIRST EXECUTABLE STATEMENT  DBSQAD
  Bquad = 0.0D0
  IF ( K<1.OR.K>20 ) THEN
    CALL XERMSG('SLATEC','DBSQAD','K DOES NOT SATISFY 1.LE.K.LE.20',2,1)
    RETURN
  ELSEIF ( N<K ) THEN
    CALL XERMSG('SLATEC','DBSQAD','N DOES NOT SATISFY N.GE.K',2,1)
    RETURN
  ELSE
    aa = MIN(X1,X2)
    bb = MAX(X1,X2)
    IF ( aa>=T(K) ) THEN
      np1 = N + 1
      IF ( bb<=T(np1) ) THEN
        IF ( aa==bb ) RETURN
        npk = N + K
        !     SELECTION OF 2, 6, OR 10 POINT GAUSS FORMULA
        jf = 0
        mf = 1
        IF ( K>4 ) THEN
          jf = 1
          mf = 3
          IF ( K>12 ) THEN
            jf = 4
            mf = 5
          END IF
        END IF
        !
        DO i = 1, mf
          summ(i) = 0.0D0
        END DO
        ilo = 1
        inbv = 1
        CALL DINTRV(T,npk,aa,ilo,il1,mflag)
        CALL DINTRV(T,npk,bb,ilo,il2,mflag)
        IF ( il2>=np1 ) il2 = N
        DO left = il1, il2
          ta = T(left)
          tb = T(left+1)
          IF ( ta/=tb ) THEN
            a = MAX(aa,ta)
            b = MIN(bb,tb)
            bma = 0.5D0*(b-a)
            bpa = 0.5D0*(b+a)
            DO m = 1, mf
              c1 = bma*gpts(jf+m)
              gx = -c1 + bpa
              y2 = DBVALU(T,Bcoef,N,K,0,gx,inbv,Work)
              gx = c1 + bpa
              y1 = DBVALU(T,Bcoef,N,K,0,gx,inbv,Work)
              summ(m) = summ(m) + (y1+y2)*bma
            END DO
          END IF
        END DO
        q = 0.0D0
        DO m = 1, mf
          q = q + gwts(jf+m)*summ(m)
        END DO
        IF ( X1>X2 ) q = -q
        Bquad = q
        RETURN
      END IF
    END IF
  END IF
  !
  !
  CALL XERMSG('SLATEC','DBSQAD',&
    'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)',2,1)
  RETURN
END SUBROUTINE DBSQAD
