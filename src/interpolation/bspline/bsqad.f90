!** BSQAD
SUBROUTINE BSQAD(T,Bcoef,N,K,X1,X2,Bquad,Work)
  !>
  !  Compute the integral of a K-th order B-spline using the
  !            B-representation.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2A2A1, E3, K6
  !***
  ! **Type:**      SINGLE PRECISION (BSQAD-S, DBSQAD-D)
  !***
  ! **Keywords:**  INTEGRAL OF B-SPLINES, QUADRATURE
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BSQAD computes the integral on (X1,X2) of a K-th order
  !         B-spline using the B-representation (T,BCOEF,N,K).  Orders
  !         K as high as 20 are permitted by applying a 2, 6, or 10
  !         point Gauss formula on subintervals of (X1,X2) which are
  !         formed by included (distinct) knots.
  !
  !         If orders K greater than 20 are needed, use BFQAD with
  !         F(X) = 1.
  !
  !     Description of Arguments
  !         Input
  !           T      - knot array of length N+K
  !           BCOEF  - B-spline coefficient array of length N
  !           N      - length of coefficient array
  !           K      - order of B-spline, 1 .LE. K .LE. 20
  !           X1,X2  - end points of quadrature interval in
  !                    T(K) .LE. X .LE. T(N+1)
  !
  !         Output
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
  ! **Routines called:**  BVALU, INTRV, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : XERMSG
  !
  INTEGER :: K, N
  REAL(SP) :: Bcoef(N), Bquad, T(N+K), Work(3*K), X1, X2, y1, y2
  INTEGER :: i, il1, il2, ilo, inbv, jf, left, m, mf, mflag, npk, np1
  REAL(SP) :: a, aa, b, bb, bma, bpa, c1, gx, q, summ(5), ta, tb

  !
  REAL(SP), PARAMETER :: gpts(9) = [ 5.77350269189625764E-01, 2.38619186083196909E-01, &
    6.61209386466264514E-01, 9.32469514203152028E-01, 1.48874338981631211E-01, &
    4.33395394129247191E-01, 6.79409568299024406E-01, 8.65063366688984511E-01, &
    9.73906528517171720E-01 ]
  REAL(SP), PARAMETER :: gwts(9) = [ 1.00000000000000000E+00, 4.67913934572691047E-01, &
    3.60761573048138608E-01, 1.71324492379170345E-01, 2.95524224714752870E-01, &
    2.69266719309996355E-01, 2.19086362515982044E-01, 1.49451349150580593E-01, &
    6.66713443086881376E-02 ]
  !
  !* FIRST EXECUTABLE STATEMENT  BSQAD
  Bquad = 0.0E0
  IF ( K<1.OR.K>20 ) THEN
    CALL XERMSG('BSQAD','K DOES NOT SATISFY 1.LE.K.LE.20',2,1)
    RETURN
  ELSEIF ( N<K ) THEN
    CALL XERMSG('BSQAD','N DOES NOT SATISFY N.GE.K',2,1)
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
          summ(i) = 0.0E0
        END DO
        ilo = 1
        inbv = 1
        CALL INTRV(T,npk,aa,ilo,il1,mflag)
        CALL INTRV(T,npk,bb,ilo,il2,mflag)
        IF ( il2>=np1 ) il2 = N
        DO left = il1, il2
          ta = T(left)
          tb = T(left+1)
          IF ( ta/=tb ) THEN
            a = MAX(aa,ta)
            b = MIN(bb,tb)
            bma = 0.5E0*(b-a)
            bpa = 0.5E0*(b+a)
            DO m = 1, mf
              c1 = bma*gpts(jf+m)
              gx = -c1 + bpa
              y2 = BVALU(T,Bcoef,N,K,0,gx,inbv,Work)
              gx = c1 + bpa
              y1 = BVALU(T,Bcoef,N,K,0,gx,inbv,Work)
              summ(m) = summ(m) + (y1+y2)*bma
            END DO
          END IF
        END DO
        q = 0.0E0
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
  CALL XERMSG('BSQAD',&
    'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)',2,1)
  RETURN
END SUBROUTINE BSQAD
