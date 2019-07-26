!** BINT4
PURE SUBROUTINE BINT4(X,Y,Ndata,Ibcl,Ibcr,Fbcl,Fbcr,Kntopt,T,Bcoef,N,K,W)
  !> Compute the B-representation of a cubic spline which interpolates given data.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E1A
  !***
  ! **Type:**      SINGLE PRECISION (BINT4-S, DBINT4-D)
  !***
  ! **Keywords:**  B-SPLINE, CUBIC SPLINES, DATA FITTING, INTERPOLATION
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BINT4 computes the B representation (T,BCOEF,N,K) of a
  !         cubic spline (K=4) which interpolates data (X(I)),Y(I))),
  !         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the
  !         specification of the spline first or second derivative at
  !         both X(1) and X(NDATA).  When this data is not specified
  !         by the problem, it is common practice to use a natural
  !         spline by setting second derivatives at X(1) and X(NDATA)
  !         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined on
  !         T(4) <= X <= T(N+1) with (ordered) interior knots at X(I))
  !         values where N=NDATA+2.  The knots T(1), T(2), T(3) lie to
  !         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4)
  !         lie to the right of T(N+1)=X(NDATA) in increasing order.  If
  !         no extrapolation outside (X(1),X(NDATA)) is anticipated, the
  !         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)=
  !         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2
  !         selects a knot placement for T(1), T(2), T(3) to make the
  !         first 7 knots symmetric about T(4)=X(1) and similarly for
  !         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3
  !         allows the user to make his own selection, in increasing
  !         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2),
  !         T(N+3), T(N+4) to the right of X(NDATA) in the work array
  !         W(1) through W(6).  In any case, the interpolation on
  !         T(4) <= X <= T(N+1) by using function BVALU is unique
  !         for given boundary conditions.
  !
  !     Description of Arguments
  !         Input
  !           X      - X vector of abscissae of length NDATA, distinct
  !                    and in increasing order
  !           Y      - Y vector of ordinates of length NDATA
  !           NDATA  - number of data points, NDATA >= 2
  !           IBCL   - selection parameter for left boundary condition
  !                    IBCL = 1 constrain the first derivative at
  !                             X(1) to FBCL
  !                         = 2 constrain the second derivative at
  !                             X(1) to FBCL
  !           IBCR   - selection parameter for right boundary condition
  !                    IBCR = 1 constrain first derivative at
  !                             X(NDATA) to FBCR
  !                    IBCR = 2 constrain second derivative at
  !                             X(NDATA) to FBCR
  !           FBCL   - left boundary values governed by IBCL
  !           FBCR   - right boundary values governed by IBCR
  !           KNTOPT - knot selection parameter
  !                    KNTOPT = 1 sets knot multiplicity at T(4) and
  !                               T(N+1) to 4
  !                           = 2 sets a symmetric placement of knots
  !                               about T(4) and T(N+1)
  !                           = 3 sets TNP)=WNP) and T(N+1+I)=w(3+I),I=1,3
  !                               where WNP),I=1,6 is supplied by the user
  !           W      - work array of dimension at least 5*(NDATA+2)
  !                    if KNTOPT=3, then W(1),W(2),W(3) are knot values to
  !                    the left of X(1) and W(4),W(5),W(6) are knot
  !                    values to the right of X(NDATA) in increasing
  !                    order to be supplied by the user
  !
  !         Output
  !           T      - knot array of length N+4
  !           BCOEF  - B-spline coefficient array of length N
  !           N      - number of coefficients, N=NDATA+2
  !           K      - order of spline, K=4
  !
  !     Error Conditions
  !         Improper  input is a fatal error
  !         Singular system of equations is a fatal error
  !
  !***
  ! **References:**  D. E. Amos, Computation with splines and B-splines,
  !                 Report SAND78-1968, Sandia Laboratories, March 1979.
  !               Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !               Carl de Boor, A Practical Guide to Splines, Applied
  !                 Mathematics Series 27, Springer-Verlag, New York,
  !                 1978.
  !***
  ! **Routines called:**  BNFAC, BNSLV, BSPVD, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : eps_sp
  !
  INTEGER, INTENT(IN) :: Ibcl, Ibcr, Kntopt, Ndata
  INTEGER, INTENT(OUT) :: K, N
  REAL(SP), INTENT(IN) :: Fbcl, Fbcr, X(Ndata), Y(Ndata)
  REAL(SP), INTENT(OUT) :: Bcoef(Ndata+2), T(Ndata+6), W(5,Ndata+2)
  !
  INTEGER :: i, iflag, ilb, ileft, it, iub, iw, iwp, j, jw, ndm, np, nwrow
  REAL(SP) :: tol, txn, tx1, vnikx(4,4), wdtol, work(15), xl
  !* FIRST EXECUTABLE STATEMENT  BINT4
  wdtol = eps_sp
  tol = SQRT(wdtol)
  IF( Ndata<2 ) THEN
    ERROR STOP 'BINT4 : NDATA IS LESS THAN 2'
  ELSE
    ndm = Ndata - 1
    DO i = 1, ndm
      IF( X(i)>=X(i+1) ) GOTO 50
    END DO
    IF( Ibcl<1 .OR. Ibcl>2 ) THEN
      ERROR STOP 'BINT4 : IBCL IS NOT 1 OR 2'
    ELSEIF( Ibcr<1 .OR. Ibcr>2 ) THEN
      ERROR STOP 'BINT4 : IBCR IS NOT 1 OR 2'
    ELSEIF( Kntopt<1 .OR. Kntopt>3 ) THEN
      ERROR STOP 'BINT4 : KNTOPT IS NOT 1, 2, OR 3'
    ELSE
      K = 4
      N = Ndata + 2
      np = N + 1
      DO i = 1, Ndata
        T(i+3) = X(i)
      END DO
      SELECT CASE (Kntopt)
        CASE (2)
          !     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS
          IF( Ndata>3 ) THEN
            tx1 = X(1) + X(1)
            txn = X(Ndata) + X(Ndata)
            DO i = 1, 3
              T(4-i) = tx1 - X(i+1)
              T(np+i) = txn - X(Ndata-i)
            END DO
          ELSE
            xl = (X(Ndata)-X(1))/3._SP
            DO i = 1, 3
              T(4-i) = T(5-i) - xl
              T(np+i) = T(np+i-1) + xl
            END DO
          END IF
        CASE (3)
          !     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE
          !     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3
          DO i = 1, 3
            T(4-i) = W(4-i,1)
            jw = MAX(1,i-1)
            iw = MOD(i+2,5) + 1
            T(np+i) = W(iw,jw)
            IF( T(4-i)>T(5-i) ) GOTO 100
            IF( T(np+i)<T(np+i-1) ) GOTO 100
          END DO
        CASE DEFAULT
          !     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA)
          DO i = 1, 3
            T(4-i) = X(1)
            T(np+i) = X(Ndata)
          END DO
      END SELECT
      !
      DO i = 1, 5
        DO j = 1, N
          W(i,j) = 0._SP
        END DO
      END DO
      !     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR
      !     RIGHT LIMITS
      it = Ibcl + 1
      CALL BSPVD(T,K,it,X(1),K,4,vnikx,work)
      iw = 0
      IF( ABS(vnikx(3,1))<tol ) iw = 1
      DO j = 1, 3
        W(j+1,4-j) = vnikx(4-j,it)
        W(j,4-j) = vnikx(4-j,1)
      END DO
      Bcoef(1) = Y(1)
      Bcoef(2) = Fbcl
      !     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1
      ileft = 4
      IF( ndm>=2 ) THEN
        DO i = 2, ndm
          ileft = ileft + 1
          CALL BSPVD(T,K,1,X(i),ileft,4,vnikx,work)
          DO j = 1, 3
            W(j+1,3+i-j) = vnikx(4-j,1)
          END DO
          Bcoef(i+1) = Y(i)
        END DO
      END IF
      !     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR
      !     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1))
      it = Ibcr + 1
      CALL BSPVD(T,K,it,X(Ndata),ileft,4,vnikx,work)
      jw = 0
      IF( ABS(vnikx(2,1))<tol ) jw = 1
      DO j = 1, 3
        W(j+1,3+Ndata-j) = vnikx(5-j,it)
        W(j+2,3+Ndata-j) = vnikx(5-j,1)
      END DO
      Bcoef(N-1) = Fbcr
      Bcoef(N) = Y(Ndata)
      !     SOLVE SYSTEM OF EQUATIONS
      ilb = 2 - jw
      iub = 2 - iw
      nwrow = 5
      iwp = iw + 1
      CALL BNFAC(W(iwp,1),nwrow,N,ilb,iub,iflag)
      IF( iflag==2 ) THEN
        !
        !
        ERROR STOP 'BINT4 : THE SYSTEM OF EQUATIONS IS SINGULAR'
        RETURN
      ELSE
        CALL BNSLV(W(iwp,1),nwrow,N,ilb,iub,Bcoef)
        RETURN
      END IF
    END IF
    50 ERROR STOP 'BINT4 : X VALUES ARE NOT DISTINCT OR NOT ORDERED'
    RETURN
  END IF
  100  ERROR STOP 'BINT4 : KNOT INPUT THROUGH W ARRAY IS NOT ORDERED PROPERLY'

END SUBROUTINE BINT4