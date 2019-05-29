!** POLYVL
SUBROUTINE POLYVL(Nder,Xx,Yfit,Yp,N,X,C,Work,Ierr)
  !>
  !  Calculate the value of a polynomial and its first NDER
  !            derivatives where the polynomial was produced by a previous
  !            call to POLINT.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E3
  !***
  ! **Type:**      SINGLE PRECISION (POLYVL-S, DPOLVL-D)
  !***
  ! **Keywords:**  POLYNOMIAL EVALUATION
  !***
  ! **Author:**  Huddleston, R. E., (SNLL)
  !***
  ! **Description:**
  !
  !     Written by Robert E. Huddleston, Sandia Laboratories, Livermore
  !
  !     Abstract -
  !        Subroutine POLYVL calculates the value of the polynomial and
  !     its first NDER derivatives where the polynomial was produced by
  !     a previous call to POLINT.
  !        The variable N and the arrays X and C must not be altered
  !     between the call to POLINT and the call to POLYVL.
  !
  !     ******  Dimensioning Information *******
  !
  !     YP   must be dimensioned by at least NDER
  !     X    must be dimensioned by at least N (see the abstract )
  !     C    must be dimensioned by at least N (see the abstract )
  !     WORK must be dimensioned by at least 2*N if NDER is .GT. 0.
  !
  !     *** Note ***
  !       If NDER=0, neither YP nor WORK need to be dimensioned variables.
  !       If NDER=1, YP does not need to be a dimensioned variable.
  !
  !
  !     *****  Input parameters
  !
  !     NDER - the number of derivatives to be evaluated
  !
  !     XX   - the argument at which the polynomial and its derivatives
  !            are to be evaluated.
  !
  !     N    - *****
  !            *       N, X, and C must not be altered between the call
  !     X    - *       to POLINT and the call to POLYVL.
  !     C    - *****
  !
  !
  !     *****  Output Parameters
  !
  !     YFIT - the value of the polynomial at XX
  !
  !     YP   - the derivatives of the polynomial at XX.  The derivative of
  !            order J at XX is stored in  YP(J), J = 1,...,NDER.
  !
  !     IERR - Output error flag with the following possible values.
  !          = 1  indicates normal execution
  !
  !     ***** Storage Parameters
  !
  !     WORK  = this is an array to provide internal working storage for
  !             POLYVL.  It must be dimensioned by at least 2*N if NDER is
  !             .GT. 0.  If NDER=0, WORK does not need to be a dimensioned
  !             variable.
  !
  !***
  ! **References:**  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: Ierr, N, Nder
  REAL :: C(N), Work(2*N), X(N), Yp(Nder), Xx, Yfit
  INTEGER i, im1, izero, k, km1, km1pi, km2pn, km2pni, m, mm, ndr, nmkp1, npkm1
  REAL :: fac, pione, pitwo, pone, ptwo, xk
  !* FIRST EXECUTABLE STATEMENT  POLYVL
  Ierr = 1
  IF ( Nder<=0 ) THEN
    !
    !     *****   CODING FOR THE CASE NDER = 0
    !
    pione = 1.0
    pone = C(1)
    Yfit = pone
    IF ( N==1 ) RETURN
    DO k = 2, N
      pitwo = (Xx-X(k-1))*pione
      pione = pitwo
      ptwo = pone + pitwo*C(k)
      pone = ptwo
    END DO
    Yfit = ptwo
    RETURN
    !
    !     *****   END OF NDER = 0 CASE
    !
  ELSEIF ( N>1 ) THEN
    !
    !     *****  END OF THE CASE  N = 1 AND  NDER .GT. 0
    !
    IF ( Nder<N ) THEN
      izero = 0
      ndr = Nder
    ELSE
      !
      !     *****  SET FLAGS FOR NUMBER OF DERIVATIVES AND FOR DERIVATIVES
      !            IN EXCESS OF THE DEGREE (N-1) OF THE POLYNOMIAL.
      !
      izero = 1
      ndr = N - 1
    END IF
    m = ndr + 1
    mm = m
    !
    !     *****  START OF THE CASE NDER .GT. 0  AND N .GT. 1
    !     *****  THE POLYNOMIAL AND ITS DERIVATIVES WILL BE EVALUATED AT XX
    !
    DO k = 1, ndr
      Yp(k) = C(k+1)
    END DO
    !
    !     *****  THE FOLLOWING SECTION OF CODE IS EASIER TO READ IF ONE
    !            BREAKS WORK INTO TWO ARRAYS W AND V. THE CODE WOULD THEN
    !            READ
    !                W(1) = 1.
    !                PONE = C(1)
    !               *DO   K = 2,N
    !               *   V(K-1) =  XX - X(K-1)
    !               *   W(K)   =  V(K-1)*W(K-1)
    !               *   PTWO   =  PONE + W(K)*C(K)
    !               *   PONE   =  PWO
    !
    !               YFIT = PTWO
    !
    Work(1) = 1.0
    pone = C(1)
    DO k = 2, N
      km1 = k - 1
      npkm1 = N + k - 1
      Work(npkm1) = Xx - X(km1)
      Work(k) = Work(npkm1)*Work(km1)
      ptwo = pone + Work(k)*C(k)
      pone = ptwo
    END DO
    Yfit = ptwo
    !
    !     ** AT THIS POINT THE POLYNOMIAL HAS BEEN EVALUATED AND INFORMATION
    !        FOR THE DERIVATIVE EVALUATIONS HAVE BEEN STORED IN THE ARRAY
    !        WORK
    IF ( N/=2 ) THEN
      IF ( m==N ) mm = ndr
      !
      !     ***** EVALUATE THE DERIVATIVES AT XX
      !
      !                  ******  DO K=2,MM   (FOR MOST CASES, MM = NDER + 1)
      !                  *  ******  DO I=2,N-K+1
      !                  *  *       W(I) = V(K-2+I)*W(I-1) + W(I)
      !                  *  *       YP(K-1) = YP(K-1) + W(I)*C(K-1+I)
      !                  ******  CONTINUE
      !
      DO k = 2, mm
        nmkp1 = N - k + 1
        km1 = k - 1
        km2pn = k - 2 + N
        DO i = 2, nmkp1
          km2pni = km2pn + i
          im1 = i - 1
          km1pi = km1 + i
          Work(i) = Work(km2pni)*Work(im1) + Work(i)
          Yp(km1) = Yp(km1) + Work(i)*C(km1pi)
        END DO
      END DO
      IF ( ndr/=1 ) THEN
        fac = 1.0
        DO k = 2, ndr
          xk = k
          fac = xk*fac
          Yp(k) = fac*Yp(k)
        END DO
      END IF
    END IF
    !
    !     ***** END OF DERIVATIVE EVALUATIONS
    !
    IF ( izero==0 ) RETURN
    !
    !     *****  SET EXCESS DERIVATIVES TO ZERO.
    !
    DO k = N, Nder
      Yp(k) = 0.0
    END DO
  ELSE
    Yfit = C(1)
    !
    !     *****  CODING FOR THE CASE  N=1 AND NDER .GT. 0
    !
    DO k = 1, Nder
      Yp(k) = 0.0
    END DO
    RETURN
  END IF
END SUBROUTINE POLYVL
