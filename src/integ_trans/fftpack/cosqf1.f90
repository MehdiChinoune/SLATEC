!** COSQF1
SUBROUTINE COSQF1(N,X,W,Xh)
  !>
  !  Compute the forward cosine transform with odd wave numbers.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (COSQF1-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine COSQF1 computes the fast Fourier transform of quarter
  !  wave data. That is, COSQF1 computes the coefficients in a cosine
  !  series representation with only odd wave numbers.  The transform
  !  is defined below at Output Parameter X
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTF

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N
  REAL(SP) :: W(N), X(N), Xh(2*N+15)
  INTEGER i, k, kc, modn, np2, ns2
  REAL(SP) :: xim1
  !* FIRST EXECUTABLE STATEMENT  COSQF1
  ns2 = (N+1)/2
  np2 = N + 2
  DO k = 2, ns2
    kc = np2 - k
    Xh(k) = X(k) + X(kc)
    Xh(kc) = X(k) - X(kc)
  END DO
  modn = MOD(N,2)
  IF ( modn==0 ) Xh(ns2+1) = X(ns2+1) + X(ns2+1)
  DO k = 2, ns2
    kc = np2 - k
    X(k) = W(k-1)*Xh(kc) + W(kc-1)*Xh(k)
    X(kc) = W(k-1)*Xh(k) - W(kc-1)*Xh(kc)
  END DO
  IF ( modn==0 ) X(ns2+1) = W(ns2)*Xh(ns2+1)
  CALL RFFTF(N,X,Xh)
  DO i = 3, N, 2
    xim1 = X(i-1) - X(i)
    X(i) = X(i-1) + X(i)
    X(i-1) = xim1
  END DO
END SUBROUTINE COSQF1
