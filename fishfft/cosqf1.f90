!*==COSQF1.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK COSQF1
SUBROUTINE COSQF1(N,X,W,Xh)
  IMPLICIT NONE
  !*--COSQF15
  !*** Start of declarations inserted by SPAG
  INTEGER i, k, kc, modn, N, np2, ns2
  REAL W, X, Xh, xim1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  COSQF1
  !***SUBSIDIARY
  !***PURPOSE  Compute the forward cosine transform with odd wave numbers.
  !***LIBRARY   SLATEC (FFTPACK)
  !***CATEGORY  J1A3
  !***TYPE      SINGLE PRECISION (COSQF1-S)
  !***KEYWORDS  FFTPACK, FOURIER TRANSFORM
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***DESCRIPTION
  !
  !  Subroutine COSQF1 computes the fast Fourier transform of quarter
  !  wave data. That is, COSQF1 computes the coefficients in a cosine
  !  series representation with only odd wave numbers.  The transform
  !  is defined below at Output Parameter X
  !
  !***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***ROUTINES CALLED  RFFTF
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  COSQF1
  DIMENSION X(*), W(*), Xh(*)
  !***FIRST EXECUTABLE STATEMENT  COSQF1
  ns2 = (N+1)/2
  np2 = N + 2
  DO k = 2, ns2
    kc = np2 - k
    Xh(k) = X(k) + X(kc)
    Xh(kc) = X(k) - X(kc)
  ENDDO
  modn = MOD(N,2)
  IF ( modn==0 ) Xh(ns2+1) = X(ns2+1) + X(ns2+1)
  DO k = 2, ns2
    kc = np2 - k
    X(k) = W(k-1)*Xh(kc) + W(kc-1)*Xh(k)
    X(kc) = W(k-1)*Xh(k) - W(kc-1)*Xh(kc)
  ENDDO
  IF ( modn==0 ) X(ns2+1) = W(ns2)*Xh(ns2+1)
  CALL RFFTF(N,X,Xh)
  DO i = 3, N, 2
    xim1 = X(i-1) - X(i)
    X(i) = X(i-1) + X(i)
    X(i-1) = xim1
  ENDDO
END SUBROUTINE COSQF1
