!** COSQB1
SUBROUTINE COSQB1(N,X,W,Xh)
  IMPLICIT NONE
  !>
  !***
  !  Compute the unnormalized inverse of COSQF1.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (COSQB1-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine COSQB1 computes the fast Fourier transform of quarter
  !  wave data. That is, COSQB1 computes a sequence from its
  !  representation in terms of a cosine series with odd wave numbers.
  !  The transform is defined below at output parameter X.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTB

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER i, k, kc, modn, N, np2, ns2
  REAL W(*), X(*), Xh(*), xim1
  !* FIRST EXECUTABLE STATEMENT  COSQB1
  ns2 = (N+1)/2
  np2 = N + 2
  DO i = 3, N, 2
    xim1 = X(i-1) + X(i)
    X(i) = X(i) - X(i-1)
    X(i-1) = xim1
  ENDDO
  X(1) = X(1) + X(1)
  modn = MOD(N,2)
  IF ( modn==0 ) X(N) = X(N) + X(N)
  CALL RFFTB(N,X,Xh)
  DO k = 2, ns2
    kc = np2 - k
    Xh(k) = W(k-1)*X(kc) + W(kc-1)*X(k)
    Xh(kc) = W(k-1)*X(k) - W(kc-1)*X(kc)
  ENDDO
  IF ( modn==0 ) X(ns2+1) = W(ns2)*(X(ns2+1)+X(ns2+1))
  DO k = 2, ns2
    kc = np2 - k
    X(k) = Xh(k) + Xh(kc)
    X(kc) = Xh(k) - Xh(kc)
  ENDDO
  X(1) = X(1) + X(1)
END SUBROUTINE COSQB1
