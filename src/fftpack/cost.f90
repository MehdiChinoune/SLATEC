!** COST
SUBROUTINE COST(N,X,Wsave)
  IMPLICIT NONE
  !>
  !***
  !  Compute the cosine transform of a real, even sequence.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (COST-S)
  !***
  ! **Keywords:**  COSINE FOURIER TRANSFORM, FFTPACK
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine COST computes the discrete Fourier cosine transform
  !  of an even sequence X(I).  The transform is defined below at output
  !  parameter X.
  !
  !  COST is the unnormalized inverse of itself since a call of COST
  !  followed by another call of COST will multiply the input sequence
  !  X by 2*(N-1).  The transform is defined below at output parameter X.
  !
  !  The array WSAVE which is used by subroutine COST must be
  !  initialized by calling subroutine COSTI(N,WSAVE).
  !
  !  Input Parameters
  !
  !  N       the length of the sequence X.  N must be greater than 1.
  !          The method is most efficient when N-1 is a product of
  !          small primes.
  !
  !  X       an array which contains the sequence to be transformed
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15
  !          in the program that calls COST.  The WSAVE array must be
  !          initialized by calling subroutine COSTI(N,WSAVE), and a
  !          different WSAVE array must be used for each different
  !          value of N.  This initialization does not have to be
  !          repeated so long as N remains unchanged.  Thus subsequent
  !          transforms can be obtained faster than the first.
  !
  !  Output Parameters
  !
  !  X       For I=1,...,N
  !
  !             X(I) = X(1)+(-1)**(I-1)*X(N)
  !
  !               + the sum from K=2 to K=N-1
  !
  !                 2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
  !
  !               A call of COST followed by another call of
  !               COST will multiply the sequence X by 2*(N-1).
  !               Hence COST is the unnormalized inverse
  !               of itself.
  !
  !  WSAVE   contains initialization calculations which must not be
  !          destroyed between calls of COST.
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
  !           changing dummy array size declarations (1) to (*)
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  REAL c1, t1, t2, tx2, Wsave(*), X(*), x1h, x1p3, xi, xim2
  INTEGER i, k, kc, modn, N, nm1, np1, ns2
  !* FIRST EXECUTABLE STATEMENT  COST
  nm1 = N - 1
  np1 = N + 1
  ns2 = N/2
  IF ( N<2 ) RETURN
  IF ( N==2 ) THEN
    x1h = X(1) + X(2)
    X(2) = X(1) - X(2)
    X(1) = x1h
    RETURN
  ELSEIF ( N>3 ) THEN
    c1 = X(1) - X(N)
    X(1) = X(1) + X(N)
    DO k = 2, ns2
      kc = np1 - k
      t1 = X(k) + X(kc)
      t2 = X(k) - X(kc)
      c1 = c1 + Wsave(kc)*t2
      t2 = Wsave(k)*t2
      X(k) = t1 - t2
      X(kc) = t1 + t2
    ENDDO
    modn = MOD(N,2)
    IF ( modn/=0 ) X(ns2+1) = X(ns2+1) + X(ns2+1)
    CALL RFFTF(nm1,X,Wsave(N+1))
    xim2 = X(2)
    X(2) = c1
    DO i = 4, N, 2
      xi = X(i)
      X(i) = X(i-2) - X(i-1)
      X(i-1) = xim2
      xim2 = xi
    ENDDO
    IF ( modn/=0 ) X(N) = xim2
    RETURN
  ENDIF
  x1p3 = X(1) + X(3)
  tx2 = X(2) + X(2)
  X(2) = X(1) - X(3)
  X(1) = x1p3 + tx2
  X(3) = x1p3 - tx2
  RETURN
END SUBROUTINE COST
