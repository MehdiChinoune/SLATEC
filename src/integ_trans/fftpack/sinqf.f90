!** SINQF
SUBROUTINE SINQF(N,X,Wsave)
  !>
  !  Compute the forward sine transform with odd wave numbers.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (SINQF-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine SINQF computes the fast Fourier transform of quarter
  !  wave data.  That is, SINQF computes the coefficients in a sine
  !  series representation with only odd wave numbers.  The transform
  !  is defined below at output parameter X.
  !
  !  SINQB is the unnormalized inverse of SINQF since a call of SINQF
  !  followed by a call of SINQB will multiply the input sequence X
  !  by 4*N.
  !
  !  The array WSAVE which is used by subroutine SINQF must be
  !  initialized by calling subroutine SINQI(N,WSAVE).
  !
  !  Input Parameters
  !
  !  N       the length of the array X to be transformed.  The method
  !          is most efficient when N is a product of small primes.
  !
  !  X       an array which contains the sequence to be transformed
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15
  !          in the program that calls SINQF.  The WSAVE array must be
  !          initialized by calling subroutine SINQI(N,WSAVE), and a
  !          different WSAVE array must be used for each different
  !          value of N.  This initialization does not have to be
  !          repeated so long as N remains unchanged.  Thus subsequent
  !          transforms can be obtained faster than the first.
  !
  !  Output Parameters
  !
  !  X       For I=1,...,N
  !
  !               X(I) = (-1)**(I-1)*X(N)
  !
  !                  + the sum from K=1 to K=N-1 of
  !
  !                  2*X(K)*SIN((2*I-1)*K*PI/(2*N))
  !
  !               A call of SINQF followed by a call of
  !               SINQB will multiply the sequence X by 4*N.
  !               Therefore SINQB is the unnormalized inverse
  !               of SINQF.
  !
  !  WSAVE   contains initialization calculations which must not
  !          be destroyed between calls of SINQF or SINQB.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  COSQF

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*)
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: N
  REAL(SP) :: Wsave(N), X(N)
  INTEGER :: k, kc, ns2
  REAL(SP) :: xhold
  !* FIRST EXECUTABLE STATEMENT  SINQF
  IF ( N==1 ) RETURN
  ns2 = N/2
  DO k = 1, ns2
    kc = N - k
    xhold = X(k)
    X(k) = X(kc+1)
    X(kc+1) = xhold
  END DO
  CALL COSQF(N,X,Wsave)
  DO k = 2, N, 2
    X(k) = -X(k)
  END DO
END SUBROUTINE SINQF
