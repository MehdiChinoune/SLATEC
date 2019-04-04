!** COSTI
SUBROUTINE COSTI(N,Wsave)
  IMPLICIT NONE
  !>
  !***
  !  Initialize a work array for COST.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (COSTI-S)
  !***
  ! **Keywords:**  COSINE FOURIER TRANSFORM, FFTPACK
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine COSTI initializes the array WSAVE which is used in
  !  subroutine COST.  The prime factorization of N together with
  !  a tabulation of the trigonometric functions are computed and
  !  stored in WSAVE.
  !
  !  Input Parameter
  !
  !  N       the length of the sequence to be transformed.  The method
  !          is most efficient when N-1 is a product of small primes.
  !
  !  Output Parameter
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15.
  !          Different WSAVE arrays are required for different values
  !          of N.  The contents of WSAVE must not be changed between
  !          calls of COST.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTI

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing references to intrinsic function FLOAT
  !               to REAL, and
  !           (c) changing definition of variable PI by using
  !               FORTRAN intrinsic function ATAN instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  REAL dt, fk, pi, Wsave(*)
  INTEGER k, kc, N, nm1, np1, ns2
  !* FIRST EXECUTABLE STATEMENT  COSTI
  IF ( N<=3 ) RETURN
  pi = 4.*ATAN(1.)
  nm1 = N - 1
  np1 = N + 1
  ns2 = N/2
  dt = pi/nm1
  fk = 0.
  DO k = 2, ns2
    kc = np1 - k
    fk = fk + 1.
    Wsave(k) = 2.*SIN(fk*dt)
    Wsave(kc) = 2.*COS(fk*dt)
  END DO
  CALL RFFTI(nm1,Wsave(N+1))
END SUBROUTINE COSTI
