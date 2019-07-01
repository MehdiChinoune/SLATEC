!** EZFFTI
PURE SUBROUTINE EZFFTI(N,Wsave)
  !> Initialize a work array for EZFFTF and EZFFTB.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A1
  !***
  ! **Type:**      SINGLE PRECISION (EZFFTI-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine EZFFTI initializes the work array WSAVE which is used in
  !  both EZFFTF and EZFFTB.  The prime factorization of N together with
  !  a tabulation of the trigonometric functions are computed and
  !  stored in WSAVE.
  !
  !  Input Parameter
  !
  !  N       the length of the sequence to be transformed.
  !
  !  Output Parameter
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15.
  !          The same work array can be used for both EZFFTF and EZFFTB
  !          as long as N remains unchanged.  Different WSAVE arrays
  !          are required for different values of N.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  EZFFT1

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER, INTENT(IN) :: N
  REAL(SP), INTENT(OUT) :: Wsave(3*N+15)
  INTEGER :: ifac(15)
  !* FIRST EXECUTABLE STATEMENT  EZFFTI
  IF( N==1 ) RETURN
  ifac = INT( Wsave(3*N+1:3*N+15) )
  CALL EZFFT1(N,Wsave(2*N+1),ifac)
  Wsave(3*N+1:3*N+15) = ifac

END SUBROUTINE EZFFTI