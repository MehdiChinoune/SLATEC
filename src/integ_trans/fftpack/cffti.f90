!** CFFTI
SUBROUTINE CFFTI(N,Wsave)
  !>
  !  Initialize a work array for CFFTF and CFFTB.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A2
  !***
  ! **Type:**      COMPLEX (RFFTI-S, CFFTI-C)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  ********************************************************************
  !  *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
  !  ********************************************************************
  !  *                                                                  *
  !  *   This routine uses non-standard Fortran 77 constructs and will  *
  !  *   be removed from the library at a future date.  You are         *
  !  *   requested to use CFFTI1.                                       *
  !  *                                                                  *
  !  ********************************************************************
  !
  !  Subroutine CFFTI initializes the array WSAVE which is used in
  !  both CFFTF and CFFTB.  The prime factorization of N together with
  !  a tabulation of the trigonometric functions are computed and
  !  stored in WSAVE.
  !
  !  Input Parameter
  !
  !  N       the length of the sequence to be transformed
  !
  !  Output Parameter
  !
  !  WSAVE   a work array which must be dimensioned at least 4*N+15.
  !          The same work array can be used for both CFFTF and CFFTB
  !          as long as N remains unchanged.  Different WSAVE arrays
  !          are required for different values of N.  The contents of
  !          WSAVE must not be changed between calls of CFFTF or CFFTB.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  CFFTI1

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900131  Routine changed from user-callable to subsidiary
  !           because of non-standard Fortran 77 arguments in the
  !           call to CFFTB1.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: N
  REAL(SP) :: Wsave(4*N+15)
  INTEGER :: ifac(15)
  !* FIRST EXECUTABLE STATEMENT  CFFTI
  IF ( N==1 ) RETURN
  ifac = INT( Wsave(4*N+1:4*N+15) )
  CALL CFFTI1(N,Wsave(2*N+1),ifac)
  Wsave(4*N+1:4*N+15) = ifac
END SUBROUTINE CFFTI
