!DECK CFFTI
SUBROUTINE CFFTI(N,Wsave)
  IMPLICIT NONE
  INTEGER iw1, iw2, N
  REAL Wsave
  !***BEGIN PROLOGUE  CFFTI
  !***SUBSIDIARY
  !***PURPOSE  Initialize a work array for CFFTF and CFFTB.
  !***LIBRARY   SLATEC (FFTPACK)
  !***CATEGORY  J1A2
  !***TYPE      COMPLEX (RFFTI-S, CFFTI-C)
  !***KEYWORDS  FFTPACK, FOURIER TRANSFORM
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***DESCRIPTION
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
  !***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***ROUTINES CALLED  CFFTI1
  !***REVISION HISTORY  (YYMMDD)
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
  !***END PROLOGUE  CFFTI
  DIMENSION Wsave(*)
  !***FIRST EXECUTABLE STATEMENT  CFFTI
  IF ( N==1 ) RETURN
  iw1 = N + N + 1
  iw2 = iw1 + N + N
  CALL CFFTI1(N,Wsave(iw1),Wsave(iw2))
END SUBROUTINE CFFTI
