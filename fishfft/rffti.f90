!*==RFFTI.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK RFFTI
      SUBROUTINE RFFTI(N,Wsave)
      IMPLICIT NONE
!*--RFFTI5
!*** Start of declarations inserted by SPAG
      INTEGER N
      REAL Wsave
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  RFFTI
!***SUBSIDIARY
!***PURPOSE  Initialize a work array for RFFTF and RFFTB.
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (RFFTI-S, CFFTI-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!   ********************************************************************
!   *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
!   ********************************************************************
!   *                                                                  *
!   *   This routine uses non-standard Fortran 77 constructs and will  *
!   *   be removed from the library at a future date.  You are         *
!   *   requested to use RFFTI1.                                       *
!   *                                                                  *
!   ********************************************************************
!
!   Subroutine RFFTI initializes the array WSAVE which is used in
!   both RFFTF and RFFTB.  The prime factorization of N together with
!   a tabulation of the trigonometric functions are computed and
!   stored in WSAVE.
!
!   Input Argument
!
!   N       the length of the sequence to be transformed.
!
!   Output Argument
!
!   WSAVE   a work array which must be dimensioned at least 2*N+15.
!           The same work array can be used for both RFFTF and RFFTB
!           as long as N remains unchanged.  Different WSAVE arrays
!           are required for different values of N.  The contents of
!           WSAVE must not be changed between calls of RFFTF or RFFTB.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTI1
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
!***END PROLOGUE  RFFTI
      DIMENSION Wsave(*)
!***FIRST EXECUTABLE STATEMENT  RFFTI
      IF ( N==1 ) RETURN
      CALL RFFTI1(N,Wsave(N+1),Wsave(2*N+1))
      END SUBROUTINE RFFTI
