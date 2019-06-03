!** RFFTF
SUBROUTINE RFFTF(N,R,Wsave)
  !>
  !  Compute the forward transform of a REAL(SP), periodic sequence.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A1
  !***
  ! **Type:**      SINGLE PRECISION (RFFTF-S, CFFTF-C)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !   ********************************************************************
  !   *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
  !   ********************************************************************
  !   *                                                                  *
  !   *   This routine uses non-standard Fortran 77 constructs and will  *
  !   *   be removed from the library at a future date.  You are         *
  !   *   requested to use RFFTF1.                                       *
  !   *                                                                  *
  !   ********************************************************************
  !
  !   Subroutine RFFTF computes the Fourier coefficients of a real
  !   periodic sequence (Fourier analysis).  The transform is defined
  !   below at output parameter R.
  !
  !   Input Arguments
  !
  !   N       the length of the array R to be transformed.  The method
  !           is most efficient when N is a product of small primes.
  !           N may change so long as different work arrays are provided.
  !
  !   R       a real array of length N which contains the sequence
  !           to be transformed.
  !
  !   WSAVE   a work array which must be dimensioned at least 2*N+15
  !           in the program that calls RFFTF.  The WSAVE array must be
  !           initialized by calling subroutine RFFTI, and a different
  !           WSAVE array must be used for each different value of N.
  !           This initialization does not have to be repeated so long as
  !           remains unchanged.  Thus subsequent transforms can be
  !           obtained faster than the first.  Moreover, the same WSAVE
  !           array can be used by RFFTF and RFFTB as long as N remains
  !           unchanged.
  !
  !   Output Argument
  !
  !   R       R(1) = the sum from I=1 to I=N of R(I)
  !
  !           If N is even set L = N/2; if N is odd set L = (N+1)/2
  !
  !             then for K = 2,...,L
  !
  !                R(2*K-2) = the sum from I = 1 to I = N of
  !
  !                     R(I)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                R(2*K-1) = the sum from I = 1 to I = N of
  !
  !                    -R(I)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !           If N is even
  !
  !                R(N) = the sum from I = 1 to I = N of
  !
  !                     (-1)**(I-1)*R(I)
  !
  !   Note:  This transform is unnormalized since a call of RFFTF
  !          followed by a call of RFFTB will multiply the input
  !          sequence by N.
  !
  !   WSAVE  contains results which must not be destroyed between
  !          calls of RFFTF or RFFTB.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTF1

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
  REAL(SP) :: R(N), Wsave(2*N+15)
  INTEGER :: ifac(15)
  !* FIRST EXECUTABLE STATEMENT  RFFTF
  IF ( N==1 ) RETURN
  ifac = INT( Wsave(2*N+1:2*N+15) )
  CALL RFFTF1(N,R,Wsave,Wsave(N+1),ifac)
END SUBROUTINE RFFTF
