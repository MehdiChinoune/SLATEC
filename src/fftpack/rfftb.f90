!** RFFTB
SUBROUTINE RFFTB(N,R,Wsave)
  IMPLICIT NONE
  !>
  !***
  !  Compute the backward fast Fourier transform of a real
  !            coefficient array.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A1
  !***
  ! **Type:**      SINGLE PRECISION (RFFTB-S, CFFTB-C)
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
  !   *   requested to use RFFTB1.                                       *
  !   *                                                                  *
  !   ********************************************************************
  !
  !   Subroutine RFFTB computes the real periodic sequence from its
  !   Fourier coefficients (Fourier synthesis).  The transform is defined
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
  !           in the program that calls RFFTB.  The WSAVE array must be
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
  !   R       For N even and for I = 1,...,N
  !
  !                R(I) = R(1)+(-1)**(I-1)*R(N)
  !
  !                     plus the sum from K=2 to K=N/2 of
  !
  !                      2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                     -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !           For N odd and for I = 1,...,N
  !
  !                R(I) = R(1) plus the sum from K=2 to K=(N+1)/2 of
  !
  !                     2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
  !
  !                    -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
  !
  !   Note:  This transform is unnormalized since a call of RFFTF
  !          followed by a call of RFFTB will multiply the input
  !          sequence by N.
  !
  !   WSAVE  contains results which must not be destroyed between
  !          calls of RFFTB or RFFTF.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTB1

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
  
  INTEGER N
  REAL R(*), Wsave(*)
  INTEGER :: ifac(15)
  !* FIRST EXECUTABLE STATEMENT  RFFTB
  IF ( N==1 ) RETURN
  ifac = INT( Wsave(2*N+1:2*N+15) )
  CALL RFFTB1(N,R,Wsave,Wsave(N+1),ifac)
  Wsave(2*N+1:2*N+15) = ifac
END SUBROUTINE RFFTB
