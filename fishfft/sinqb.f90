!*==SINQB.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK SINQB
      SUBROUTINE SINQB(N,X,Wsave)
      IMPLICIT NONE
!*--SINQB5
!*** Start of declarations inserted by SPAG
      INTEGER k , kc , N , ns2
      REAL Wsave , X , xhold
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SINQB
!***PURPOSE  Compute the unnormalized inverse of SINQF.
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (SINQB-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine SINQB computes the fast Fourier transform of quarter
!  wave data.  That is, SINQB computes a sequence from its
!  representation in terms of a sine series with odd wave numbers.
!  the transform is defined below at output parameter X.
!
!  SINQF is the unnormalized inverse of SINQB since a call of SINQB
!  followed by a call of SINQF will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine SINQB must be
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
!          in the program that calls SINQB.  The WSAVE array must be
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
!               X(I)= the sum from K=1 to K=N of
!
!                 4*X(K)*SIN((2*K-1)*I*PI/(2*N))
!
!               a call of SINQB followed by a call of
!               SINQF will multiply the sequence X by 4*N.
!               Therefore SINQF is the unnormalized inverse
!               of SINQB.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of SINQB or SINQF.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  COSQB
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   861211  REVISION DATE from Version 3.2
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINQB
      DIMENSION X(*) , Wsave(*)
!***FIRST EXECUTABLE STATEMENT  SINQB
      IF ( N>1 ) THEN
        ns2 = N/2
        DO k = 2 , N , 2
          X(k) = -X(k)
        ENDDO
        CALL COSQB(N,X,Wsave)
        DO k = 1 , ns2
          kc = N - k
          xhold = X(k)
          X(k) = X(kc+1)
          X(kc+1) = xhold
        ENDDO
        GOTO 99999
      ENDIF
      X(1) = 4.*X(1)
      RETURN
99999 END SUBROUTINE SINQB
