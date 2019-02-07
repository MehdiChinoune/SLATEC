!*==COSQF.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK COSQF
SUBROUTINE COSQF(N,X,Wsave)
  IMPLICIT NONE
  !*--COSQF5
  !*** Start of declarations inserted by SPAG
  INTEGER N
  REAL sqrt2 , tsqx , Wsave , X
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  COSQF
  !***PURPOSE  Compute the forward cosine transform with odd wave numbers.
  !***LIBRARY   SLATEC (FFTPACK)
  !***CATEGORY  J1A3
  !***TYPE      SINGLE PRECISION (COSQF-S)
  !***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***DESCRIPTION
  !
  !  Subroutine COSQF computes the fast Fourier transform of quarter
  !  wave data. That is, COSQF computes the coefficients in a cosine
  !  series representation with only odd wave numbers.  The transform
  !  is defined below at Output Parameter X
  !
  !  COSQF is the unnormalized inverse of COSQB since a call of COSQF
  !  followed by a call of COSQB will multiply the input sequence X
  !  by 4*N.
  !
  !  The array WSAVE which is used by subroutine COSQF must be
  !  initialized by calling subroutine COSQI(N,WSAVE).
  !
  !
  !  Input Parameters
  !
  !  N       the length of the array X to be transformed.  The method
  !          is most efficient when N is a product of small primes.
  !
  !  X       an array which contains the sequence to be transformed
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15
  !          in the program that calls COSQF.  The WSAVE array must be
  !          initialized by calling subroutine COSQI(N,WSAVE), and a
  !          different WSAVE array must be used for each different
  !          value of N.  This initialization does not have to be
  !          repeated so long as N remains unchanged.  Thus subsequent
  !          transforms can be obtained faster than the first.
  !
  !  Output Parameters
  !
  !  X       For I=1,...,N
  !
  !               X(I) = X(1) plus the sum from K=2 to K=N of
  !
  !                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
  !
  !               A call of COSQF followed by a call of
  !               COSQB will multiply the sequence X by 4*N.
  !               Therefore COSQB is the unnormalized inverse
  !               of COSQF.
  !
  !  WSAVE   contains initialization calculations which must not
  !          be destroyed between calls of COSQF or COSQB.
  !
  !***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***ROUTINES CALLED  COSQF1
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variable SQRT2 by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  COSQF
  DIMENSION X(*) , Wsave(*)
  !***FIRST EXECUTABLE STATEMENT  COSQF
  sqrt2 = SQRT(2.)
  IF ( N<2 ) THEN
  ELSEIF ( N==2 ) THEN
    tsqx = sqrt2*X(2)
    X(2) = X(1) - tsqx
    X(1) = X(1) + tsqx
  ELSE
    CALL COSQF1(N,X,Wsave,Wsave(N+1))
    GOTO 99999
  ENDIF
  RETURN
  99999 END SUBROUTINE COSQF
