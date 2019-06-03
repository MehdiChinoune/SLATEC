!** COSQB
SUBROUTINE COSQB(N,X,Wsave)
  !>
  !  Compute the unnormalized inverse cosine transform.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (COSQB-S)
  !***
  ! **Keywords:**  FFTPACK, INVERSE COSINE FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine COSQB computes the fast Fourier transform of quarter
  !  wave data. That is, COSQB computes a sequence from its
  !  representation in terms of a cosine series with odd wave numbers.
  !  The transform is defined below at output parameter X.
  !
  !  COSQB is the unnormalized inverse of COSQF since a call of COSQB
  !  followed by a call of COSQF will multiply the input sequence X
  !  by 4*N.
  !
  !  The array WSAVE which is used by subroutine COSQB must be
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
  !          in the program that calls COSQB.  The WSAVE array must be
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
  !               X(I)= the sum from K=1 to K=N of
  !
  !                  2*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
  !
  !               A call of COSQB followed by a call of
  !               COSQF will multiply the sequence X by 4*N.
  !               Therefore COSQF is the unnormalized inverse
  !               of COSQB.
  !
  !  WSAVE   contains initialization calculations which must not
  !          be destroyed between calls of COSQB or COSQF.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  COSQB1

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variable TSQRT2 by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   861211  REVISION DATE from Version 3.2
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N
  REAL(SP) :: Wsave(3*N+15), X(N)
  REAL(SP) :: tsqrt2, x1
  !* FIRST EXECUTABLE STATEMENT  COSQB
  tsqrt2 = 2.*SQRT(2.)
  IF ( N<2 ) THEN
    X(1) = 4.*X(1)
    RETURN
  ELSEIF ( N==2 ) THEN
    x1 = 4.*(X(1)+X(2))
    X(2) = tsqrt2*(X(1)-X(2))
    X(1) = x1
    RETURN
  ELSE
    CALL COSQB1(N,X,Wsave,Wsave(N+1))
  END IF
END SUBROUTINE COSQB
