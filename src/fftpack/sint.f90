!** SINT
SUBROUTINE SINT(N,X,Wsave)
  IMPLICIT NONE
  !>
  !***
  !  Compute the sine transform of a real, odd sequence.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A3
  !***
  ! **Type:**      SINGLE PRECISION (SINT-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine SINT computes the discrete Fourier sine transform
  !  of an odd sequence X(I).  The transform is defined below at
  !  output parameter X.
  !
  !  SINT is the unnormalized inverse of itself since a call of SINT
  !  followed by another call of SINT will multiply the input sequence
  !  X by 2*(N+1).
  !
  !  The array WSAVE which is used by subroutine SINT must be
  !  initialized by calling subroutine SINTI(N,WSAVE).
  !
  !  Input Parameters
  !
  !  N       the length of the sequence to be transformed.  The method
  !          is most efficient when N+1 is the product of small primes.
  !
  !  X       an array which contains the sequence to be transformed
  !
  !
  !  WSAVE   a work array with dimension at least INT(3.5*N+16)
  !          in the program that calls SINT.  The WSAVE array must be
  !          initialized by calling subroutine SINTI(N,WSAVE), and a
  !          different WSAVE array must be used for each different
  !          value of N.  This initialization does not have to be
  !          repeated so long as N remains unchanged.  Thus subsequent
  !          transforms can be obtained faster than the first.
  !
  !  Output Parameters
  !
  !  X       For I=1,...,N
  !
  !               X(I)= the sum from K=1 to K=N
  !
  !                    2*X(K)*SIN(K*I*PI/(N+1))
  !
  !               A call of SINT followed by another call of
  !               SINT will multiply the sequence X by 2*(N+1).
  !               Hence SINT is the unnormalized inverse
  !               of itself.
  !
  !  WSAVE   contains initialization calculations which must not be
  !          destroyed between calls of SINT.
  !
  !***
  ! **References:**  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
  !                 Computations (G. Rodrigue, ed.), Academic Press,
  !                 1982, pp. 51-83.
  !***
  ! **Routines called:**  RFFTF

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variable SQRT3 by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER i, k, kc, kw, modn, N, nf, np1, ns2
  REAL sqrt3, t1, t2, Wsave, X, xh
  DIMENSION X(*), Wsave(*)
  !* FIRST EXECUTABLE STATEMENT  SINT
  sqrt3 = SQRT(3.)
  IF ( N<2 ) THEN
    X(1) = X(1) + X(1)
    RETURN
  ELSEIF ( N==2 ) THEN
    xh = sqrt3*(X(1)+X(2))
    X(2) = sqrt3*(X(1)-X(2))
    X(1) = xh
    RETURN
  ELSE
    np1 = N + 1
    ns2 = N/2
    Wsave(1) = 0.
    kw = np1
    DO k = 1, ns2
      kw = kw + 1
      kc = np1 - k
      t1 = X(k) - X(kc)
      t2 = Wsave(kw)*(X(k)+X(kc))
      Wsave(k+1) = t1 + t2
      Wsave(kc+1) = t2 - t1
    ENDDO
    modn = MOD(N,2)
    IF ( modn/=0 ) Wsave(ns2+2) = 4.*X(ns2+1)
    nf = np1 + ns2 + 1
    CALL RFFTF(np1,Wsave,Wsave(nf))
    X(1) = .5*Wsave(1)
    DO i = 3, N, 2
      X(i-1) = -Wsave(i)
      X(i) = X(i-2) + Wsave(i-1)
    ENDDO
    IF ( modn/=0 ) RETURN
    X(N) = -Wsave(N+1)
  ENDIF
END SUBROUTINE SINT
