!** EZFFTF
PURE SUBROUTINE EZFFTF(N,R,Azero,A,B,Wsave)
  !> Compute a simplified REAL(SP), periodic, fast Fourier forward transform.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Category:**  J1A1
  !***
  ! **Type:**      SINGLE PRECISION (EZFFTF-S)
  !***
  ! **Keywords:**  FFTPACK, FOURIER TRANSFORM
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Description:**
  !
  !  Subroutine EZFFTF computes the Fourier coefficients of a real
  !  periodic sequence (Fourier analysis).  The transform is defined
  !  below at Output Parameters AZERO, A and B.  EZFFTF is a simplified
  !  but slower version of RFFTF.
  !
  !  Input Parameters
  !
  !  N       the length of the array R to be transformed.  The method
  !          is most efficient when N is the product of small primes.
  !
  !  R       a real array of length N which contains the sequence
  !          to be transformed.  R is not destroyed.
  !
  !
  !  WSAVE   a work array which must be dimensioned at least 3*N+15
  !          in the program that calls EZFFTF.  The WSAVE array must be
  !          initialized by calling subroutine EZFFTI(N,WSAVE), and a
  !          different WSAVE array must be used for each different
  !          value of N.  This initialization does not have to be
  !          repeated so long as N remains unchanged.  Thus subsequent
  !          transforms can be obtained faster than the first.
  !          The same WSAVE array can be used by EZFFTF and EZFFTB.
  !
  !  Output Parameters
  !
  !  AZERO   the sum from I=1 to I=N of R(I)/N
  !
  !  A,B     for N even B(N/2)=0. and A(N/2) is the sum from I=1 to
  !          I=N of (-1)**(I-1)*R(I)/N
  !
  !          for N even define KMAX=N/2-1
  !          for N odd  define KMAX=(N-1)/2
  !
  !          then for  K=1,...,KMAX
  !
  !               A(K) equals the sum from I=1 to I=N of
  !
  !                    2./N*R(I)*COS(K*(I-1)*2*PI/N)
  !
  !               B(K) equals the sum from I=1 to I=N of
  !
  !                    2./N*R(I)*SIN(K*(I-1)*2*PI/N)
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
  !           (b) changing references to intrinsic function FLOAT
  !               to REAL.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER, INTENT(IN) :: N
  REAL(SP), INTENT(IN) :: R(N)
  REAL(SP), INTENT(INOUT) :: Wsave(3*N+15)
  REAL(SP), INTENT(OUT) :: A(N/2), Azero, B(N/2)
  REAL(SP) :: cf, cfm
  INTEGER :: i, ns2, ns2m
  !* FIRST EXECUTABLE STATEMENT  EZFFTF
  IF( N<2 ) THEN
    Azero = R(1)
    RETURN
  ELSEIF( N==2 ) THEN
    Azero = 0.5_SP*(R(1)+R(2))
    A(1) = 0.5_SP*(R(1)-R(2))
    RETURN
  ELSE
    DO i = 1, N
      Wsave(i) = R(i)
    END DO
    CALL RFFTF(N,Wsave,Wsave(N+1))
    cf = 2._SP/N
    cfm = -cf
    Azero = 0.5_SP*cf*Wsave(1)
    ns2 = (N+1)/2
    ns2m = ns2 - 1
    DO i = 1, ns2m
      A(i) = cf*Wsave(2*i)
      B(i) = cfm*Wsave(2*i+1)
    END DO
    IF( MOD(N,2)==0 ) A(ns2) = 0.5_SP*cf*Wsave(N)
  END IF

END SUBROUTINE EZFFTF