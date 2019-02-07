!*==POLCOF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK POLCOF
SUBROUTINE POLCOF(Xx,N,X,C,D,Work)
  IMPLICIT NONE
  !*--POLCOF5
  !*** Start of declarations inserted by SPAG
  REAL C , D , pone , ptwo , Work , X , Xx
  INTEGER i , im1 , k , km1 , km1pi , km2n , km2npi , N , nm1 , nmkp1 , &
    npkm1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  POLCOF
  !***PURPOSE  Compute the coefficients of the polynomial fit (including
  !            Hermite polynomial fits) produced by a previous call to
  !            POLINT.
  !***LIBRARY   SLATEC
  !***CATEGORY  E1B
  !***TYPE      SINGLE PRECISION (POLCOF-S, DPOLCF-D)
  !***KEYWORDS  COEFFICIENTS, POLYNOMIAL
  !***AUTHOR  Huddleston, R. E., (SNLL)
  !***DESCRIPTION
  !
  !     Written by Robert E. Huddleston, Sandia Laboratories, Livermore
  !
  !     Abstract
  !        Subroutine POLCOF computes the coefficients of the polynomial
  !     fit (including Hermite polynomial fits ) produced by a previous
  !     call to POLINT. The coefficients of the polynomial, expanded about
  !     XX, are stored in the array D. The expansion is of the form
  !     P(Z) = D(1) + D(2)*(Z-XX) +D(3)*((Z-XX)**2) + ... +
  !                                                  D(N)*((Z-XX)**(N-1)).
  !     Between the call to POLINT and the call to POLCOF the variable N
  !     and the arrays X and C must not be altered.
  !
  !     *****  INPUT PARAMETERS
  !
  !     XX   - The point about which the Taylor expansion is to be made.
  !
  !     N    - ****
  !            *     N, X, and C must remain unchanged between the
  !     X    - *     call to POLINT or the call to POLCOF.
  !     C    - ****
  !
  !     *****  OUTPUT PARAMETER
  !
  !     D    - The array of coefficients for the Taylor expansion as
  !            explained in the abstract
  !
  !     *****  STORAGE PARAMETER
  !
  !     WORK - This is an array to provide internal working storage. It
  !            must be dimensioned by at least 2*N in the calling program.
  !
  !
  !     **** Note - There are two methods for evaluating the fit produced
  !     by POLINT. You may call POLYVL to perform the task, or you may
  !     call POLCOF to obtain the coefficients of the Taylor expansion and
  !     then write your own evaluation scheme. Due to the inherent errors
  !     in the computations of the Taylor expansion from the Newton
  !     coefficients produced by POLINT, much more accuracy may be
  !     expected by calling POLYVL as opposed to writing your own scheme.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   890213  DATE WRITTEN
  !   891024  Corrected KEYWORD section.  (WRB)
  !   891024  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  POLCOF
  !
  DIMENSION X(*) , C(*) , D(*) , Work(*)
  !***FIRST EXECUTABLE STATEMENT  POLCOF
  DO k = 1 , N
    D(k) = C(k)
  ENDDO
  IF ( N==1 ) RETURN
  Work(1) = 1.0
  pone = C(1)
  nm1 = N - 1
  DO k = 2 , N
    km1 = k - 1
    npkm1 = N + k - 1
    Work(npkm1) = Xx - X(km1)
    Work(k) = Work(npkm1)*Work(km1)
    ptwo = pone + Work(k)*C(k)
    pone = ptwo
  ENDDO
  D(1) = ptwo
  IF ( N==2 ) RETURN
  DO k = 2 , nm1
    km1 = k - 1
    km2n = k - 2 + N
    nmkp1 = N - k + 1
    DO i = 2 , nmkp1
      km2npi = km2n + i
      im1 = i - 1
      km1pi = km1 + i
      Work(i) = Work(km2npi)*Work(im1) + Work(i)
      D(k) = D(k) + Work(i)*D(km1pi)
    ENDDO
  ENDDO
END SUBROUTINE POLCOF
