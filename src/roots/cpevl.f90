!** CPEVL
SUBROUTINE CPEVL(N,M,A,Z,C,B,Kbd)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CPZERO
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (CPEVL-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !        Evaluate a complex polynomial and its derivatives.
  !        Optionally compute error bounds for these values.
  !
  !   INPUT...
  !        N = Degree of the polynomial
  !        M = Number of derivatives to be calculated,
  !            M=0 evaluates only the function
  !            M=1 evaluates the function and first derivative, etc.
  !             if M .GT. N+1 function and all N derivatives will be
  !                calculated.
  !       A = Complex vector containing the N+1 coefficients of polynomial
  !               A(I)= coefficient of Z**(N+1-I)
  !        Z = Complex point at which the evaluation is to take place.
  !        C = Array of 2(M+1) words into which values are placed.
  !        B = Array of 2(M+1) words only needed if bounds are to be
  !              calculated.  It is not used otherwise.
  !        KBD = A logical variable, e.g. .TRUE. or .FALSE. which is
  !              to be set .TRUE. if bounds are to be computed.
  !
  !  OUTPUT...
  !        C =  C(I+1) contains the complex value of the I-th
  !              derivative at Z, I=0,...,M
  !        B =  B(I) contains the bounds on the real and imaginary parts
  !              of C(I) if they were requested.
  !
  !***
  ! **See also:**  CPZERO
  !***
  ! **Routines called:**  I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  REAL r, s
  INTEGER i, I1MACH, j, M, mini, N, np1
  COMPLEX A(*), C(*), Z, ci, cim1, B(*), bi, bim1, t
  LOGICAL Kbd
  REAL :: d1 = 0.0
  !* FIRST EXECUTABLE STATEMENT  CPEVL
  IF ( d1==0.0 ) d1 = REAL(I1MACH(10))**(1-I1MACH(11))
  np1 = N + 1
  DO j = 1, np1
    ci = 0.0
    cim1 = A(j)
    bi = 0.0
    bim1 = 0.0
    mini = MIN(M+1,N+2-j)
    DO i = 1, mini
      IF ( j/=1 ) ci = C(i)
      IF ( i/=1 ) cim1 = C(i-1)
      C(i) = cim1 + Z*ci
      IF ( Kbd ) THEN
        IF ( j/=1 ) bi = B(i)
        IF ( i/=1 ) bim1 = B(i-1)
        t = bi + (3.*d1+4.*d1*d1)*ZA(ci)
        r = REAL(ZA(Z)*CMPLX(REAL(t),-AIMAG(t)))
        s = AIMAG(ZA(Z)*t)
        B(i) = (1.+8.*d1)*(bim1+d1*ZA(cim1)+CMPLX(r,s))
        IF ( j==1 ) B(i) = 0.0
      ENDIF
    ENDDO
  ENDDO
CONTAINS
  COMPLEX FUNCTION ZA(q)
    COMPLEX, INTENT(IN) :: q
    ZA = CMPLX(ABS(REAL(q)),ABS(AIMAG(q)))
  END FUNCTION ZA
END SUBROUTINE CPEVL
