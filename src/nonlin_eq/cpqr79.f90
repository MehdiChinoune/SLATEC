!** CPQR79
SUBROUTINE CPQR79(Ndeg,Coeff,Root,Ierr,Work)
  !>
  !  Find the zeros of a polynomial with complex coefficients.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F1A1B
  !***
  ! **Type:**      COMPLEX (RPQR79-S, CPQR79-C)
  !***
  ! **Keywords:**  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
  !***
  ! **Author:**  Vandevender, W. H., (SNLA)
  !***
  ! **Description:**
  !
  !   Abstract
  !       This routine computes all zeros of a polynomial of degree NDEG
  !       with complex coefficients by computing the eigenvalues of the
  !       companion matrix.
  !
  !   Description of Parameters
  !       The user must dimension all arrays appearing in the call list
  !            COEFF(NDEG+1), ROOT(NDEG), WORK(2*NDEG*(NDEG+1))
  !
  !    --Input--
  !      NDEG    degree of polynomial
  !
  !      COEFF   COMPLEX coefficients in descending order.  i.e.,
  !              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
  !
  !      WORK    REAL work array of dimension at least 2*NDEG*(NDEG+1)
  !
  !   --Output--
  !      ROOT    COMPLEX vector of roots
  !
  !      IERR    Output Error Code
  !           - Normal Code
  !          0  means the roots were computed.
  !           - Abnormal Codes
  !          1  more than 30 QR iterations on some eigenvalue of the
  !             companion matrix
  !          2  COEFF(1)=0.0
  !          3  NDEG is invalid (less than or equal to 0)
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  COMQR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   791201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   911010  Code reworked and simplified.  (RWC and WRB)
  USE service, ONLY : XERMSG
  USE linear, ONLY : COMQR
  INTEGER km1
  COMPLEX Coeff(*), Root(*), scalee, c
  REAL Work(*)
  INTEGER Ndeg, Ierr, k, khr, khi, kwr, kwi, kad, kj
  !* FIRST EXECUTABLE STATEMENT  CPQR79
  Ierr = 0
  IF ( ABS(Coeff(1))==0.0 ) THEN
    Ierr = 2
    CALL XERMSG('SLATEC','CPQR79','LEADING COEFFICIENT IS ZERO.',2,1)
    RETURN
  END IF
  !
  IF ( Ndeg<=0 ) THEN
    Ierr = 3
    CALL XERMSG('SLATEC','CPQR79','DEGREE INVALID.',3,1)
    RETURN
  END IF
  !
  IF ( Ndeg==1 ) THEN
    Root(1) = -Coeff(2)/Coeff(1)
    RETURN
  END IF
  !
  scalee = 1.0E0/Coeff(1)
  khr = 1
  khi = khr + Ndeg*Ndeg
  kwr = khi + khi - khr
  kwi = kwr + Ndeg
  !
  DO k = 1, kwr
    Work(k) = 0.0E0
  END DO
  !
  DO k = 1, Ndeg
    kad = (k-1)*Ndeg + 1
    c = scalee*Coeff(k+1)
    Work(kad) = -REAL(c)
    kj = khi + kad - 1
    Work(kj) = -AIMAG(c)
    IF ( k/=Ndeg ) Work(kad+k) = 1.0E0
  END DO
  !
  CALL COMQR(Ndeg,Ndeg,1,Ndeg,Work(khr),Work(khi),Work(kwr),Work(kwi),Ierr)
  !
  IF ( Ierr/=0 ) THEN
    Ierr = 1
    CALL XERMSG('SLATEC','CPQR79','NO CONVERGENCE IN 30 QR ITERATIONS.',1,1)
    RETURN
  END IF
  !
  DO k = 1, Ndeg
    km1 = k - 1
    Root(k) = CMPLX(Work(kwr+km1),Work(kwi+km1))
  END DO
END SUBROUTINE CPQR79
