!** RPQR79
SUBROUTINE RPQR79(Ndeg,Coeff,Root,Ierr,Work)
  !>
  !  Find the zeros of a polynomial with real coefficients.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F1A1A
  !***
  ! **Type:**      SINGLE PRECISION (RPQR79-S, CPQR79-C)
  !***
  ! **Keywords:**  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
  !***
  ! **Author:**  Vandevender, W. H., (SNLA)
  !***
  ! **Description:**
  !
  !   Abstract
  !       This routine computes all zeros of a polynomial of degree NDEG
  !       with real coefficients by computing the eigenvalues of the
  !       companion matrix.
  !
  !   Description of Parameters
  !       The user must dimension all arrays appearing in the call list
  !            COEFF(NDEG+1), ROOT(NDEG), WORK(NDEG*(NDEG+2))
  !
  !    --Input--
  !      NDEG    degree of polynomial
  !
  !      COEFF   REAL coefficients in descending order.  i.e.,
  !              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
  !
  !      WORK    REAL work array of dimension at least NDEG*(NDEG+2)
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
  ! **Routines called:**  HQR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800601  DATE WRITTEN
  !   890505  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   911010  Code reworked and simplified.  (RWC and WRB)
  USE service, ONLY : XERMSG
  USE linear, ONLY : HQR
  INTEGER :: Ndeg, Ierr
  REAL :: Coeff(Ndeg+1), Work(Ndeg*(Ndeg+2))
  COMPLEX :: Root(Ndeg)
  INTEGER :: km1, kwend, k, kh, kwr, kwi, kcol
  REAL :: scalee
  !* FIRST EXECUTABLE STATEMENT  RPQR79
  Ierr = 0
  IF ( ABS(Coeff(1))==0.0 ) THEN
    Ierr = 2
    CALL XERMSG('SLATEC','RPQR79','LEADING COEFFICIENT IS ZERO.',2,1)
    RETURN
  END IF
  !
  IF ( Ndeg<=0 ) THEN
    Ierr = 3
    CALL XERMSG('SLATEC','RPQR79','DEGREE INVALID.',3,1)
    RETURN
  END IF
  !
  IF ( Ndeg==1 ) THEN
    Root(1) = CMPLX(-Coeff(2)/Coeff(1),0.0)
    RETURN
  END IF
  !
  scalee = 1.0E0/Coeff(1)
  kh = 1
  kwr = kh + Ndeg*Ndeg
  kwi = kwr + Ndeg
  kwend = kwi + Ndeg - 1
  !
  DO k = 1, kwend
    Work(k) = 0.0E0
  END DO
  !
  DO k = 1, Ndeg
    kcol = (k-1)*Ndeg + 1
    Work(kcol) = -Coeff(k+1)*scalee
    IF ( k/=Ndeg ) Work(kcol+k) = 1.0E0
  END DO
  !
  CALL HQR(Ndeg,Ndeg,1,Ndeg,Work(kh),Work(kwr),Work(kwi),Ierr)
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
END SUBROUTINE RPQR79
