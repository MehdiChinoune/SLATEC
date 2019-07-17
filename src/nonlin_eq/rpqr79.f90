!** RPQR79
PURE SUBROUTINE RPQR79(Ndeg,Coeff,Root,Ierr)
  !> Find the zeros of a polynomial with real coefficients.
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
  USE lapack, ONLY : SHSEQR

  INTEGER, INTENT(IN) :: Ndeg
  INTEGER, INTENT(OUT) :: Ierr
  REAL(SP), INTENT(IN) :: Coeff(Ndeg+1)
  COMPLEX(SP), INTENT(OUT) :: Root(Ndeg)
  INTEGER :: kwend, k, kh, kwr, kwi
  REAL(SP) :: scalee
  REAL(SP) :: z(1,Ndeg), h(Ndeg,Ndeg), wrkr(Ndeg), wrki(Ndeg), wrk(Ndeg)
  !* FIRST EXECUTABLE STATEMENT  RPQR79
  Ierr = 0
  IF( ABS(Coeff(1))==0.0 ) THEN
    Ierr = 2
    ERROR STOP 'RPQR79 : LEADING COEFFICIENT IS ZERO.'
    RETURN
  END IF
  !
  IF( Ndeg<=0 ) THEN
    Ierr = 3
    ERROR STOP 'RPQR79 : DEGREE INVALID.'
    RETURN
  END IF
  !
  IF( Ndeg==1 ) THEN
    Root(1) = CMPLX(-Coeff(2)/Coeff(1),0._SP,SP)
    RETURN
  END IF
  !
  scalee = 1._SP/Coeff(1)
  kh = 1
  kwr = kh + Ndeg*Ndeg
  kwi = kwr + Ndeg
  kwend = kwi + Ndeg - 1
  !
  h = 0._SP
  DO k = 1, Ndeg
    h(1,k) = -Coeff(k+1)*scalee
    IF( k==Ndeg ) Cycle
    h(k+1,k) = 1._SP
  END DO
  !
  CALL SHSEQR('E','N',Ndeg,1,Ndeg,h,Ndeg,wrkr,wrki,z,1,wrk,Ndeg,Ierr)
  !
  IF( Ierr/=0 ) THEN
    Ierr = 1
    ERROR STOP 'CPQR79 : NO CONVERGENCE IN 30 QR ITERATIONS.'
    RETURN
  END IF
  !
  Root(1:Ndeg) = CMPLX(wrkr(1:Ndeg),wrki(1:Ndeg),SP)

END SUBROUTINE RPQR79