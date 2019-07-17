!** CPQR79
PURE SUBROUTINE CPQR79(Ndeg,Coeff,Root,Ierr)
  !> Find the zeros of a polynomial with complex coefficients.
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
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   911010  Code reworked and simplified.  (RWC and WRB)
  USE lapack, ONLY : CHSEQR

  INTEGER, INTENT(IN) :: Ndeg
  INTEGER, INTENT(OUT) :: Ierr
  COMPLEX(SP), INTENT(IN) :: Coeff(Ndeg+1)
  COMPLEX(SP), INTENT(OUT) :: Root(Ndeg)
  INTEGER :: k, khr, khi
  COMPLEX(SP) :: scalee
  COMPLEX(SP) :: h(Ndeg,Ndeg), z(1,Ndeg), wrk(Ndeg)
  !* FIRST EXECUTABLE STATEMENT  CPQR79
  Ierr = 0
  IF( ABS(Coeff(1))==0._SP ) THEN
    Ierr = 2
    ERROR STOP 'CPQR79 : LEADING COEFFICIENT IS ZERO.'
    RETURN
  END IF
  !
  IF( Ndeg<=0 ) THEN
    Ierr = 3
    ERROR STOP 'CPQR79 : DEGREE INVALID.'
    RETURN
  END IF
  !
  IF( Ndeg==1 ) THEN
    Root(1) = -Coeff(2)/Coeff(1)
    RETURN
  END IF
  !
  scalee = 1._SP/Coeff(1)
  khr = 1
  khi = khr + Ndeg*Ndeg
  !
  h = CMPLX( 0._SP, 0._SP, SP )
  DO k = 1, Ndeg
    h(1,k) = -scalee*Coeff(k+1)
    IF( k==Ndeg ) CYCLE
    h(k+1,k) = 1._SP
  END DO
  !
  CALL CHSEQR('E','N',Ndeg,1,Ndeg,h,Ndeg,Root,z,1,wrk,Ndeg,Ierr)
  !
  IF( Ierr/=0 ) THEN
    Ierr = 1
    ERROR STOP 'CPQR79 : NO CONVERGENCE IN 30 QR ITERATIONS.'
    RETURN
  END IF

END SUBROUTINE CPQR79