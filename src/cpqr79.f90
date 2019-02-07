!*==CPQR79.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CPQR79
SUBROUTINE CPQR79(Ndeg,Coeff,Root,Ierr,Work)
  IMPLICIT NONE
  !*--CPQR795
  !*** Start of declarations inserted by SPAG
  INTEGER km1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CPQR79
  !***PURPOSE  Find the zeros of a polynomial with complex coefficients.
  !***LIBRARY   SLATEC
  !***CATEGORY  F1A1B
  !***TYPE      COMPLEX (RPQR79-S, CPQR79-C)
  !***KEYWORDS  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
  !***AUTHOR  Vandevender, W. H., (SNLA)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  COMQR, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   791201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   911010  Code reworked and simplified.  (RWC and WRB)
  !***END PROLOGUE  CPQR79
  COMPLEX Coeff(*) , Root(*) , scale , c
  REAL Work(*)
  INTEGER Ndeg , Ierr , k , khr , khi , kwr , kwi , kad , kj
  !***FIRST EXECUTABLE STATEMENT  CPQR79
  Ierr = 0
  IF ( ABS(Coeff(1))==0.0 ) THEN
    Ierr = 2
    CALL XERMSG('SLATEC','CPQR79','LEADING COEFFICIENT IS ZERO.',2,1)
    RETURN
  ENDIF
  !
  IF ( Ndeg<=0 ) THEN
    Ierr = 3
    CALL XERMSG('SLATEC','CPQR79','DEGREE INVALID.',3,1)
    RETURN
  ENDIF
  !
  IF ( Ndeg==1 ) THEN
    Root(1) = -Coeff(2)/Coeff(1)
    RETURN
  ENDIF
  !
  scale = 1.0E0/Coeff(1)
  khr = 1
  khi = khr + Ndeg*Ndeg
  kwr = khi + khi - khr
  kwi = kwr + Ndeg
  !
  DO k = 1 , kwr
    Work(k) = 0.0E0
  ENDDO
  !
  DO k = 1 , Ndeg
    kad = (k-1)*Ndeg + 1
    c = scale*Coeff(k+1)
    Work(kad) = -REAL(c)
    kj = khi + kad - 1
    Work(kj) = -AIMAG(c)
    IF ( k/=Ndeg ) Work(kad+k) = 1.0E0
  ENDDO
  !
  CALL COMQR(Ndeg,Ndeg,1,Ndeg,Work(khr),Work(khi),Work(kwr),Work(kwi),Ierr)
  !
  IF ( Ierr/=0 ) THEN
    Ierr = 1
    CALL XERMSG('SLATEC','CPQR79','NO CONVERGENCE IN 30 QR ITERATIONS.',1,1)
    RETURN
  ENDIF
  !
  DO k = 1 , Ndeg
    km1 = k - 1
    Root(k) = CMPLX(Work(kwr+km1),Work(kwi+km1))
  ENDDO
END SUBROUTINE CPQR79
