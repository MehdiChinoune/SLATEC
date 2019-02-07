!*==PPVAL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PPVAL
FUNCTION PPVAL(Ldc,C,Xi,Lxi,K,Ideriv,X,Inppv)
  IMPLICIT NONE
  !*--PPVAL5
  !*** Start of declarations inserted by SPAG
  REAL PPVAL
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PPVAL
  !***PURPOSE  Calculate the value of the IDERIV-th derivative of the
  !            B-spline from the PP-representation.
  !***LIBRARY   SLATEC
  !***CATEGORY  E3, K6
  !***TYPE      SINGLE PRECISION (PPVAL-S, DPPVAL-D)
  !***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract
  !         PPVAL is the PPVALU function of the reference.
  !
  !         PPVAL calculates (at X) the value of the IDERIV-th
  !         derivative of the B-spline from the PP-representation
  !         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in
  !         the interval XI(J) .LE. X .LT. XI(J+1) is evaluated, J=1,LXI.
  !         Right limiting values at X=XI(J) are obtained.  PPVAL will
  !         extrapolate beyond XI(1) and XI(LXI+1).
  !
  !         To obtain left limiting values (left derivatives) at XI(J),
  !         replace LXI by J-1 and set X=XI(J),J=2,LXI+1.
  !
  !     Description of Arguments
  !         Input
  !          LDC     - leading dimension of C matrix, LDC .GE. K
  !          C       - matrix of dimension at least (K,LXI) containing
  !                    right derivatives at break points XI(*).
  !          XI      - break point vector of length LXI+1
  !          LXI     - number of polynomial pieces
  !          K       - order of B-spline, K .GE. 1
  !          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
  !                    IDERIV=0 gives the B-spline value
  !          X       - argument, XI(1) .LE. X .LE. XI(LXI+1)
  !          INPPV   - an initialization parameter which must be set
  !                    to 1 the first time PPVAL is called.
  !
  !         Output
  !          INPPV   - INPPV contains information for efficient process-
  !                    ing after the initial call and INPPV must not
  !                    be changed by the user.  Distinct splines require
  !                    distinct INPPV parameters.
  !          PPVAL   - value of the IDERIV-th derivative at X
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***REFERENCES  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***ROUTINES CALLED  INTRV, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  PPVAL
  !
  INTEGER i , Ideriv , Inppv , j , K , Ldc , Lxi , ndummy
  REAL C , dx , fltk , X , Xi
  DIMENSION Xi(*) , C(Ldc,*)
  !***FIRST EXECUTABLE STATEMENT  PPVAL
  PPVAL = 0.0E0
  IF ( K<1 ) THEN
    CALL XERMSG('SLATEC','PPVAL','K DOES NOT SATISFY K.GE.1',2,1)
    RETURN
  ELSEIF ( Ldc<K ) THEN
    !
    !
    CALL XERMSG('SLATEC','PPVAL','LDC DOES NOT SATISFY LDC.GE.K',2,1)
    RETURN
  ELSEIF ( Lxi<1 ) THEN
    CALL XERMSG('SLATEC','PPVAL','LXI DOES NOT SATISFY LXI.GE.1',2,1)
    RETURN
  ELSEIF ( Ideriv<0.OR.Ideriv>=K ) THEN
    CALL XERMSG('SLATEC','PPVAL','IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K',&
      2,1)
    GOTO 99999
  ELSE
    i = K - Ideriv
    fltk = i
    CALL INTRV(Xi,Lxi,X,Inppv,i,ndummy)
    dx = X - Xi(i)
    j = K
    DO
      PPVAL = (PPVAL/fltk)*dx + C(j,i)
      j = j - 1
      fltk = fltk - 1.0E0
      IF ( fltk<=0.0E0 ) EXIT
    ENDDO
  ENDIF
  RETURN
  99999 END FUNCTION PPVAL
