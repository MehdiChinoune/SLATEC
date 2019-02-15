!DECK DPOCH
REAL(8) FUNCTION DPOCH(A,X)
  IMPLICIT NONE
  INTEGER i, ia, n
  !***BEGIN PROLOGUE  DPOCH
  !***PURPOSE  Evaluate a generalization of Pochhammer's symbol.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C1, C7A
  !***TYPE      DOUBLE PRECISION (POCH-S, DPOCH-D)
  !***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate a double precision generalization of Pochhammer's symbol
  ! (A)-sub-X = GAMMA(A+X)/GAMMA(A) for double precision A and X.
  ! For X a non-negative integer, POCH(A,X) is just Pochhammer's symbol.
  ! This is a preliminary version that does not handle wrong arguments
  ! properly and may not properly handle the case when the result is
  ! computed to less than half of double precision.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D9LGMC, DFAC, DGAMMA, DGAMR, DLGAMS, DLNREL, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  DPOCH
  REAL(8) :: A, X, absa, absax, alnga, alngax, ax, b, pi, &
    sgnga, sgngax, DFAC, DLNREL, D9LGMC, DGAMMA, &
    DGAMR, DCOT
  EXTERNAL DGAMMA
  SAVE pi
  DATA pi/3.141592653589793238462643383279503D0/
  !***FIRST EXECUTABLE STATEMENT  DPOCH
  ax = A + X
  IF ( ax<=0.0D0 ) THEN
    IF ( AINT(ax)==ax ) THEN
      !
      IF ( A>0.0D0.OR.AINT(A)/=A ) CALL XERMSG('SLATEC','DPOCH',&
        'A+X IS NON-POSITIVE INTEGER BUT A IS NOT',2,2)
      !
      ! WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
      !
      DPOCH = 1.0D0
      IF ( X==0.D0 ) RETURN
      !
      n = X
      IF ( MIN(A+X,A)<(-20.0D0) ) THEN
        !
        DPOCH = (-1.0D0)&
          **n*EXP((A-0.5D0)*DLNREL(X/(A-1.0D0))+X*LOG(-A+1.0D0-X)&
          -X+D9LGMC(-A+1.0D0)-D9LGMC(-A-X+1.D0))
        RETURN
      ELSE
        !
        ia = A
        DPOCH = (-1.0D0)**n*DFAC(-ia)/DFAC(-ia-n)
        RETURN
      ENDIF
    ENDIF
  ENDIF
  !
  ! A+X IS NOT ZERO OR A NEGATIVE INTEGER.
  !
  DPOCH = 0.0D0
  IF ( A<=0.0D0.AND.AINT(A)==A ) RETURN
  !
  n = ABS(X)
  IF ( REAL(n, 8)/=X.OR.n>20 ) THEN
    !
    absax = ABS(A+X)
    absa = ABS(A)
    IF ( MAX(absax,absa)<=20.0D0 ) THEN
      DPOCH = DGAMMA(A+X)*DGAMR(A)
      RETURN
      !
    ELSEIF ( ABS(X)>0.5D0*absa ) THEN
      !
      CALL DLGAMS(A+X,alngax,sgngax)
      CALL DLGAMS(A,alnga,sgnga)
      DPOCH = sgngax*sgnga*EXP(alngax-alnga)
      GOTO 99999
    ENDIF
  ELSE
    !
    ! X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
    !
    DPOCH = 1.0D0
    IF ( n==0 ) RETURN
    DO i = 1, n
      DPOCH = DPOCH*(A+i-1)
    ENDDO
    RETURN
  ENDIF
  !
  ! ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
  ! A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
  ! GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
  ! SIN(PI*A)/SIN(PI*(A+X))
  !
  b = A
  IF ( b<0.0D0 ) b = -A - X + 1.0D0
  DPOCH = EXP((b-0.5D0)*DLNREL(X/b)+X*LOG(b+X)-X+D9LGMC(b+X)-D9LGMC(b))
  IF ( A<0.0D0.AND.DPOCH/=0.0D0 )&
    DPOCH = DPOCH/(COS(pi*X)+DCOT(pi*A)*SIN(pi*X))
  RETURN
  !
  99999 CONTINUE
  END FUNCTION DPOCH
