!*==DXRED.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DXRED
SUBROUTINE DXRED(X,Ix,Ierror)
  IMPLICIT NONE
  !*--DXRED5
  !*** Start of declarations inserted by SPAG
  INTEGER i, Ierror, ixa, ixa1, ixa2
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DXRED
  !***PURPOSE  To provide double-precision floating-point arithmetic
  !            with an extended exponent range.
  !***LIBRARY   SLATEC
  !***CATEGORY  A3D
  !***TYPE      DOUBLE PRECISION (XRED-S, DXRED-D)
  !***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
  !***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
  !           Smith, John M., (NBS and George Mason University)
  !***DESCRIPTION
  !     DOUBLE PRECISION X
  !     INTEGER IX
  !
  !                  IF
  !                  RADIX**(-2L) .LE. (ABS(X),IX) .LE. RADIX**(2L)
  !                  THEN DXRED TRANSFORMS (X,IX) SO THAT IX=0.
  !                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE,
  !                  THEN DXRED TAKES NO ACTION.
  !                  THIS SUBROUTINE IS USEFUL IF THE
  !                  RESULTS OF EXTENDED-RANGE CALCULATIONS
  !                  ARE TO BE USED IN SUBSEQUENT ORDINARY
  !                  DOUBLE-PRECISION CALCULATIONS.
  !
  !***SEE ALSO  DXSET
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DXBLK2
  !***REVISION HISTORY  (YYMMDD)
  !   820712  DATE WRITTEN
  !   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  !***END PROLOGUE  DXRED
  REAL(8) :: X
  INTEGER Ix
  REAL(8) :: RADix, RADixl, RAD2l, DLG10r, xa
  INTEGER L, L2, KMAx
  COMMON /DXBLK2/ RADix, RADixl, RAD2l, DLG10r, L, L2, KMAx
  SAVE /DXBLK2/
  !
  !***FIRST EXECUTABLE STATEMENT  DXRED
  Ierror = 0
  IF ( X==0.0D0 ) THEN
    Ix = 0
  ELSE
    xa = ABS(X)
    IF ( Ix/=0 ) THEN
      ixa = ABS(Ix)
      ixa1 = ixa/L2
      ixa2 = MOD(ixa,L2)
      IF ( Ix>0 ) THEN
        !
        DO WHILE ( xa>=1.0D0 )
          xa = xa/RAD2l
          ixa1 = ixa1 + 1
        ENDDO
        xa = xa*RADix**ixa2
        IF ( ixa1/=0 ) THEN
          DO i = 1, ixa1
            IF ( xa>1.0D0 ) GOTO 99999
            xa = xa*RAD2l
          ENDDO
        ENDIF
      ELSE
        DO WHILE ( xa<=1.0D0 )
          xa = xa*RAD2l
          ixa1 = ixa1 + 1
        ENDDO
        xa = xa/RADix**ixa2
        IF ( ixa1/=0 ) THEN
          DO i = 1, ixa1
            IF ( xa<1.0D0 ) GOTO 99999
            xa = xa/RAD2l
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    IF ( xa<=RAD2l ) THEN
      IF ( xa>1.0D0 ) THEN
        X = SIGN(xa,X)
        Ix = 0
      ELSEIF ( RAD2l*xa>=1.0D0 ) THEN
        X = SIGN(xa,X)
        Ix = 0
      ENDIF
    ENDIF
  ENDIF
  99999 CONTINUE
  END SUBROUTINE DXRED
