!DECK SPLPFL
SUBROUTINE SPLPFL(Mrelas,Nvars,Ienter,Ileave,Ibasis,Ind,Ibb,Theta,Dirnrm,&
    Rprnrm,Csc,Ww,Bl,Bu,Erp,Rprim,Primal,Finite,Zerolv)
  IMPLICIT NONE
  INTEGER i, Ienter, Ileave, j, Mrelas, n20005, n20036, Nvars
  !***BEGIN PROLOGUE  SPLPFL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SPLPFL-S, DPLPFL-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/.
  !
  !     THIS SUBPROGRAM IS PART OF THE SPLP( ) PACKAGE.
  !     IT IMPLEMENTS THE PROCEDURE (CHOOSE VARIABLE TO LEAVE BASIS).
  !     REVISED 811130-1045
  !     REVISED YYMMDD-HHMM
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  SPLPFL
  INTEGER Ibasis(*), Ind(*), Ibb(*)
  REAL Csc(*), Ww(*), Bl(*), Bu(*), Erp(*), Rprim(*), Primal(*), &
    bound, Dirnrm, ratio, Rprnrm, Theta, zero
  LOGICAL Finite, Zerolv
  !***FIRST EXECUTABLE STATEMENT  SPLPFL
  zero = 0.E0
  !
  !     SEE IF THE ENTERING VARIABLE IS RESTRICTING THE STEP LENGTH
  !     BECAUSE OF AN UPPER BOUND.
  Finite = .FALSE.
  j = Ibasis(Ienter)
  IF ( Ind(j)==3 ) THEN
    Theta = Bu(j) - Bl(j)
    IF ( j<=Nvars ) Theta = Theta/Csc(j)
    Finite = .TRUE.
    Ileave = Ienter
  ENDIF
  !
  !     NOW USE THE BASIC VARIABLES TO POSSIBLY RESTRICT THE STEP
  !     LENGTH EVEN FURTHER.
  i = 1
  n20005 = Mrelas
  DO WHILE ( (n20005-i)>=0 )
    j = Ibasis(i)
    !
    !     IF THIS IS A FREE VARIABLE, DO NOT USE IT TO
    !     RESTRICT THE STEP LENGTH.
    IF ( Ind(j)==4 ) THEN
      i = i + 1
      !
      !     IF DIRECTION COMPONENT IS ABOUT ZERO, IGNORE IT FOR COMPUTING
      !     THE STEP LENGTH.
    ELSEIF ( ABS(Ww(i))<=Dirnrm*Erp(i) ) THEN
      i = i + 1
    ELSEIF ( Ww(i)<=zero ) THEN
      !
      !     IF THE VARIABLE IS LESS THAN ITS LOWER BOUND, IT CAN
      !     INCREASE ONLY TO ITS LOWER BOUND.
      IF ( Primal(i+Nvars)<zero ) THEN
        ratio = Rprim(i)/Ww(i)
        IF ( ratio<zero ) ratio = zero
        IF ( .NOT.Finite ) THEN
          Ileave = i
          Theta = ratio
          Finite = .TRUE.
        ELSEIF ( ratio<Theta ) THEN
          Ileave = i
          !
          !     IF THE BASIC VARIABLE IS FEASIBLE AND IS NOT AT ITS UPPER BOUND,
          !     THEN IT CAN INCREASE TO ITS UPPER BOUND.
          Theta = ratio
        ENDIF
      ELSEIF ( Ind(j)==3.AND.Primal(i+Nvars)==zero ) THEN
        bound = Bu(j) - Bl(j)
        IF ( j<=Nvars ) bound = bound/Csc(j)
        ratio = (bound-Rprim(i))/(-Ww(i))
        IF ( .NOT.Finite ) THEN
          Ileave = -i
          Theta = ratio
          Finite = .TRUE.
        ELSEIF ( ratio<Theta ) THEN
          Ileave = -i
          Theta = ratio
        ENDIF
      ENDIF
      i = i + 1
      !
      !     IF RPRIM(I) IS ESSENTIALLY ZERO, SET RATIO TO ZERO AND EXIT LOOP.
    ELSEIF ( ABS(Rprim(i))>Rprnrm*Erp(i) ) THEN
      !
      !     THE VALUE OF RPRIM(I) WILL DECREASE ONLY TO ITS LOWER BOUND OR
      !     ONLY TO ITS UPPER BOUND.  IF IT DECREASES TO ITS
      !     UPPER BOUND, THEN RPRIM(I) HAS ALREADY BEEN TRANSLATED
      !     TO ITS UPPER BOUND AND NOTHING NEEDS TO BE DONE TO IBB(J).
      IF ( Rprim(i)>zero ) THEN
        ratio = Rprim(i)/Ww(i)
        IF ( .NOT.Finite ) THEN
          Ileave = i
          Theta = ratio
          Finite = .TRUE.
        ELSEIF ( ratio<Theta ) THEN
          Ileave = i
          !
          !     THE VALUE RPRIM(I).LT.ZERO WILL NOT RESTRICT THE STEP.
          !
          !     THE DIRECTION COMPONENT IS NEGATIVE, THEREFORE THE VARIABLE WILL
          !     INCREASE.
          Theta = ratio
        ENDIF
      ENDIF
      i = i + 1
    ELSE
      Theta = zero
      Ileave = i
      Finite = .TRUE.
      EXIT
    ENDIF
  ENDDO
  !
  !     IF STEP LENGTH IS FINITE, SEE IF STEP LENGTH IS ABOUT ZERO.
  IF ( Finite ) THEN
    Zerolv = .TRUE.
    i = 1
    n20036 = Mrelas
    DO WHILE ( (n20036-i)>=0 )
      Zerolv = Zerolv .AND. ABS(Theta*Ww(i))<=Erp(i)*Rprnrm
      IF ( .NOT.Zerolv ) EXIT
      i = i + 1
    ENDDO
  ENDIF
END SUBROUTINE SPLPFL
