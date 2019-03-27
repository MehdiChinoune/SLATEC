!** DINTRV
SUBROUTINE DINTRV(Xt,Lxt,X,Ilo,Ileft,Mflag)
  IMPLICIT NONE
  !>
  !***
  !  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
  !            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
  !            the X interval.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E3, K6
  !***
  ! **Type:**      DOUBLE PRECISION (INTRV-S, DINTRV-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract    **** a double precision routine ****
  !         DINTRV is the INTERV routine of the reference.
  !
  !         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
  !         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
  !         the X interval.  Precisely,
  !
  !                      X .LT. XT(1)                1         -1
  !         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I , MFLAG=0
  !           XT(LXT) .LE. X                         LXT        1,
  !
  !         That is, when multiplicities are present in the break point
  !         to the left of X, the largest index is taken for ILEFT.
  !
  !     Description of Arguments
  !
  !         Input      XT,X are double precision
  !          XT      - XT is a knot or break point vector of length LXT
  !          LXT     - length of the XT vector
  !          X       - argument
  !          ILO     - an initialization parameter which must be set
  !                    to 1 the first time the spline array XT is
  !                    processed by DINTRV.
  !
  !         Output
  !          ILO     - ILO contains information for efficient process-
  !                    ing after the initial call and ILO must not be
  !                    changed by the user.  Distinct splines require
  !                    distinct ILO parameters.
  !          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
  !          MFLAG   - signals when X lies out of bounds
  !
  !     Error Conditions
  !         None
  !
  !***
  ! **References:**  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  !
  INTEGER ihi, Ileft, Ilo, istep, Lxt, Mflag, middle
  REAL(8) :: X, Xt(*)
  !* FIRST EXECUTABLE STATEMENT  DINTRV
  ihi = Ilo + 1
  IF ( ihi>=Lxt ) THEN
    IF ( X>=Xt(Lxt) ) GOTO 300
    IF ( Lxt<=1 ) GOTO 100
    Ilo = Lxt - 1
    ihi = Lxt
  ENDIF
  !
  IF ( X>=Xt(ihi) ) THEN
    !- ** NOW X .GE. XT(ILO) . FIND UPPER BOUND
    istep = 1
    DO
      Ilo = ihi
      ihi = Ilo + istep
      IF ( ihi>=Lxt ) THEN
        IF ( X>=Xt(Lxt) ) GOTO 300
        ihi = Lxt
        EXIT
      ELSE
        IF ( X<Xt(ihi) ) EXIT
        istep = istep*2
      ENDIF
    ENDDO
  ELSE
    IF ( X>=Xt(Ilo) ) GOTO 200
    !
    !- ** NOW X .LT. XT(IHI) . FIND LOWER BOUND
    istep = 1
    DO
      ihi = Ilo
      Ilo = ihi - istep
      IF ( Ilo<=1 ) THEN
        Ilo = 1
        IF ( X>=Xt(1) ) EXIT
        GOTO 100
      ELSE
        IF ( X>=Xt(Ilo) ) EXIT
        istep = istep*2
      ENDIF
    ENDDO
  ENDIF
  DO
    !
    !- ** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
    middle = (Ilo+ihi)/2
    IF ( middle==Ilo ) GOTO 200
    !     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
    IF ( X<Xt(middle) ) THEN
      ihi = middle
    ELSE
      Ilo = middle
    ENDIF
  ENDDO
  !- ** SET OUTPUT AND RETURN
  100  Mflag = -1
  Ileft = 1
  RETURN
  200  Mflag = 0
  Ileft = Ilo
  RETURN
  300  Mflag = 1
  Ileft = Lxt
END SUBROUTINE DINTRV
