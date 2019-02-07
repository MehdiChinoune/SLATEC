!*==PCHCE.f90  processed by SPAG 6.72Dc at 11:00 on  6 Feb 2019
!DECK PCHCE
SUBROUTINE PCHCE(Ic,Vc,N,X,H,Slope,D,Incfd,Ierr)
  IMPLICIT NONE
  !*--PCHCE5
  !***BEGIN PROLOGUE  PCHCE
  !***SUBSIDIARY
  !***PURPOSE  Set boundary conditions for PCHIC
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      SINGLE PRECISION (PCHCE-S, DPCHCE-D)
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !          PCHCE:  PCHIC End Derivative Setter.
  !
  !    Called by PCHIC to set end derivatives as requested by the user.
  !    It must be called after interior derivative values have been set.
  !                      -----
  !
  !    To facilitate two-dimensional applications, includes an increment
  !    between successive values of the D-array.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  IC(2), N, IERR
  !        REAL  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
  !
  !        CALL  PCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
  !
  !   Parameters:
  !
  !     IC -- (input) integer array of length 2 specifying desired
  !           boundary conditions:
  !           IC(1) = IBEG, desired condition at beginning of data.
  !           IC(2) = IEND, desired condition at end of data.
  !           ( see prologue to PCHIC for details. )
  !
  !     VC -- (input) real array of length 2 specifying desired boundary
  !           values.  VC(1) need be set only if IC(1) = 2 or 3 .
  !                    VC(2) need be set only if IC(2) = 2 or 3 .
  !
  !     N -- (input) number of data points.  (assumes N.GE.2)
  !
  !     X -- (input) real array of independent variable values.  (the
  !           elements of X are assumed to be strictly increasing.)
  !
  !     H -- (input) real array of interval lengths.
  !     SLOPE -- (input) real array of data slopes.
  !           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
  !                  H(I) =  X(I+1)-X(I),
  !              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
  !
  !     D -- (input) real array of derivative values at the data points.
  !           The value corresponding to X(I) must be stored in
  !                D(1+(I-1)*INCFD),  I=1(1)N.
  !          (output) the value of D at X(1) and/or X(N) is changed, if
  !           necessary, to produce the requested boundary conditions.
  !           no other entries in D are changed.
  !
  !     INCFD -- (input) increment between successive values in D.
  !           This argument is provided primarily for 2-D applications.
  !
  !     IERR -- (output) error flag.
  !           Normal return:
  !              IERR = 0  (no errors).
  !           Warning errors:
  !              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
  !                        monotonicity.
  !              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
  !                        adjusted for monotonicity.
  !              IERR = 3  if both of the above are true.
  !
  !    -------
  !    WARNING:  This routine does no validity-checking of arguments.
  !    -------
  !
  !  Fortran intrinsics used:  ABS.
  !
  !***SEE ALSO  PCHIC
  !***ROUTINES CALLED  PCHDF, PCHST, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870707  Minor corrections made to prologue..
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR section in prologue.  (WRB)
  !   930503  Improved purpose.  (FNF)
  !***END PROLOGUE  PCHCE
  !
  !  Programming notes:
  !     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
  !        either argument is zero, +1 if they are of the same sign, and
  !        -1 if they are of opposite sign.
  !     2. One could reduce the number of arguments and amount of local
  !        storage, at the expense of reduced code clarity, by passing in
  !        the array WK (rather than splitting it into H and SLOPE) and
  !        increasing its length enough to incorporate STEMP and XTEMP.
  !     3. The two monotonicity checks only use the sufficient conditions.
  !        Thus, it is possible (but unlikely) for a boundary condition to
  !        be changed, even though the original interpolant was monotonic.
  !        (At least the result is a continuous function of the data.)
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER Ic(2), N, Incfd, Ierr
  REAL Vc(2), X(*), H(*), Slope(*), D(Incfd,*)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER ibeg, iend, ierf, index, j, k
  REAL half, stemp(3), three, two, xtemp(4), zero
  SAVE zero, half, two, three
  REAL PCHDF, PCHST
  !
  !  INITIALIZE.
  !
  DATA zero/0./, half/0.5/, two/2./, three/3./
  !
  !***FIRST EXECUTABLE STATEMENT  PCHCE
  ibeg = Ic(1)
  iend = Ic(2)
  Ierr = 0
  !
  !  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
  !
  IF ( ABS(ibeg)>N ) ibeg = 0
  IF ( ABS(iend)>N ) iend = 0
  !
  !  TREAT BEGINNING BOUNDARY CONDITION.
  !
  IF ( ibeg/=0 ) THEN
    k = ABS(ibeg)
    IF ( k==1 ) THEN
      !        BOUNDARY VALUE PROVIDED.
      D(1,1) = Vc(1)
    ELSEIF ( k==2 ) THEN
      !        BOUNDARY SECOND DERIVATIVE PROVIDED.
      D(1,1) = half*((three*Slope(1)-D(1,2))-half*Vc(1)*H(1))
    ELSEIF ( k<5 ) THEN
      !        USE K-POINT DERIVATIVE FORMULA.
      !        PICK UP FIRST K POINTS, IN REVERSE ORDER.
      DO j = 1, k
        index = k - j + 1
        !           INDEX RUNS FROM K DOWN TO 1.
        xtemp(j) = X(index)
        IF ( j<k ) stemp(j) = Slope(index-1)
      ENDDO
      !                 -----------------------------
      D(1,1) = PCHDF(k,xtemp,stemp,ierf)
      !                 -----------------------------
      IF ( ierf/=0 ) GOTO 100
    ELSE
      !        USE 'NOT A KNOT' CONDITION.
      D(1,1) = (three*(H(1)*Slope(2)+H(2)*Slope(1))-two*(H(1)+H(2))*D(1,2)&
        -H(1)*D(1,3))/H(2)
    ENDIF
    !
    IF ( ibeg<=0 ) THEN
      !
      !  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
      !
      IF ( Slope(1)==zero ) THEN
        IF ( D(1,1)/=zero ) THEN
          D(1,1) = zero
          Ierr = Ierr + 1
        ENDIF
      ELSEIF ( PCHST(D(1,1),Slope(1))<zero ) THEN
        D(1,1) = zero
        Ierr = Ierr + 1
      ELSEIF ( ABS(D(1,1))>three*ABS(Slope(1)) ) THEN
        D(1,1) = three*Slope(1)
        Ierr = Ierr + 1
      ENDIF
    ENDIF
  ENDIF
  !
  !  TREAT END BOUNDARY CONDITION.
  !
  IF ( iend/=0 ) THEN
    k = ABS(iend)
    IF ( k==1 ) THEN
      !        BOUNDARY VALUE PROVIDED.
      D(1,N) = Vc(2)
    ELSEIF ( k==2 ) THEN
      !        BOUNDARY SECOND DERIVATIVE PROVIDED.
      D(1,N) = half*((three*Slope(N-1)-D(1,N-1))+half*Vc(2)*H(N-1))
    ELSEIF ( k<5 ) THEN
      !        USE K-POINT DERIVATIVE FORMULA.
      !        PICK UP LAST K POINTS.
      DO j = 1, k
        index = N - k + j
        !           INDEX RUNS FROM N+1-K UP TO N.
        xtemp(j) = X(index)
        IF ( j<k ) stemp(j) = Slope(index)
      ENDDO
      !                 -----------------------------
      D(1,N) = PCHDF(k,xtemp,stemp,ierf)
      !                 -----------------------------
      IF ( ierf/=0 ) GOTO 100
    ELSE
      !        USE 'NOT A KNOT' CONDITION.
      D(1,N) = (three*(H(N-1)*Slope(N-2)+H(N-2)*Slope(N-1))&
        -two*(H(N-1)+H(N-2))*D(1,N-1)-H(N-1)*D(1,N-2))/H(N-2)
    ENDIF
    !
    IF ( iend<=0 ) THEN
      !
      !  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
      !
      IF ( Slope(N-1)==zero ) THEN
        IF ( D(1,N)/=zero ) THEN
          D(1,N) = zero
          Ierr = Ierr + 2
        ENDIF
      ELSEIF ( PCHST(D(1,N),Slope(N-1))<zero ) THEN
        D(1,N) = zero
        Ierr = Ierr + 2
      ELSEIF ( ABS(D(1,N))>three*ABS(Slope(N-1)) ) THEN
        D(1,N) = three*Slope(N-1)
        Ierr = Ierr + 2
      ENDIF
    ENDIF
  ENDIF
  !
  !  NORMAL RETURN.
  !
  RETURN
  !
  !  ERROR RETURN.
  !
  !     ERROR RETURN FROM PCHDF.
  !   *** THIS CASE SHOULD NEVER OCCUR ***
  100  Ierr = -1
  CALL XERMSG('SLATEC','PCHCE','ERROR RETURN FROM PCHDF',Ierr,1)
  !------------- LAST LINE OF PCHCE FOLLOWS ------------------------------
END SUBROUTINE PCHCE
