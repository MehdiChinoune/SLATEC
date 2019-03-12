!DECK PCHIC
SUBROUTINE PCHIC(Ic,Vc,Switch,N,X,F,D,Incfd,Wk,Nwk,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PCHIC
  !***PURPOSE  Set derivatives needed to determine a piecewise monotone
  !            piecewise cubic Hermite interpolant to given data.
  !            User control is available over boundary conditions and/or
  !            treatment of points where monotonicity switches direction.
  !***LIBRARY   SLATEC (PCHIP)
  !***CATEGORY  E1A
  !***TYPE      SINGLE PRECISION (PCHIC-S, DPCHIC-D)
  !***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
  !             PCHIP, PIECEWISE CUBIC INTERPOLATION,
  !             SHAPE-PRESERVING INTERPOLATION
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***DESCRIPTION
  !
  !         PCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
  !
  !     Sets derivatives needed to determine a piecewise monotone piece-
  !     wise cubic interpolant to the data given in X and F satisfying the
  !     boundary conditions specified by IC and VC.
  !
  !     The treatment of points where monotonicity switches direction is
  !     controlled by argument SWITCH.
  !
  !     To facilitate two-dimensional applications, includes an increment
  !     between successive values of the F- and D-arrays.
  !
  !     The resulting piecewise cubic Hermite function may be evaluated
  !     by PCHFE or PCHFD.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  IC(2), N, NWK, IERR
  !        REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
  !
  !        CALL  PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
  !
  !   Parameters:
  !
  !     IC -- (input) integer array of length 2 specifying desired
  !           boundary conditions:
  !           IC(1) = IBEG, desired condition at beginning of data.
  !           IC(2) = IEND, desired condition at end of data.
  !
  !           IBEG = 0  for the default boundary condition (the same as
  !                     used by PCHIM).
  !           If IBEG.NE.0, then its sign indicates whether the boundary
  !                     derivative is to be adjusted, if necessary, to be
  !                     compatible with monotonicity:
  !              IBEG.GT.0  if no adjustment is to be performed.
  !              IBEG.LT.0  if the derivative is to be adjusted for
  !                     monotonicity.
  !
  !           Allowable values for the magnitude of IBEG are:
  !           IBEG = 1  if first derivative at X(1) is given in VC(1).
  !           IBEG = 2  if second derivative at X(1) is given in VC(1).
  !           IBEG = 3  to use the 3-point difference formula for D(1).
  !                     (Reverts to the default b.c. if N.LT.3 .)
  !           IBEG = 4  to use the 4-point difference formula for D(1).
  !                     (Reverts to the default b.c. if N.LT.4 .)
  !           IBEG = 5  to set D(1) so that the second derivative is con-
  !              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
  !              This option is somewhat analogous to the "not a knot"
  !              boundary condition provided by PCHSP.
  !
  !          NOTES (IBEG):
  !           1. An error return is taken if ABS(IBEG).GT.5 .
  !           2. Only in case  IBEG.LE.0  is it guaranteed that the
  !              interpolant will be monotonic in the first interval.
  !              If the returned value of D(1) lies between zero and
  !              3*SLOPE(1), the interpolant will be monotonic.  This
  !              is **NOT** checked if IBEG.GT.0 .
  !           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
  !              tonicity, a warning error is returned.
  !
  !           IEND may take on the same values as IBEG, but applied to
  !           derivative at X(N).  In case IEND = 1 or 2, the value is
  !           given in VC(2).
  !
  !          NOTES (IEND):
  !           1. An error return is taken if ABS(IEND).GT.5 .
  !           2. Only in case  IEND.LE.0  is it guaranteed that the
  !              interpolant will be monotonic in the last interval.
  !              If the returned value of D(1+(N-1)*INCFD) lies between
  !              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
  !              This is **NOT** checked if IEND.GT.0 .
  !           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
  !              achieve monotonicity, a warning error is returned.
  !
  !     VC -- (input) real array of length 2 specifying desired boundary
  !           values, as indicated above.
  !           VC(1) need be set only if IC(1) = 1 or 2 .
  !           VC(2) need be set only if IC(2) = 1 or 2 .
  !
  !     SWITCH -- (input) indicates desired treatment of points where
  !           direction of monotonicity switches:
  !           Set SWITCH to zero if interpolant is required to be mono-
  !           tonic in each interval, regardless of monotonicity of data.
  !             NOTES:
  !              1. This will cause D to be set to zero at all switch
  !                 points, thus forcing extrema there.
  !              2. The result of using this option with the default boun-
  !                 dary conditions will be identical to using PCHIM, but
  !                 will generally cost more compute time.
  !                 This option is provided only to facilitate comparison
  !                 of different switch and/or boundary conditions.
  !           Set SWITCH nonzero to use a formula based on the 3-point
  !              difference formula in the vicinity of switch points.
  !           If SWITCH is positive, the interpolant on each interval
  !              containing an extremum is controlled to not deviate from
  !              the data by more than SWITCH*DFLOC, where DFLOC is the
  !              maximum of the change of F on this interval and its two
  !              immediate neighbors.
  !           If SWITCH is negative, no such control is to be imposed.
  !
  !     N -- (input) number of data points.  (Error return if N.LT.2 .)
  !
  !     X -- (input) real array of independent variable values.  The
  !           elements of X must be strictly increasing:
  !                X(I-1) .LT. X(I),  I = 2(1)N.
  !           (Error return if not.)
  !
  !     F -- (input) real array of dependent variable values to be inter-
  !           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
  !
  !     D -- (output) real array of derivative values at the data points.
  !           These values will determine a monotone cubic Hermite func-
  !           tion on each subinterval on which the data are monotonic,
  !           except possibly adjacent to switches in monotonicity.
  !           The value corresponding to X(I) is stored in
  !                D(1+(I-1)*INCFD),  I=1(1)N.
  !           No other entries in D are changed.
  !
  !     INCFD -- (input) increment between successive values in F and D.
  !           This argument is provided primarily for 2-D applications.
  !           (Error return if  INCFD.LT.1 .)
  !
  !     WK -- (scratch) real array of working storage.  The user may wish
  !           to know that the returned values are:
  !              WK(I)     = H(I)     = X(I+1) - X(I) ;
  !              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
  !           for  I = 1(1)N-1.
  !
  !     NWK -- (input) length of work array.
  !           (Error return if  NWK.LT.2*(N-1) .)
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
  !           "Recoverable" errors:
  !              IERR = -1  if N.LT.2 .
  !              IERR = -2  if INCFD.LT.1 .
  !              IERR = -3  if the X-array is not strictly increasing.
  !              IERR = -4  if ABS(IBEG).GT.5 .
  !              IERR = -5  if ABS(IEND).GT.5 .
  !              IERR = -6  if both of the above are true.
  !              IERR = -7  if NWK.LT.2*(N-1) .
  !             (The D-array has not been changed in any of these cases.)
  !               NOTE:  The above errors are checked in the order listed,
  !                   and following arguments have **NOT** been validated.
  !
  !***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
  !                 Package, Report UCRL-87285, Lawrence Livermore Nation-
  !                 al Laboratory, July 1982.  [Poster presented at the
  !                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
  !               2. F. N. Fritsch and J. Butland, A method for construc-
  !                 ting local monotone piecewise cubic interpolants, SIAM
  !                 Journal on Scientific and Statistical Computing 5, 2
  !                 (June 1984), pp. 300-304.
  !               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
  !                 cubic interpolation, SIAM Journal on Numerical Ana-
  !                 lysis 17, 2 (April 1980), pp. 238-246.
  !***ROUTINES CALLED  PCHCE, PCHCI, PCHCS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820804  Converted to SLATEC library version.
  !   870813  Updated Reference 2.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890703  Corrected category record.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920429  Revised format and order of references.  (WRB,FNF)
  !***END PROLOGUE  PCHIC
  !  Programming notes:
  !
  !     To produce a double precision version, simply:
  !        a. Change PCHIC to DPCHIC wherever it occurs,
  !        b. Change PCHCE to DPCHCE wherever it occurs,
  !        c. Change PCHCI to DPCHCI wherever it occurs,
  !        d. Change PCHCS to DPCHCS wherever it occurs,
  !        e. Change the real declarations to double precision, and
  !        f. Change the constant  ZERO  to double precision.
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER Ic(2), N, Incfd, Nwk, Ierr
  REAL Vc(2), Switch, X(*), F(Incfd,*), D(Incfd,*), Wk(Nwk)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER i, ibeg, iend, nless1
  REAL zero
  SAVE zero
  DATA zero/0./
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  !***FIRST EXECUTABLE STATEMENT  PCHIC
  IF ( N<2 ) THEN
    !
    !  ERROR RETURNS.
    !
    !     N.LT.2 RETURN.
    Ierr = -1
    CALL XERMSG('SLATEC','PCHIC','NUMBER OF DATA POINTS LESS THAN TWO',Ierr,&
      1)
    RETURN
  ELSEIF ( Incfd<1 ) THEN
    !
    !     INCFD.LT.1 RETURN.
    Ierr = -2
    CALL XERMSG('SLATEC','PCHIC','INCREMENT LESS THAN ONE',Ierr,1)
    RETURN
  ELSE
    DO i = 2, N
      IF ( X(i)<=X(i-1) ) GOTO 100
    ENDDO
    !
    ibeg = Ic(1)
    iend = Ic(2)
    Ierr = 0
    IF ( ABS(ibeg)>5 ) Ierr = Ierr - 1
    IF ( ABS(iend)>5 ) Ierr = Ierr - 2
    IF ( Ierr<0 ) THEN
      !
      !     IC OUT OF RANGE RETURN.
      Ierr = Ierr - 3
      CALL XERMSG('SLATEC','PCHIC','IC OUT OF RANGE',Ierr,1)
      RETURN
    ELSE
      !
      !  FUNCTION DEFINITION IS OK -- GO ON.
      !
      nless1 = N - 1
      IF ( Nwk<2*nless1 ) THEN
        !
        !     NWK .LT. 2*(N-1)  RETURN.
        Ierr = -7
        CALL XERMSG('SLATEC','PCHIC','WORK ARRAY TOO SMALL',Ierr,1)
        RETURN
      ELSE
        !
        !  SET UP H AND SLOPE ARRAYS.
        !
        DO i = 1, nless1
          Wk(i) = X(i+1) - X(i)
          Wk(nless1+i) = (F(1,i+1)-F(1,i))/Wk(i)
        ENDDO
        !
        !  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
        !
        IF ( nless1>1 ) THEN
          !
          !  NORMAL CASE  (N .GE. 3) .
          !
          !
          !  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
          !
          !     --------------------------------------
          CALL PCHCI(N,Wk(1),Wk(N),D,Incfd)
          !     --------------------------------------
          !
          !  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
          !
          IF ( Switch/=zero ) THEN
            !     ----------------------------------------------------
            CALL PCHCS(Switch,N,Wk(1),Wk(N),D,Incfd,Ierr)
            !     ----------------------------------------------------
            IF ( Ierr/=0 ) THEN
              !
              !     ERROR RETURN FROM PCHCS.
              Ierr = -8
              CALL XERMSG('SLATEC','PCHIC','ERROR RETURN FROM PCHCS',Ierr,1)
              RETURN
            ENDIF
          ENDIF
        ELSE
          D(1,1) = Wk(2)
          D(1,N) = Wk(2)
        ENDIF
        !
        !  SET END CONDITIONS.
        !
        IF ( (ibeg/=0).OR.(iend/=0) ) THEN
          !     -------------------------------------------------------
          CALL PCHCE(Ic,Vc,N,X,Wk(1),Wk(N),D,Incfd,Ierr)
          !     -------------------------------------------------------
          IF ( Ierr<0 ) THEN
            !
            !     ERROR RETURN FROM PCHCE.
            !   *** THIS CASE SHOULD NEVER OCCUR ***
            Ierr = -9
            CALL XERMSG('SLATEC','PCHIC','ERROR RETURN FROM PCHCE',Ierr,1)
            RETURN
          ENDIF
        ENDIF
        !
        !  NORMAL RETURN.
        !
        RETURN
      ENDIF
    ENDIF
  ENDIF
  !
  !     X-ARRAY NOT STRICTLY INCREASING.
  100  Ierr = -3
  CALL XERMSG('SLATEC','PCHIC','X-ARRAY NOT STRICTLY INCREASING',Ierr,1)
  RETURN
  !------------- LAST LINE OF PCHIC FOLLOWS ------------------------------
  RETURN
END SUBROUTINE PCHIC
