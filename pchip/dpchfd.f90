!DECK DPCHFD
SUBROUTINE DPCHFD(N,X,F,D,Incfd,Skip,Ne,Xe,Fe,De,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DPCHFD
  !***PURPOSE  Evaluate a piecewise cubic Hermite function and its first
  !            derivative at an array of points.  May be used by itself
  !            for Hermite interpolation, or as an evaluator for DPCHIM
  !            or DPCHIC. If only function values are required, use
  !            DPCHFE instead.
  !***LIBRARY   SLATEC (PCHIP)
  !***CATEGORY  E3, H1
  !***TYPE      DOUBLE PRECISION (PCHFD-S, DPCHFD-D)
  !***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
  !             HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC EVALUATION
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***DESCRIPTION
  !
  !          DPCHFD:  Piecewise Cubic Hermite Function and Derivative
  !                  evaluator
  !
  !     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
  !     gether with its first derivative, at the points  XE(J), J=1(1)NE.
  !
  !     If only function values are required, use DPCHFE, instead.
  !
  !     To provide compatibility with DPCHIM and DPCHIC, includes an
  !     increment between successive values of the F- and D-arrays.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  N, NE, IERR
  !        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
  !                          DE(NE)
  !        LOGICAL  SKIP
  !
  !        CALL  DPCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
  !
  !   Parameters:
  !
  !     N -- (input) number of data points.  (Error return if N.LT.2 .)
  !
  !     X -- (input) real*8 array of independent variable values.  The
  !           elements of X must be strictly increasing:
  !                X(I-1) .LT. X(I),  I = 2(1)N.
  !           (Error return if not.)
  !
  !     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
  !           the value corresponding to X(I).
  !
  !     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
  !           is the value corresponding to X(I).
  !
  !     INCFD -- (input) increment between successive values in F and D.
  !           (Error return if  INCFD.LT.1 .)
  !
  !     SKIP -- (input/output) logical variable which should be set to
  !           .TRUE. if the user wishes to skip checks for validity of
  !           preceding parameters, or to .FALSE. otherwise.
  !           This will save time in case these checks have already
  !           been performed (say, in DPCHIM or DPCHIC).
  !           SKIP will be set to .TRUE. on normal return.
  !
  !     NE -- (input) number of evaluation points.  (Error return if
  !           NE.LT.1 .)
  !
  !     XE -- (input) real*8 array of points at which the functions are to
  !           be evaluated.
  !
  !
  !          NOTES:
  !           1. The evaluation will be most efficient if the elements
  !              of XE are increasing relative to X;
  !              that is,   XE(J) .GE. X(I)
  !              implies    XE(K) .GE. X(I),  all K.GE.J .
  !           2. If any of the XE are outside the interval [X(1),X(N)],
  !              values are extrapolated from the nearest extreme cubic,
  !              and a warning error is returned.
  !
  !     FE -- (output) real*8 array of values of the cubic Hermite
  !           function defined by  N, X, F, D  at the points  XE.
  !
  !     DE -- (output) real*8 array of values of the first derivative of
  !           the same function at the points  XE.
  !
  !     IERR -- (output) error flag.
  !           Normal return:
  !              IERR = 0  (no errors).
  !           Warning error:
  !              IERR.GT.0  means that extrapolation was performed at
  !                 IERR points.
  !           "Recoverable" errors:
  !              IERR = -1  if N.LT.2 .
  !              IERR = -2  if INCFD.LT.1 .
  !              IERR = -3  if the X-array is not strictly increasing.
  !              IERR = -4  if NE.LT.1 .
  !           (Output arrays have not been changed in any of these cases.)
  !               NOTE:  The above errors are checked in the order listed,
  !                   and following arguments have **NOT** been validated.
  !              IERR = -5  if an error has occurred in the lower-level
  !                         routine DCHFDV.  NB: this should never happen.
  !                         Notify the author **IMMEDIATELY** if it does.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DCHFDV, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811020  DATE WRITTEN
  !   820803  Minor cosmetic changes for release 1.
  !   870707  Corrected XERROR calls for d.p. name(s).
  !   890206  Corrected XERROR calls.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DPCHFD
  !  Programming notes:
  !
  !     1. To produce a single precision version, simply:
  !        a. Change DPCHFD to PCHFD, and DCHFDV to CHFDV, wherever they
  !           occur,
  !        b. Change the double precision declaration to real,
  !
  !     2. Most of the coding between the call to DCHFDV and the end of
  !        the IR-loop could be eliminated if it were permissible to
  !        assume that XE is ordered relative to X.
  !
  !     3. DCHFDV does not assume that X1 is less than X2.  thus, it would
  !        be possible to write a version of DPCHFD that assumes a strict-
  !        ly decreasing X-array by simply running the IR-loop backwards
  !        (and reversing the order of appropriate tests).
  !
  !     4. The present code has a minor bug, which I have decided is not
  !        worth the effort that would be required to fix it.
  !        If XE contains points in [X(N-1),X(N)], followed by points .LT.
  !        X(N-1), followed by points .GT.X(N), the extrapolation points
  !        will be counted (at least) twice in the total returned in IERR.
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER N, Incfd, Ne, Ierr
  REAL(8) :: X(*), F(Incfd,*), D(Incfd,*), Xe(*), Fe(*), De(*)
  LOGICAL Skip
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER i, ierc, ir, j, jfirst, next(2), nj
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHFD
  IF ( .NOT.(Skip) ) THEN
    !
    IF ( N<2 ) THEN
      !
      !  ERROR RETURNS.
      !
      !     N.LT.2 RETURN.
      Ierr = -1
      CALL XERMSG('SLATEC','DPCHFD','NUMBER OF DATA POINTS LESS THAN TWO',&
        Ierr,1)
      RETURN
    ELSEIF ( Incfd<1 ) THEN
      !
      !     INCFD.LT.1 RETURN.
      Ierr = -2
      CALL XERMSG('SLATEC','DPCHFD','INCREMENT LESS THAN ONE',Ierr,1)
      RETURN
    ELSE
      DO i = 2, N
        IF ( X(i)<=X(i-1) ) GOTO 500
      ENDDO
    ENDIF
  ENDIF
  !
  !  FUNCTION DEFINITION IS OK, GO ON.
  !
  IF ( Ne<1 ) THEN
    !
    !     NE.LT.1 RETURN.
    Ierr = -4
    CALL XERMSG('SLATEC','DPCHFD',&
      'NUMBER OF EVALUATION POINTS LESS THAN ONE',Ierr,1)
    RETURN
  ELSE
    Ierr = 0
    Skip = .TRUE.
    !
    !  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
    !                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
    jfirst = 1
    ir = 2
  ENDIF
  !
  !     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
  !
  100 CONTINUE
  IF ( jfirst>Ne ) GOTO 400
  !
  !     LOCATE ALL POINTS IN INTERVAL.
  !
  DO j = jfirst, Ne
    IF ( Xe(j)>=X(ir) ) GOTO 200
  ENDDO
  j = Ne + 1
  GOTO 300
  !
  !     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
  !
  200 CONTINUE
  IF ( ir==N ) j = Ne + 1
  !
  300  nj = j - jfirst
  !
  !     SKIP EVALUATION IF NO POINTS IN INTERVAL.
  !
  IF ( nj/=0 ) THEN
    !
    !     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
    !
    !       ----------------------------------------------------------------
    CALL DCHFDV(X(ir-1),X(ir),F(1,ir-1),F(1,ir),D(1,ir-1),D(1,ir),nj,&
      Xe(jfirst),Fe(jfirst),De(jfirst),next,ierc)
    !       ----------------------------------------------------------------
    IF ( ierc<0 ) GOTO 600
    !
    IF ( next(2)/=0 ) THEN
      !        IF (NEXT(2) .GT. 0)  THEN
      !           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
      !           RIGHT OF X(IR).
      !
      IF ( ir<N ) GOTO 600
      !           IF (IR .EQ. N)  THEN
      !              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
      !           ELSE
      !              WE SHOULD NEVER HAVE GOTTEN HERE.
      Ierr = Ierr + next(2)
    ENDIF
    !           ENDIF
    !        ENDIF
    !
    IF ( next(1)/=0 ) THEN
      !        IF (NEXT(1) .GT. 0)  THEN
      !           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
      !           LEFT OF X(IR-1).
      !
      IF ( ir>2 ) THEN
        !           ELSE
        !              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
        !              EVALUATION INTERVAL.
        !
        !              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
        DO i = jfirst, j - 1
          IF ( Xe(i)<X(ir-1) ) GOTO 320
        ENDDO
        !              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
        !                     IN DCHFDV.
        GOTO 600
      ELSE
        !           IF (IR .EQ. 2)  THEN
        !              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
        Ierr = Ierr + next(1)
        GOTO 350
      ENDIF
      !
      !              RESET J.  (THIS WILL BE THE NEW JFIRST.)
      320      j = i
      !
      !              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
      DO i = 1, ir - 1
        IF ( Xe(j)<X(i) ) EXIT
      ENDDO
      !              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
      !
      !              AT THIS POINT, EITHER  XE(J) .LT. X(1)
      !                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
      !              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
      !              CYCLING.
      ir = MAX(1,i-1)
    ENDIF
    !           ENDIF
    !        ENDIF
    !
    350    jfirst = j
  ENDIF
  !
  !     END OF IR-LOOP.
  !
  ir = ir + 1
  IF ( ir<=N ) GOTO 100
  !
  !  NORMAL RETURN.
  !
  400  RETURN
  !
  !     X-ARRAY NOT STRICTLY INCREASING.
  500  Ierr = -3
  CALL XERMSG('SLATEC','DPCHFD','X-ARRAY NOT STRICTLY INCREASING',Ierr,1)
  RETURN
  !
  !     ERROR RETURN FROM DCHFDV.
  !   *** THIS CASE SHOULD NEVER OCCUR ***
  600  Ierr = -5
  CALL XERMSG('SLATEC','DPCHFD','ERROR RETURN FROM DCHFDV -- FATAL',Ierr,2)
  !------------- LAST LINE OF DPCHFD FOLLOWS -----------------------------
END SUBROUTINE DPCHFD
