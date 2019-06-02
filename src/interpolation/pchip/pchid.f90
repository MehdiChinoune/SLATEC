!** PCHID
REAL FUNCTION PCHID(N,X,F,D,Incfd,Skip,Ia,Ib,Ierr)
  !>
  !  Evaluate the definite integral of a piecewise cubic
  !            Hermite function over an interval whose endpoints are data
  !            points.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Category:**  E3, H2A1B2
  !***
  ! **Type:**      SINGLE PRECISION (PCHID-S, DPCHID-D)
  !***
  ! **Keywords:**  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
  !             QUADRATURE
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***
  ! **Description:**
  !
  !          PCHID:  Piecewise Cubic Hermite Integrator, Data limits
  !
  !     Evaluates the definite integral of the cubic Hermite function
  !     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
  !
  !     To provide compatibility with PCHIM and PCHIC, includes an
  !     increment between successive values of the F- and D-arrays.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  N, IA, IB, IERR
  !        REAL  X(N), F(INCFD,N), D(INCFD,N)
  !        LOGICAL  SKIP
  !
  !        VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
  !
  !   Parameters:
  !
  !     VALUE -- (output) value of the requested integral.
  !
  !     N -- (input) number of data points.  (Error return if N.LT.2 .)
  !
  !     X -- (input) real array of independent variable values.  The
  !           elements of X must be strictly increasing:
  !                X(I-1) .LT. X(I),  I = 2(1)N.
  !           (Error return if not.)
  !
  !     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
  !           the value corresponding to X(I).
  !
  !     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
  !           the value corresponding to X(I).
  !
  !     INCFD -- (input) increment between successive values in F and D.
  !           (Error return if  INCFD.LT.1 .)
  !
  !     SKIP -- (input/output) logical variable which should be set to
  !           .TRUE. if the user wishes to skip checks for validity of
  !           preceding parameters, or to .FALSE. otherwise.
  !           This will save time in case these checks have already
  !           been performed (say, in PCHIM or PCHIC).
  !           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
  !
  !     IA,IB -- (input) indices in X-array for the limits of integration.
  !           both must be in the range [1,N].  (Error return if not.)
  !           No restrictions on their relative values.
  !
  !     IERR -- (output) error flag.
  !           Normal return:
  !              IERR = 0  (no errors).
  !           "Recoverable" errors:
  !              IERR = -1  if N.LT.2 .
  !              IERR = -2  if INCFD.LT.1 .
  !              IERR = -3  if the X-array is not strictly increasing.
  !              IERR = -4  if IA or IB is out of range.
  !                (VALUE will be zero in any of these cases.)
  !               NOTE:  The above errors are checked in the order listed,
  !                   and following arguments have **NOT** been validated.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820723  DATE WRITTEN
  !   820804  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890703  Corrected category record.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   930504  Corrected to set VALUE=0 when IERR.ne.0.  (FNF)
  USE service, ONLY : XERMSG
  !
  !  Programming notes:
  !  1. This routine uses a special formula that is valid only for
  !     integrals whose limits coincide with data values.  This is
  !     mathematically equivalent to, but much more efficient than,
  !     calls to CHFIE.
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER N, Incfd, Ia, Ib, Ierr
  REAL X(N), F(Incfd,N), D(Incfd,N)
  LOGICAL Skip
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER i, iup, low
  REAL h, summ, value
  !
  !  INITIALIZE.
  !
  REAL, PARAMETER :: zero = 0., half = 0.5, six = 6.
  !* FIRST EXECUTABLE STATEMENT  PCHID
  value = zero
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  IF ( .NOT.(Skip) ) THEN
    !
    IF ( N<2 ) THEN
      !
      !  ERROR RETURNS.
      !
      !     N.LT.2 RETURN.
      Ierr = -1
      CALL XERMSG('PCHID','NUMBER OF DATA POINTS LESS THAN TWO',Ierr,1)
      GOTO 100
    ELSEIF ( Incfd<1 ) THEN
      !
      !     INCFD.LT.1 RETURN.
      Ierr = -2
      CALL XERMSG('PCHID','INCREMENT LESS THAN ONE',Ierr,1)
      GOTO 100
    ELSE
      DO i = 2, N
        IF ( X(i)<=X(i-1) ) GOTO 200
      END DO
    END IF
  END IF
  !
  !  FUNCTION DEFINITION IS OK, GO ON.
  !
  Skip = .TRUE.
  IF ( (Ia<1).OR.(Ia>N) ) GOTO 300
  IF ( (Ib<1).OR.(Ib>N) ) GOTO 300
  Ierr = 0
  !
  !  COMPUTE INTEGRAL VALUE.
  !
  IF ( Ia/=Ib ) THEN
    low = MIN(Ia,Ib)
    iup = MAX(Ia,Ib) - 1
    summ = zero
    DO i = low, iup
      h = X(i+1) - X(i)
      summ = summ + h*((F(1,i)+F(1,i+1))+(D(1,i)-D(1,i+1))*(h/six))
    END DO
    value = half*summ
    IF ( Ia>Ib ) value = -value
  END IF
  !
  !  NORMAL RETURN.
  !
  100  PCHID = value
  RETURN
  !
  !     X-ARRAY NOT STRICTLY INCREASING.
  200  Ierr = -3
  CALL XERMSG('PCHID','X-ARRAY NOT STRICTLY INCREASING',Ierr,1)
  GOTO 100
  !
  !     IA OR IB OUT OF RANGE RETURN.
  300  Ierr = -4
  CALL XERMSG('PCHID','IA OR IB OUT OF RANGE',Ierr,1)
  GOTO 100
  !------------- LAST LINE OF PCHID FOLLOWS ------------------------------
END FUNCTION PCHID
