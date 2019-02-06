!*==PCHIA.f90  processed by SPAG 6.72Dc at 11:00 on  6 Feb 2019
!DECK PCHIA
      REAL FUNCTION PCHIA(N,X,F,D,Incfd,Skip,A,B,Ierr)
      IMPLICIT NONE
!*--PCHIA5
!***BEGIN PROLOGUE  PCHIA
!***PURPOSE  Evaluate the definite integral of a piecewise cubic
!            Hermite function over an arbitrary interval.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H2A1B2
!***TYPE      SINGLE PRECISION (PCHIA-S, DPCHIA-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
!             QUADRATURE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [A, B].
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
!        REAL  VALUE, PCHIA
!        LOGICAL  SKIP
!
!        VALUE = PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
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
!           SKIP will be set to .TRUE. on return with IERR.GE.0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -4  in case of an error return from PCHID (which
!                         should never occur).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CHFIE, PCHID, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820730  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected double precision conversion instructions.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   930503  Corrected to set VALUE=0 when IERR.lt.0.  (FNF)
!   930504  Changed CHFIV to CHFIE.  (FNF)
!***END PROLOGUE  PCHIA
!
!  Programming notes:
!  1. The error flag from PCHID is tested, because a logic flaw
!     could conceivably result in IERD=-4, which should be reported.
!**End
!
!  DECLARE ARGUMENTS.
!
      INTEGER N , Incfd , Ierr
      REAL X(*) , F(Incfd,*) , D(Incfd,*) , A , B
      LOGICAL Skip
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER i , ia , ib , ierd , il , ir
      REAL value , xa , xb , zero
      SAVE zero
      REAL CHFIE , PCHID
!
!  INITIALIZE.
!
      DATA zero/0./
!***FIRST EXECUTABLE STATEMENT  PCHIA
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
          CALL XERMSG('SLATEC','PCHIA','NUMBER OF DATA POINTS LESS THAN TWO',
     &                Ierr,1)
          GOTO 100
        ELSEIF ( Incfd<1 ) THEN
!
!     INCFD.LT.1 RETURN.
          Ierr = -2
          CALL XERMSG('SLATEC','PCHIA','INCREMENT LESS THAN ONE',Ierr,1)
          GOTO 100
        ELSE
          DO i = 2 , N
            IF ( X(i)<=X(i-1) ) GOTO 200
          ENDDO
        ENDIF
      ENDIF
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
      Skip = .TRUE.
      Ierr = 0
      IF ( (A<X(1)).OR.(A>X(N)) ) Ierr = Ierr + 1
      IF ( (B<X(1)).OR.(B>X(N)) ) Ierr = Ierr + 2
!
!  COMPUTE INTEGRAL VALUE.
!
      IF ( A/=B ) THEN
        xa = MIN(A,B)
        xb = MAX(A,B)
        IF ( xb<=X(2) ) THEN
!           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!                   --------------------------------------
          value = CHFIE(X(1),X(2),F(1,1),F(1,2),D(1,1),D(1,2),A,B)
!                   --------------------------------------
        ELSEIF ( xa>=X(N-1) ) THEN
!           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!                   -----------------------------------------
          value = CHFIE(X(N-1),X(N),F(1,N-1),F(1,N),D(1,N-1),D(1,N),A,B)
!                   -----------------------------------------
        ELSE
!           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
!      ......LOCATE IA AND IB SUCH THAT
!               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
          ia = 1
          DO i = 1 , N - 1
            IF ( xa>X(i) ) ia = i + 1
          ENDDO
!             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
!             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
!
          ib = N
          DO i = N , ia , -1
            IF ( xb<X(i) ) ib = i - 1
          ENDDO
!             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
!             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
!
!     ......COMPUTE THE INTEGRAL.
          IF ( ib<ia ) THEN
!              THIS MEANS IB = IA-1 AND
!                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
!                      ------------------------------------------
            value = CHFIE(X(ib),X(ia),F(1,ib),F(1,ia),D(1,ib),D(1,ia),A,B)
!                      ------------------------------------------
          ELSE
!
!              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
!                (Case (IB .EQ. IA) is taken care of by initialization
!                 of VALUE to ZERO.)
            IF ( ib>ia ) THEN
!                         ---------------------------------------------
              value = PCHID(N,X,F,D,Incfd,Skip,ia,ib,ierd)
!                         ---------------------------------------------
              IF ( ierd<0 ) THEN
!
!     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
                Ierr = -4
                CALL XERMSG('SLATEC','PCHIA','TROUBLE IN PCHID',Ierr,1)
                GOTO 100
              ENDIF
            ENDIF
!
!              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
            IF ( xa<X(ia) ) THEN
              il = MAX(1,ia-1)
              ir = il + 1
!                                 -------------------------------------
              value = value + CHFIE(X(il),X(ir),F(1,il),F(1,ir),D(1,il),D(1,ir),
     &                xa,X(ia))
!                                 -------------------------------------
            ENDIF
!
!              THEN ADD ON INTEGRAL OVER (X(IB),XB).
            IF ( xb>X(ib) ) THEN
              ir = MIN(ib+1,N)
              il = ir - 1
!                                 -------------------------------------
              value = value + CHFIE(X(il),X(ir),F(1,il),F(1,ir),D(1,il),D(1,ir),
     &                X(ib),xb)
!                                 -------------------------------------
            ENDIF
!
!              FINALLY, ADJUST SIGN IF NECESSARY.
            IF ( A>B ) value = -value
          ENDIF
        ENDIF
      ENDIF
!
!  NORMAL RETURN.
!
 100  PCHIA = value
      RETURN
!
!     X-ARRAY NOT STRICTLY INCREASING.
 200  Ierr = -3
      CALL XERMSG('SLATEC','PCHIA','X-ARRAY NOT STRICTLY INCREASING',Ierr,1)
      GOTO 100
!------------- LAST LINE OF PCHIA FOLLOWS ------------------------------
      END FUNCTION PCHIA
