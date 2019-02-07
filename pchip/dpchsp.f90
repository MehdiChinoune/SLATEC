!*==DPCHSP.f90  processed by SPAG 6.72Dc at 11:00 on  6 Feb 2019
!DECK DPCHSP
SUBROUTINE DPCHSP(Ic,Vc,N,X,F,D,Incfd,Wk,Nwk,Ierr)
  IMPLICIT NONE
  !*--DPCHSP5
  !***BEGIN PROLOGUE  DPCHSP
  !***PURPOSE  Set derivatives needed to determine the Hermite represen-
  !            tation of the cubic spline interpolant to given data, with
  !            specified boundary conditions.
  !***LIBRARY   SLATEC (PCHIP)
  !***CATEGORY  E1A
  !***TYPE      DOUBLE PRECISION (PCHSP-S, DPCHSP-D)
  !***KEYWORDS  CUBIC HERMITE INTERPOLATION, PCHIP,
  !             PIECEWISE CUBIC INTERPOLATION, SPLINE INTERPOLATION
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             P.O. Box 808  (L-316)
  !             Livermore, CA  94550
  !             FTS 532-4275, (510) 422-4275
  !***DESCRIPTION
  !
  !          DPCHSP:   Piecewise Cubic Hermite Spline
  !
  !     Computes the Hermite representation of the cubic spline inter-
  !     polant to the data given in X and F satisfying the boundary
  !     conditions specified by IC and VC.
  !
  !     To facilitate two-dimensional applications, includes an increment
  !     between successive values of the F- and D-arrays.
  !
  !     The resulting piecewise cubic Hermite function may be evaluated
  !     by DPCHFE or DPCHFD.
  !
  !     NOTE:  This is a modified version of C. de Boor's cubic spline
  !            routine CUBSPL.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        PARAMETER  (INCFD = ...)
  !        INTEGER  IC(2), N, NWK, IERR
  !        DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
  !
  !        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
  !
  !   Parameters:
  !
  !     IC -- (input) integer array of length 2 specifying desired
  !           boundary conditions:
  !           IC(1) = IBEG, desired condition at beginning of data.
  !           IC(2) = IEND, desired condition at end of data.
  !
  !           IBEG = 0  to set D(1) so that the third derivative is con-
  !              tinuous at X(2).  This is the "not a knot" condition
  !              provided by de Boor's cubic spline routine CUBSPL.
  !              < This is the default boundary condition. >
  !           IBEG = 1  if first derivative at X(1) is given in VC(1).
  !           IBEG = 2  if second derivative at X(1) is given in VC(1).
  !           IBEG = 3  to use the 3-point difference formula for D(1).
  !                     (Reverts to the default b.c. if N.LT.3 .)
  !           IBEG = 4  to use the 4-point difference formula for D(1).
  !                     (Reverts to the default b.c. if N.LT.4 .)
  !          NOTES:
  !           1. An error return is taken if IBEG is out of range.
  !           2. For the "natural" boundary condition, use IBEG=2 and
  !              VC(1)=0.
  !
  !           IEND may take on the same values as IBEG, but applied to
  !           derivative at X(N).  In case IEND = 1 or 2, the value is
  !           given in VC(2).
  !
  !          NOTES:
  !           1. An error return is taken if IEND is out of range.
  !           2. For the "natural" boundary condition, use IEND=2 and
  !              VC(2)=0.
  !
  !     VC -- (input) real*8 array of length 2 specifying desired boundary
  !           values, as indicated above.
  !           VC(1) need be set only if IC(1) = 1 or 2 .
  !           VC(2) need be set only if IC(2) = 1 or 2 .
  !
  !     N -- (input) number of data points.  (Error return if N.LT.2 .)
  !
  !     X -- (input) real*8 array of independent variable values.  The
  !           elements of X must be strictly increasing:
  !                X(I-1) .LT. X(I),  I = 2(1)N.
  !           (Error return if not.)
  !
  !     F -- (input) real*8 array of dependent variable values to be
  !           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
  !           X(I).
  !
  !     D -- (output) real*8 array of derivative values at the data
  !           points.  These values will determine the cubic spline
  !           interpolant with the requested boundary conditions.
  !           The value corresponding to X(I) is stored in
  !                D(1+(I-1)*INCFD),  I=1(1)N.
  !           No other entries in D are changed.
  !
  !     INCFD -- (input) increment between successive values in F and D.
  !           This argument is provided primarily for 2-D applications.
  !           (Error return if  INCFD.LT.1 .)
  !
  !     WK -- (scratch) real*8 array of working storage.
  !
  !     NWK -- (input) length of work array.
  !           (Error return if NWK.LT.2*N .)
  !
  !     IERR -- (output) error flag.
  !           Normal return:
  !              IERR = 0  (no errors).
  !           "Recoverable" errors:
  !              IERR = -1  if N.LT.2 .
  !              IERR = -2  if INCFD.LT.1 .
  !              IERR = -3  if the X-array is not strictly increasing.
  !              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
  !              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
  !              IERR = -6  if both of the above are true.
  !              IERR = -7  if NWK is too small.
  !               NOTE:  The above errors are checked in the order listed,
  !                   and following arguments have **NOT** been validated.
  !             (The D-array has not been changed in any of these cases.)
  !              IERR = -8  in case of trouble solving the linear system
  !                         for the interior derivative values.
  !             (The D-array may have been changed in this case.)
  !             (             Do **NOT** use it!                )
  !
  !***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
  !                 Verlag, New York, 1978, pp. 53-59.
  !***ROUTINES CALLED  DPCHDF, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   820503  DATE WRITTEN
  !   820804  Converted to SLATEC library version.
  !   870707  Corrected XERROR calls for d.p. name(s).
  !   890206  Corrected XERROR calls.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890703  Corrected category record.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920429  Revised format and order of references.  (WRB,FNF)
  !***END PROLOGUE  DPCHSP
  !  Programming notes:
  !
  !     To produce a single precision version, simply:
  !        a. Change DPCHSP to PCHSP wherever it occurs,
  !        b. Change the double precision declarations to real, and
  !        c. Change the constants ZERO, HALF, ... to single precision.
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER Ic(2) , N , Incfd , Nwk , Ierr
  REAL(8) :: Vc(2) , X(*) , F(Incfd,*) , D(Incfd,*) , Wk(2,*)
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER ibeg , iend , index , j , nm1
  REAL(8) :: g , half , one , stemp(3) , three , two , xtemp(4) , zero
  SAVE zero , half , one , two , three
  REAL(8) :: DPCHDF
  !
  DATA zero/0.D0/ , half/.5D0/ , one/1.D0/ , two/2.D0/ , three/3.D0/
  !
  !  VALIDITY-CHECK ARGUMENTS.
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHSP
  IF ( N<2 ) THEN
    !
    !  ERROR RETURNS.
    !
    !     N.LT.2 RETURN.
    Ierr = -1
    CALL XERMSG('SLATEC','DPCHSP','NUMBER OF DATA POINTS LESS THAN TWO',&
      Ierr,1)
    RETURN
  ELSE
    IF ( Incfd<1 ) THEN
      !
      !     INCFD.LT.1 RETURN.
      Ierr = -2
      CALL XERMSG('SLATEC','DPCHSP','INCREMENT LESS THAN ONE',Ierr,1)
      RETURN
    ELSE
      DO j = 2 , N
        IF ( X(j)<=X(j-1) ) GOTO 20
      ENDDO
      !
      ibeg = Ic(1)
      iend = Ic(2)
      Ierr = 0
      IF ( (ibeg<0).OR.(ibeg>4) ) Ierr = Ierr - 1
      IF ( (iend<0).OR.(iend>4) ) Ierr = Ierr - 2
      IF ( Ierr<0 ) THEN
        !
        !     IC OUT OF RANGE RETURN.
        Ierr = Ierr - 3
        CALL XERMSG('SLATEC','DPCHSP','IC OUT OF RANGE',Ierr,1)
        RETURN
        !
        !  FUNCTION DEFINITION IS OK -- GO ON.
        !
      ELSEIF ( Nwk<2*N ) THEN
        !
        !     NWK TOO SMALL RETURN.
        Ierr = -7
        CALL XERMSG('SLATEC','DPCHSP','WORK ARRAY TOO SMALL',Ierr,1)
        RETURN
      ELSE
        !
        !  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
        !  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
        DO j = 2 , N
          Wk(1,j) = X(j) - X(j-1)
          Wk(2,j) = (F(1,j)-F(1,j-1))/Wk(1,j)
        ENDDO
        !
        !  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
        !
        IF ( ibeg>N ) ibeg = 0
        IF ( iend>N ) iend = 0
        !
        !  SET UP FOR BOUNDARY CONDITIONS.
        !
        IF ( (ibeg==1).OR.(ibeg==2) ) THEN
          D(1,1) = Vc(1)
        ELSEIF ( ibeg>2 ) THEN
          !        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
          DO j = 1 , ibeg
            index = ibeg - j + 1
            !           INDEX RUNS FROM IBEG DOWN TO 1.
            xtemp(j) = X(index)
            IF ( j<ibeg ) stemp(j) = Wk(2,index)
          ENDDO
          !                 --------------------------------
          D(1,1) = DPCHDF(ibeg,xtemp,stemp,Ierr)
          !                 --------------------------------
          IF ( Ierr/=0 ) GOTO 100
          ibeg = 1
        ENDIF
        !
        IF ( (iend==1).OR.(iend==2) ) THEN
          D(1,N) = Vc(2)
        ELSEIF ( iend>2 ) THEN
          !        PICK UP LAST IEND POINTS.
          DO j = 1 , iend
            index = N - iend + j
            !           INDEX RUNS FROM N+1-IEND UP TO N.
            xtemp(j) = X(index)
            IF ( j<iend ) stemp(j) = Wk(2,index+1)
          ENDDO
          !                 --------------------------------
          D(1,N) = DPCHDF(iend,xtemp,stemp,Ierr)
          !                 --------------------------------
          IF ( Ierr/=0 ) GOTO 100
          iend = 1
        ENDIF
        !
        ! --------------------( BEGIN CODING FROM CUBSPL )--------------------
        !
        !  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
        !  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
        !  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
        !     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
        !
        !  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
        !             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
        !
        IF ( ibeg==0 ) THEN
          IF ( N==2 ) THEN
            !           NO CONDITION AT LEFT END AND N = 2.
            Wk(2,1) = one
            Wk(1,1) = one
            D(1,1) = two*Wk(2,2)
          ELSE
            !           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            Wk(2,1) = Wk(1,3)
            Wk(1,1) = Wk(1,2) + Wk(1,3)
            D(1,1) = ((Wk(1,2)+two*Wk(1,1))*Wk(2,2)*Wk(1,3)+Wk(1,2)&
              **2*Wk(2,3))/Wk(1,1)
          ENDIF
        ELSEIF ( ibeg==1 ) THEN
          !        SLOPE PRESCRIBED AT LEFT END.
          Wk(2,1) = one
          Wk(1,1) = zero
        ELSE
          !        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
          Wk(2,1) = two
          Wk(1,1) = one
          D(1,1) = three*Wk(2,2) - half*Wk(1,2)*D(1,1)
        ENDIF
        !
        !  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
        !  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
        !  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
        !
        nm1 = N - 1
        IF ( nm1>1 ) THEN
          DO j = 2 , nm1
            IF ( Wk(2,j-1)==zero ) GOTO 50
            g = -Wk(1,j+1)/Wk(2,j-1)
            D(1,j) = g*D(1,j-1) + three*(Wk(1,j)*Wk(2,j+1)+Wk(1,j+1)*Wk(2,j)&
              )
            Wk(2,j) = g*Wk(1,j-1) + two*(Wk(1,j)+Wk(1,j+1))
          ENDDO
        ENDIF
        !
        !  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
        !           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
        !
        !     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
        !     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
        !     AT THIS POINT.
        IF ( iend/=1 ) THEN
          !
          IF ( iend/=0 ) THEN
            !        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
            D(1,N) = three*Wk(2,N) + half*Wk(1,N)*D(1,N)
            Wk(2,N) = two
            IF ( Wk(2,N-1)==zero ) GOTO 50
            g = -one/Wk(2,N-1)
          ELSEIF ( N==2.AND.ibeg==0 ) THEN
            !           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = Wk(2,2)
            GOTO 10
          ELSEIF ( (N==2).OR.(N==3.AND.ibeg==0) ) THEN
            !           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
            !           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = two*Wk(2,N)
            Wk(2,N) = one
            IF ( Wk(2,N-1)==zero ) GOTO 50
            g = -one/Wk(2,N-1)
          ELSE
            !           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
            !           KNOT AT LEFT END POINT.
            g = Wk(1,N-1) + Wk(1,N)
            !           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((Wk(1,N)+two*g)*Wk(2,N)*Wk(1,N-1)+Wk(1,N)&
              **2*(F(1,N-1)-F(1,N-2))/Wk(1,N-1))/g
            IF ( Wk(2,N-1)==zero ) GOTO 50
            g = -g/Wk(2,N-1)
            Wk(2,N) = Wk(1,N-1)
          ENDIF
          !
          !  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
          !
          Wk(2,N) = g*Wk(1,N-1) + Wk(2,N)
          IF ( Wk(2,N)==zero ) GOTO 50
          D(1,N) = (g*D(1,N-1)+D(1,N))/Wk(2,N)
        ENDIF
        !
        !  CARRY OUT BACK SUBSTITUTION
        !
        10         DO j = nm1 , 1 , -1
        IF ( Wk(2,j)==zero ) GOTO 50
        D(1,j) = (D(1,j)-Wk(1,j)*D(1,j+1))/Wk(2,j)
      ENDDO
      ! --------------------(  END  CODING FROM CUBSPL )--------------------
      !
      !  NORMAL RETURN.
      !
      RETURN
      ENDIF
      !
      !     X-ARRAY NOT STRICTLY INCREASING.
      20       Ierr = -3
      CALL XERMSG('SLATEC','DPCHSP','X-ARRAY NOT STRICTLY INCREASING',Ierr,&
        1)
      RETURN
    ENDIF
    !
    !     SINGULAR SYSTEM.
    !   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
    !   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
    50     Ierr = -8
    CALL XERMSG('SLATEC','DPCHSP','SINGULAR LINEAR SYSTEM',Ierr,1)
    RETURN
  ENDIF
  !
  !     ERROR RETURN FROM DPCHDF.
  !   *** THIS CASE SHOULD NEVER OCCUR ***
  100  Ierr = -9
  CALL XERMSG('SLATEC','DPCHSP','ERROR RETURN FROM DPCHDF',Ierr,1)
  !------------- LAST LINE OF DPCHSP FOLLOWS -----------------------------
END SUBROUTINE DPCHSP
