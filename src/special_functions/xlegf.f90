!** XLEGF
SUBROUTINE XLEGF(Dnu1,Nudiff,Mu1,Mu2,Theta,Id,Pqa,Ipqa,Ierror)
  !>
  !  Compute normalized Legendre polynomials and associated
  !            Legendre functions.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C3A2, C9
  !***
  ! **Type:**      SINGLE PRECISION (XLEGF-S, DXLEGF-D)
  !***
  ! **Keywords:**  LEGENDRE FUNCTIONS
  !***
  ! **Author:**  Smith, John M., (NBS and George Mason University)
  !***
  ! **Description:**
  !
  !   XLEGF: Extended-range Single-precision Legendre Functions
  !
  !   A feature of the XLEGF subroutine for Legendre functions is
  ! the use of extended-range arithmetic, a software extension of
  ! ordinary floating-point arithmetic that greatly increases the
  ! exponent range of the representable numbers. This avoids the
  ! need for scaling the solutions to lie within the exponent range
  ! of the most restrictive manufacturer's hardware. The increased
  ! exponent range is achieved by allocating an integer storage
  ! location together with each floating-point storage location.
  !
  !   The interpretation of the pair (X,I) where X is floating-point
  ! and I is integer is X*(IR**I) where IR is the internal radix of
  ! the computer arithmetic.
  !
  !   This subroutine computes one of the following vectors:
  !
  ! 1. Legendre function of the first kind of negative order, either
  !    a. P(-MU1,NU,X), P(-MU1-1,NU,X), ..., P(-MU2,NU,X) or
  !    b. P(-MU,NU1,X), P(-MU,NU1+1,X), ..., P(-MU,NU2,X)
  ! 2. Legendre function of the second kind, either
  !    a. Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X) or
  !    b. Q(MU,NU1,X), Q(MU,NU1+1,X), ..., Q(MU,NU2,X)
  ! 3. Legendre function of the first kind of positive order, either
  !    a. P(MU1,NU,X), P(MU1+1,NU,X), ..., P(MU2,NU,X) or
  !    b. P(MU,NU1,X), P(MU,NU1+1,X), ..., P(MU,NU2,X)
  ! 4. Normalized Legendre polynomials, either
  !    a. PN(MU1,NU,X), PN(MU1+1,NU,X), ..., PN(MU2,NU,X) or
  !    b. PN(MU,NU1,X), PN(MU,NU1+1,X), ..., PN(MU,NU2,X)
  !
  ! where X = COS(THETA).
  !
  !   The input values to XLEGF are DNU1, NUDIFF, MU1, MU2, THETA,
  ! and ID. These must satisfy
  !
  !    DNU1 is REAL and greater than or equal to -0.5;
  !    NUDIFF is INTEGER and non-negative;
  !    MU1 is INTEGER and non-negative;
  !    MU2 is INTEGER and greater than or equal to MU1;
  !    THETA is REAL and in the half-open interval (0,PI/2];
  !    ID is INTEGER and equal to 1, 2, 3 or 4;
  !
  ! and  additionally either NUDIFF = 0 or MU2 = MU1.
  !
  !   If ID=1 and NUDIFF=0, a vector of type 1a above is computed
  ! with NU=DNU1.
  !
  !   If ID=1 and MU1=MU2, a vector of type 1b above is computed
  ! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
  !
  !   If ID=2 and NUDIFF=0, a vector of type 2a above is computed
  ! with NU=DNU1.
  !
  !   If ID=2 and MU1=MU2, a vector of type 2b above is computed
  ! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
  !
  !   If ID=3 and NUDIFF=0, a vector of type 3a above is computed
  ! with NU=DNU1.
  !
  !   If ID=3 and MU1=MU2, a vector of type 3b above is computed
  ! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
  !
  !   If ID=4 and NUDIFF=0, a vector of type 4a above is computed
  ! with NU=DNU1.
  !
  !   If ID=4 and MU1=MU2, a vector of type 4b above is computed
  ! with NU1=DNU1, NU2=DNU1+NUDIFF and MU=MU1.
  !
  !   In each case the vector of computed Legendre function values
  ! is returned in the extended-range vector (PQA(I),IPQA(I)). The
  ! length of this vector is either MU2-MU1+1 or NUDIFF+1.
  !
  !   Where possible, XLEGF returns IPQA(I) as zero. In this case the
  ! value of the Legendre function is contained entirely in PQA(I),
  ! so it can be used in subsequent computations without further
  ! consideration of extended-range arithmetic. If IPQA(I) is nonzero,
  ! then the value of the Legendre function is not representable in
  ! floating-point because of underflow or overflow. The program that
  ! calls XLEGF must test IPQA(I) to ensure correct usage.
  !
  !   IERROR is an error indicator. If no errors are detected, IERROR=0
  ! when control returns to the calling routine. If an error is detected,
  ! IERROR is returned as nonzero. The calling routine must check the
  ! value of IERROR.
  !
  !   If IERROR=110 or 111, invalid input was provided to XLEGF.
  !   If IERROR=101,102,103, or 104, invalid input was provided to XSET.
  !   If IERROR=105 or 106, an internal consistency error occurred in
  ! XSET (probably due to a software malfunction in the library routine
  ! I1MACH).
  !   If IERROR=107, an overflow or underflow of an extended-range number
  ! was detected in XADJ.
  !   If IERROR=108, an overflow or underflow of an extended-range number
  ! was detected in XC210.
  !
  !***
  ! **See also:**  XSET
  !***
  ! **References:**  Olver and Smith, Associated Legendre Functions on the
  !                 Cut, J Comp Phys, v 51, n 3, Sept 1983, pp 502--518.
  !               Smith, Olver and Lozier, Extended-Range Arithmetic and
  !                 Normalized Legendre Polynomials, ACM Trans on Math
  !                 Softw, v 7, n 1, March 1981, pp 93--105.
  !***
  ! **Routines called:**  XERMSG, XPMU, XPMUP, XPNRM, XPQNU, XQMU, XQNU,
  !                    XRED, XSET

  !* REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)
  USE service, ONLY : XERMSG
  INTEGER :: Id, Ierror, Mu1, Mu2, Nudiff, Ipqa(Nudiff+Mu2-Mu1+1)
  REAL :: Pqa(Nudiff+Mu2-Mu1+1), Dnu1, Theta
  INTEGER :: i, l
  REAL :: dnu2, sx, x, pi2
  !
  !* FIRST EXECUTABLE STATEMENT  XLEGF
  Ierror = 0
  CALL XSET(0,0,0.0,0,Ierror)
  IF ( Ierror/=0 ) RETURN
  pi2 = 2.*ATAN(1.)
  !
  !        ZERO OUTPUT ARRAYS
  !
  l = (Mu2-Mu1) + Nudiff + 1
  DO i = 1, l
    Pqa(i) = 0.
    Ipqa(i) = 0
  END DO
  !
  !        CHECK FOR VALID INPUT VALUES
  !
  IF ( Nudiff>=0 ) THEN
    IF ( Dnu1>=-.5 ) THEN
      IF ( Mu2>=Mu1 ) THEN
        IF ( Mu1>=0 ) THEN
          IF ( Theta<=0..OR.Theta>pi2 ) THEN
            CALL XERMSG('SLATEC','XLEGF','THETA out of range',111,1)
            Ierror = 111
            RETURN
          ELSEIF ( Id>=1.AND.Id<=4 ) THEN
            IF ( (Mu1==Mu2).OR.(Nudiff<=0) ) THEN
              !
              !        IF DNU1 IS NOT AN INTEGER, NORMALIZED P(MU,DNU,X)
              !        CANNOT BE CALCULATED.  IF DNU1 IS AN INTEGER AND
              !        MU1.GT.DNU2 THEN ALL VALUES OF P(+MU,DNU,X) AND
              !        NORMALIZED P(MU,NU,X) WILL BE ZERO.
              !
              dnu2 = Dnu1 + Nudiff
              IF ( (Id/=3).OR.(MOD(Dnu1,1.)==0.) ) THEN
                IF ( (Id==4).AND.(MOD(Dnu1,1.)/=0.) ) GOTO 100
                IF ( (Id==3.OR.Id==4).AND.Mu1>dnu2 ) RETURN
              END IF
              !
              x = COS(Theta)
              sx = 1./SIN(Theta)
              IF ( Id/=2 ) THEN
                IF ( Mu2<=Mu1 ) THEN
                  !
                  !        FIXED MU, VARIABLE NU
                  !        CALL XPQNU TO CALCULATE P(-MU,DNU1,X),....,P(-MU,DNU2,X)
                  !
                  CALL XPQNU(Dnu1,dnu2,Mu1,Theta,Id,Pqa,Ipqa,Ierror)
                  IF ( Ierror/=0 ) RETURN
                ELSE
                  !
                  !        FIXED NU, VARIABLE MU
                  !        CALL XPMU TO CALCULATE P(-MU1,NU,X),....,P(-MU2,NU,X)
                  !
                  CALL XPMU(Dnu1,dnu2,Mu1,Mu2,Theta,x,sx,Id,Pqa,Ipqa,Ierror)
                  IF ( Ierror/=0 ) RETURN
                END IF
                !
                !        IF ID = 3, TRANSFORM P(-MU,NU,X) VECTOR INTO
                !        P(MU,NU,X) VECTOR.
                !
                IF ( Id==3 ) CALL XPMUP(Dnu1,dnu2,Mu1,Mu2,Pqa,Ipqa,Ierror)
                IF ( Ierror/=0 ) RETURN
                !
                !        IF ID = 4, TRANSFORM P(-MU,NU,X) VECTOR INTO
                !        NORMALIZED P(MU,NU,X) VECTOR.
                !
                IF ( Id==4 ) CALL XPNRM(Dnu1,dnu2,Mu1,Mu2,Pqa,Ipqa,Ierror)
                IF ( Ierror/=0 ) RETURN
                !
              ELSEIF ( Mu2==Mu1 ) THEN
                !
                !        FIXED MU, VARIABLE NU
                !        CALL XQNU TO CALCULATE Q(MU,DNU1,X),....,Q(MU,DNU2,X)
                !
                CALL XQNU(Dnu1,dnu2,Mu1,Theta,x,sx,Id,Pqa,Ipqa,Ierror)
                IF ( Ierror/=0 ) RETURN
              ELSE
                !
                !        FIXED NU, VARIABLE MU
                !        CALL XQMU TO CALCULATE Q(MU1,NU,X),....,Q(MU2,NU,X)
                !
                CALL XQMU(Dnu1,dnu2,Mu1,Mu2,Theta,x,sx,Id,Pqa,Ipqa,Ierror)
                IF ( Ierror/=0 ) RETURN
              END IF
              !
              !        PLACE RESULTS IN REDUCED FORM IF POSSIBLE
              !        AND RETURN TO MAIN PROGRAM.
              !
              DO i = 1, l
                CALL XRED(Pqa(i),Ipqa(i),Ierror)
                IF ( Ierror/=0 ) RETURN
              END DO
              RETURN
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
  !
  !        *****     ERROR TERMINATION     *****
  !
  100  CALL XERMSG('SLATEC','XLEGF','DNU1, NUDIFF, MU1, MU2, or ID not valid',&
    110,1)
  Ierror = 110
  RETURN
END SUBROUTINE XLEGF
