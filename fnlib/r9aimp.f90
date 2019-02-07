!*==R9AIMP.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK R9AIMP
SUBROUTINE R9AIMP(X,Ampl,Theta)
  IMPLICIT NONE
  !*--R9AIMP5
  !*** Start of declarations inserted by SPAG
  REAL am21cs, am22cs, Ampl, ath1cs, ath2cs, CSEVL, eta, pi4, &
    R1MACH, sqrtx, Theta, X, xsml, z
  INTEGER INITS, nam21, nam22, nath1, nath2
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  R9AIMP
  !***SUBSIDIARY
  !***PURPOSE  Evaluate the Airy modulus and phase.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      SINGLE PRECISION (R9AIMP-S, D9AIMP-D)
  !***KEYWORDS  AIRY FUNCTION, FNLIB, MODULUS, PHASE, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate the Airy modulus and phase for X .LE. -1.0
  !
  ! Series for AM21       on the interval -1.25000D-01 to  0.
  !                                        with weighted error   2.89E-17
  !                                         log weighted error  16.54
  !                               significant figures required  14.15
  !                                    decimal places required  17.34
  !
  ! Series for ATH1       on the interval -1.25000D-01 to  0.
  !                                        with weighted error   2.53E-17
  !                                         log weighted error  16.60
  !                               significant figures required  15.15
  !                                    decimal places required  17.38
  !
  ! Series for AM22       on the interval -1.00000D+00 to -1.25000D-01
  !                                        with weighted error   2.99E-17
  !                                         log weighted error  16.52
  !                               significant figures required  14.57
  !                                    decimal places required  17.28
  !
  ! Series for ATH2       on the interval -1.00000D+00 to -1.25000D-01
  !                                        with weighted error   2.57E-17
  !                                         log weighted error  16.59
  !                               significant figures required  15.07
  !                                    decimal places required  17.34
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  R9AIMP
  DIMENSION am21cs(40), ath1cs(36), am22cs(33), ath2cs(32)
  LOGICAL first
  SAVE am21cs, ath1cs, am22cs, ath2cs, pi4, nam21, nath1, nam22, &
    nath2, xsml, first
  DATA am21cs(1)/.0065809191761485E0/
  DATA am21cs(2)/.0023675984685722E0/
  DATA am21cs(3)/.0001324741670371E0/
  DATA am21cs(4)/.0000157600904043E0/
  DATA am21cs(5)/.0000027529702663E0/
  DATA am21cs(6)/.0000006102679017E0/
  DATA am21cs(7)/.0000001595088468E0/
  DATA am21cs(8)/.0000000471033947E0/
  DATA am21cs(9)/.0000000152933871E0/
  DATA am21cs(10)/.0000000053590722E0/
  DATA am21cs(11)/.0000000020000910E0/
  DATA am21cs(12)/.0000000007872292E0/
  DATA am21cs(13)/.0000000003243103E0/
  DATA am21cs(14)/.0000000001390106E0/
  DATA am21cs(15)/.0000000000617011E0/
  DATA am21cs(16)/.0000000000282491E0/
  DATA am21cs(17)/.0000000000132979E0/
  DATA am21cs(18)/.0000000000064188E0/
  DATA am21cs(19)/.0000000000031697E0/
  DATA am21cs(20)/.0000000000015981E0/
  DATA am21cs(21)/.0000000000008213E0/
  DATA am21cs(22)/.0000000000004296E0/
  DATA am21cs(23)/.0000000000002284E0/
  DATA am21cs(24)/.0000000000001232E0/
  DATA am21cs(25)/.0000000000000675E0/
  DATA am21cs(26)/.0000000000000374E0/
  DATA am21cs(27)/.0000000000000210E0/
  DATA am21cs(28)/.0000000000000119E0/
  DATA am21cs(29)/.0000000000000068E0/
  DATA am21cs(30)/.0000000000000039E0/
  DATA am21cs(31)/.0000000000000023E0/
  DATA am21cs(32)/.0000000000000013E0/
  DATA am21cs(33)/.0000000000000008E0/
  DATA am21cs(34)/.0000000000000005E0/
  DATA am21cs(35)/.0000000000000003E0/
  DATA am21cs(36)/.0000000000000001E0/
  DATA am21cs(37)/.0000000000000001E0/
  DATA am21cs(38)/.0000000000000000E0/
  DATA am21cs(39)/.0000000000000000E0/
  DATA am21cs(40)/.0000000000000000E0/
  DATA ath1cs(1)/ - .07125837815669365E0/
  DATA ath1cs(2)/ - .00590471979831451E0/
  DATA ath1cs(3)/ - .00012114544069499E0/
  DATA ath1cs(4)/ - .00000988608542270E0/
  DATA ath1cs(5)/ - .00000138084097352E0/
  DATA ath1cs(6)/ - .00000026142640172E0/
  DATA ath1cs(7)/ - .00000006050432589E0/
  DATA ath1cs(8)/ - .00000001618436223E0/
  DATA ath1cs(9)/ - .00000000483464911E0/
  DATA ath1cs(10)/ - .00000000157655272E0/
  DATA ath1cs(11)/ - .00000000055231518E0/
  DATA ath1cs(12)/ - .00000000020545441E0/
  DATA ath1cs(13)/ - .00000000008043412E0/
  DATA ath1cs(14)/ - .00000000003291252E0/
  DATA ath1cs(15)/ - .00000000001399875E0/
  DATA ath1cs(16)/ - .00000000000616151E0/
  DATA ath1cs(17)/ - .00000000000279614E0/
  DATA ath1cs(18)/ - .00000000000130428E0/
  DATA ath1cs(19)/ - .00000000000062373E0/
  DATA ath1cs(20)/ - .00000000000030512E0/
  DATA ath1cs(21)/ - .00000000000015239E0/
  DATA ath1cs(22)/ - .00000000000007758E0/
  DATA ath1cs(23)/ - .00000000000004020E0/
  DATA ath1cs(24)/ - .00000000000002117E0/
  DATA ath1cs(25)/ - .00000000000001132E0/
  DATA ath1cs(26)/ - .00000000000000614E0/
  DATA ath1cs(27)/ - .00000000000000337E0/
  DATA ath1cs(28)/ - .00000000000000188E0/
  DATA ath1cs(29)/ - .00000000000000105E0/
  DATA ath1cs(30)/ - .00000000000000060E0/
  DATA ath1cs(31)/ - .00000000000000034E0/
  DATA ath1cs(32)/ - .00000000000000020E0/
  DATA ath1cs(33)/ - .00000000000000011E0/
  DATA ath1cs(34)/ - .00000000000000007E0/
  DATA ath1cs(35)/ - .00000000000000004E0/
  DATA ath1cs(36)/ - .00000000000000002E0/
  DATA am22cs(1)/ - .01562844480625341E0/
  DATA am22cs(2)/.00778336445239681E0/
  DATA am22cs(3)/.00086705777047718E0/
  DATA am22cs(4)/.00015696627315611E0/
  DATA am22cs(5)/.00003563962571432E0/
  DATA am22cs(6)/.00000924598335425E0/
  DATA am22cs(7)/.00000262110161850E0/
  DATA am22cs(8)/.00000079188221651E0/
  DATA am22cs(9)/.00000025104152792E0/
  DATA am22cs(10)/.00000008265223206E0/
  DATA am22cs(11)/.00000002805711662E0/
  DATA am22cs(12)/.00000000976821090E0/
  DATA am22cs(13)/.00000000347407923E0/
  DATA am22cs(14)/.00000000125828132E0/
  DATA am22cs(15)/.00000000046298826E0/
  DATA am22cs(16)/.00000000017272825E0/
  DATA am22cs(17)/.00000000006523192E0/
  DATA am22cs(18)/.00000000002490471E0/
  DATA am22cs(19)/.00000000000960156E0/
  DATA am22cs(20)/.00000000000373448E0/
  DATA am22cs(21)/.00000000000146417E0/
  DATA am22cs(22)/.00000000000057826E0/
  DATA am22cs(23)/.00000000000022991E0/
  DATA am22cs(24)/.00000000000009197E0/
  DATA am22cs(25)/.00000000000003700E0/
  DATA am22cs(26)/.00000000000001496E0/
  DATA am22cs(27)/.00000000000000608E0/
  DATA am22cs(28)/.00000000000000248E0/
  DATA am22cs(29)/.00000000000000101E0/
  DATA am22cs(30)/.00000000000000041E0/
  DATA am22cs(31)/.00000000000000017E0/
  DATA am22cs(32)/.00000000000000007E0/
  DATA am22cs(33)/.00000000000000002E0/
  DATA ath2cs(1)/.00440527345871877E0/
  DATA ath2cs(2)/ - .03042919452318455E0/
  DATA ath2cs(3)/ - .00138565328377179E0/
  DATA ath2cs(4)/ - .00018044439089549E0/
  DATA ath2cs(5)/ - .00003380847108327E0/
  DATA ath2cs(6)/ - .00000767818353522E0/
  DATA ath2cs(7)/ - .00000196783944371E0/
  DATA ath2cs(8)/ - .00000054837271158E0/
  DATA ath2cs(9)/ - .00000016254615505E0/
  DATA ath2cs(10)/ - .00000005053049981E0/
  DATA ath2cs(11)/ - .00000001631580701E0/
  DATA ath2cs(12)/ - .00000000543420411E0/
  DATA ath2cs(13)/ - .00000000185739855E0/
  DATA ath2cs(14)/ - .00000000064895120E0/
  DATA ath2cs(15)/ - .00000000023105948E0/
  DATA ath2cs(16)/ - .00000000008363282E0/
  DATA ath2cs(17)/ - .00000000003071196E0/
  DATA ath2cs(18)/ - .00000000001142367E0/
  DATA ath2cs(19)/ - .00000000000429811E0/
  DATA ath2cs(20)/ - .00000000000163389E0/
  DATA ath2cs(21)/ - .00000000000062693E0/
  DATA ath2cs(22)/ - .00000000000024260E0/
  DATA ath2cs(23)/ - .00000000000009461E0/
  DATA ath2cs(24)/ - .00000000000003716E0/
  DATA ath2cs(25)/ - .00000000000001469E0/
  DATA ath2cs(26)/ - .00000000000000584E0/
  DATA ath2cs(27)/ - .00000000000000233E0/
  DATA ath2cs(28)/ - .00000000000000093E0/
  DATA ath2cs(29)/ - .00000000000000037E0/
  DATA ath2cs(30)/ - .00000000000000015E0/
  DATA ath2cs(31)/ - .00000000000000006E0/
  DATA ath2cs(32)/ - .00000000000000002E0/
  DATA pi4/0.78539816339744831E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  R9AIMP
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    nam21 = INITS(am21cs,40,eta)
    nath1 = INITS(ath1cs,36,eta)
    nam22 = INITS(am22cs,33,eta)
    nath2 = INITS(ath2cs,32,eta)
    !
    xsml = -1.0/R1MACH(3)**0.3333
  ENDIF
  first = .FALSE.
  !
  IF ( X>=(-2.0) ) THEN
    !
    IF ( X>(-1.0) ) CALL XERMSG('SLATEC','R9AIMP','X MUST BE LE -1.0',1,2)
    !
    z = (16.0/X**3+9.0)/7.0
    Ampl = 0.3125 + CSEVL(z,am22cs,nam22)
    Theta = -0.625 + CSEVL(z,ath2cs,nath2)
  ELSE
    z = 1.0
    IF ( X>xsml ) z = 16.0/X**3 + 1.0
    Ampl = 0.3125 + CSEVL(z,am21cs,nam21)
    Theta = -0.625 + CSEVL(z,ath1cs,nath1)
  ENDIF
  !
  sqrtx = SQRT(-X)
  Ampl = SQRT(Ampl/sqrtx)
  Theta = pi4 - X*sqrtx*Theta
  !
END SUBROUTINE R9AIMP
