!*==ERFC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK ERFC
FUNCTION ERFC(X)
  IMPLICIT NONE
  !*--ERFC5
  !*** Start of declarations inserted by SPAG
  REAL CSEVL, erc2cs, ERFC, erfccs, erfcs, eta, R1MACH, sqeps, &
    sqrtpi, txmax, X, xmax, xsml, y
  INTEGER INITS, nterc2, nterf, nterfc
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  ERFC
  !***PURPOSE  Compute the complementary error function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C8A, L5A1E
  !***TYPE      SINGLE PRECISION (ERFC-S, DERFC-D)
  !***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ERFC(X) calculates the single precision complementary error
  ! function for single precision argument X.
  !
  ! Series for ERF        on the interval  0.          to  1.00000D+00
  !                                        with weighted error   7.10E-18
  !                                         log weighted error  17.15
  !                               significant figures required  16.31
  !                                    decimal places required  17.71
  !
  ! Series for ERFC       on the interval  0.          to  2.50000D-01
  !                                        with weighted error   4.81E-17
  !                                         log weighted error  16.32
  !                        approx significant figures required  15.0
  !
  !
  ! Series for ERC2       on the interval  2.50000D-01 to  1.00000D+00
  !                                        with weighted error   5.22E-17
  !                                         log weighted error  16.28
  !                        approx significant figures required  15.0
  !                                    decimal places required  16.96
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  !***END PROLOGUE  ERFC
  DIMENSION erfcs(13), erfccs(24), erc2cs(23)
  LOGICAL first
  SAVE erfcs, erc2cs, erfccs, sqrtpi, nterf, nterfc, nterc2, xsml, &
    xmax, sqeps, first
  DATA erfcs(1)/ - .049046121234691808E0/
  DATA erfcs(2)/ - .14226120510371364E0/
  DATA erfcs(3)/.010035582187599796E0/
  DATA erfcs(4)/ - .000576876469976748E0/
  DATA erfcs(5)/.000027419931252196E0/
  DATA erfcs(6)/ - .000001104317550734E0/
  DATA erfcs(7)/.000000038488755420E0/
  DATA erfcs(8)/ - .000000001180858253E0/
  DATA erfcs(9)/.000000000032334215E0/
  DATA erfcs(10)/ - .000000000000799101E0/
  DATA erfcs(11)/.000000000000017990E0/
  DATA erfcs(12)/ - .000000000000000371E0/
  DATA erfcs(13)/.000000000000000007E0/
  DATA erc2cs(1)/ - .069601346602309501E0/
  DATA erc2cs(2)/ - .041101339362620893E0/
  DATA erc2cs(3)/.003914495866689626E0/
  DATA erc2cs(4)/ - .000490639565054897E0/
  DATA erc2cs(5)/.000071574790013770E0/
  DATA erc2cs(6)/ - .000011530716341312E0/
  DATA erc2cs(7)/.000001994670590201E0/
  DATA erc2cs(8)/ - .000000364266647159E0/
  DATA erc2cs(9)/.000000069443726100E0/
  DATA erc2cs(10)/ - .000000013712209021E0/
  DATA erc2cs(11)/.000000002788389661E0/
  DATA erc2cs(12)/ - .000000000581416472E0/
  DATA erc2cs(13)/.000000000123892049E0/
  DATA erc2cs(14)/ - .000000000026906391E0/
  DATA erc2cs(15)/.000000000005942614E0/
  DATA erc2cs(16)/ - .000000000001332386E0/
  DATA erc2cs(17)/.000000000000302804E0/
  DATA erc2cs(18)/ - .000000000000069666E0/
  DATA erc2cs(19)/.000000000000016208E0/
  DATA erc2cs(20)/ - .000000000000003809E0/
  DATA erc2cs(21)/.000000000000000904E0/
  DATA erc2cs(22)/ - .000000000000000216E0/
  DATA erc2cs(23)/.000000000000000052E0/
  DATA erfccs(1)/0.0715179310202925E0/
  DATA erfccs(2)/ - .026532434337606719E0/
  DATA erfccs(3)/.001711153977920853E0/
  DATA erfccs(4)/ - .000163751663458512E0/
  DATA erfccs(5)/.000019871293500549E0/
  DATA erfccs(6)/ - .000002843712412769E0/
  DATA erfccs(7)/.000000460616130901E0/
  DATA erfccs(8)/ - .000000082277530261E0/
  DATA erfccs(9)/.000000015921418724E0/
  DATA erfccs(10)/ - .000000003295071356E0/
  DATA erfccs(11)/.000000000722343973E0/
  DATA erfccs(12)/ - .000000000166485584E0/
  DATA erfccs(13)/.000000000040103931E0/
  DATA erfccs(14)/ - .000000000010048164E0/
  DATA erfccs(15)/.000000000002608272E0/
  DATA erfccs(16)/ - .000000000000699105E0/
  DATA erfccs(17)/.000000000000192946E0/
  DATA erfccs(18)/ - .000000000000054704E0/
  DATA erfccs(19)/.000000000000015901E0/
  DATA erfccs(20)/ - .000000000000004729E0/
  DATA erfccs(21)/.000000000000001432E0/
  DATA erfccs(22)/ - .000000000000000439E0/
  DATA erfccs(23)/.000000000000000138E0/
  DATA erfccs(24)/ - .000000000000000048E0/
  DATA sqrtpi/1.7724538509055160E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  ERFC
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    nterf = INITS(erfcs,13,eta)
    nterfc = INITS(erfccs,24,eta)
    nterc2 = INITS(erc2cs,23,eta)
    !
    xsml = -SQRT(-LOG(sqrtpi*R1MACH(3)))
    txmax = SQRT(-LOG(sqrtpi*R1MACH(1)))
    xmax = txmax - 0.5*LOG(txmax)/txmax - 0.01
    sqeps = SQRT(2.0*R1MACH(3))
  ENDIF
  first = .FALSE.
  !
  IF ( X<=xsml ) THEN
    !
    ! ERFC(X) = 1.0 - ERF(X) FOR X .LT. XSML
    !
    ERFC = 2.
    RETURN
    !
  ELSEIF ( X>xmax ) THEN
    !
    CALL XERMSG('SLATEC','ERFC','X SO BIG ERFC UNDERFLOWS',1,1)
    ERFC = 0.
    GOTO 99999
  ELSE
    y = ABS(X)
    IF ( y<=1.0 ) THEN
      !
      ! ERFC(X) = 1.0 - ERF(X) FOR -1. .LE. X .LE. 1.
      !
      IF ( y<sqeps ) ERFC = 1.0 - 2.0*X/sqrtpi
      IF ( y>=sqeps ) ERFC = 1.0 - X*(1.0+CSEVL(2.*X*X-1.,erfcs,nterf))
      RETURN
    ENDIF
  ENDIF
  !
  ! ERFC(X) = 1.0 - ERF(X) FOR 1. .LT. ABS(X) .LE. XMAX
  !
  y = y*y
  IF ( y<=4. ) ERFC = EXP(-y)/ABS(X)*(0.5+CSEVL((8./y-5.)/3.,erc2cs,nterc2))
  IF ( y>4. ) ERFC = EXP(-y)/ABS(X)*(0.5+CSEVL(8./y-1.,erfccs,nterfc))
  IF ( X<0. ) ERFC = 2.0 - ERFC
  RETURN
  !
  99999 CONTINUE
  END FUNCTION ERFC
