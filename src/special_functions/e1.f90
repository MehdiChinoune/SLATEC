!** E1
REAL FUNCTION E1(X)
  !>
  !  Compute the exponential integral E1(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      SINGLE PRECISION (E1-S, DE1-D)
  !***
  ! **Keywords:**  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! E1 calculates the single precision exponential integral, E1(X), for
  ! positive single precision argument X and the Cauchy principal value
  ! for negative X.  If principal values are used everywhere, then, for
  ! all X,
  !
  !    E1(X) = -Ei(-X)
  ! or
  !    Ei(X) = -E1(-X).
  !
  !
  ! Series for AE11       on the interval -1.00000D-01 to  0.
  !                                        with weighted error   1.76E-17
  !                                         log weighted error  16.75
  !                               significant figures required  15.70
  !                                    decimal places required  17.55
  !
  !
  ! Series for AE12       on the interval -2.50000D-01 to -1.00000D-01
  !                                        with weighted error   5.83E-17
  !                                         log weighted error  16.23
  !                               significant figures required  15.76
  !                                    decimal places required  16.93
  !
  !
  ! Series for E11        on the interval -4.00000D+00 to -1.00000D+00
  !                                        with weighted error   1.08E-18
  !                                         log weighted error  17.97
  !                               significant figures required  19.02
  !                                    decimal places required  18.61
  !
  !
  ! Series for E12        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   3.15E-18
  !                                         log weighted error  17.50
  !                        approx significant figures required  15.8
  !                                    decimal places required  18.10
  !
  !
  ! Series for AE13       on the interval  2.50000D-01 to  1.00000D+00
  !                                        with weighted error   2.34E-17
  !                                         log weighted error  16.63
  !                               significant figures required  16.14
  !                                    decimal places required  17.33
  !
  !
  ! Series for AE14       on the interval  0.          to  2.50000D-01
  !                                        with weighted error   5.41E-17
  !                                         log weighted error  16.27
  !                               significant figures required  15.38
  !                                    decimal places required  16.97
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891115  Modified prologue description.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL eta, X, xmaxt
  INTEGER, SAVE :: ntae11, ntae12, nte11, nte12, ntae13, ntae14
  REAL, SAVE :: xmax
  REAL, PARAMETER :: ae11cs(39) = [ .12150323971606579E0, -.065088778513550150E0, &
    .004897651357459670E0, -.000649237843027216E0, .000093840434587471E0, &
    .000000420236380882E0, -.000008113374735904E0, .000002804247688663E0, &
    .000000056487164441E0, -.000000344809174450E0, .000000058209273578E0, &
    .000000038711426349E0, -.000000012453235014E0,-.000000005118504888E0, &
    .000000002148771527E0,  .000000000868459898E0,-.000000000343650105E0, &
    -.000000000179796603E0, .000000000047442060E0, .000000000040423282E0, &
    -.000000000003543928E0,-.000000000008853444E0,-.000000000000960151E0, &
    .000000000001692921E0,  .000000000000607990E0,-.000000000000224338E0, &
    -.000000000000200327E0,-.000000000000006246E0, .000000000000045571E0, &
    .000000000000016383E0, -.000000000000005561E0,-.000000000000006074E0, &
    -.000000000000000862E0, .000000000000001223E0, .000000000000000716E0, &
    -.000000000000000024E0,-.000000000000000201E0,-.000000000000000082E0, &
    .000000000000000017E0 ]
  REAL, PARAMETER :: ae12cs(25) = [ .58241749513472674E0, -.15834885090578275E0, &
    -.006764275590323141E0, .005125843950185725E0, .000435232492169391E0, &
    -.000143613366305483E0,-.000041801320556301E0,-.000002713395758640E0, &
    .000001151381913647E0,  .000000420650022012E0, .000000066581901391E0, &
    .000000000662143777E0, -.000000002844104870E0,-.000000000940724197E0, &
    -.000000000177476602E0,-.000000000015830222E0, .000000000002905732E0, &
    .000000000001769356E0,  .000000000000492735E0, .000000000000093709E0, &
    .000000000000010707E0, -.000000000000000537E0,-.000000000000000716E0, &
    -.000000000000000244E0,-.000000000000000058E0 ]
  REAL, PARAMETER :: e11cs(19) = [ -16.113461655571494026E0, 7.7940727787426802769E0, &
    -1.9554058188631419507E0, .37337293866277945612E0,-.05692503191092901938E0, &
    .00721107776966009185E0, -.00078104901449841593E0, .00007388093356262168E0, &
    -.00000620286187580820E0, .00000046816002303176E0,-.00000003209288853329E0, &
    .00000000201519974874E0, -.00000000011673686816E0, .00000000000627627066E0, &
    -.00000000000031481541E0, .00000000000001479904E0,-.00000000000000065457E0, &
    .00000000000000002733E0, -.00000000000000000108E0 ]
  REAL, PARAMETER :: e12cs(16) =  [ -0.037390214792202795E0, 0.042723986062209577E0, &
    -.1303182079849700544E0, .01441912402469889073E0, -.00134617078051068022E0, &
    .00010731029253063780E0, -.00000742999951611943E0, .00000045377325690753E0, &
    -.00000002476417211390E0, .00000000122076581374E0,-.00000000005485141480E0, &
    .00000000000226362142E0, -.00000000000008635897E0, .00000000000000306291E0, &
    -.00000000000000010148E0, .00000000000000000315E0 ]
  REAL, PARAMETER :: ae13cs(25) = [ -.60577324664060346E0, -.11253524348366090E0, &
    .013432266247902779E0, -.001926845187381145E0, .000309118337720603E0, &
    -.000053564132129618E0, .000009827812880247E0,-.000001885368984916E0, &
    .000000374943193568E0, -.000000076823455870E0, .000000016143270567E0, &
    -.000000003466802211E0, .000000000758754209E0,-.000000000168864333E0, &
    .000000000038145706E0, -.000000000008733026E0, .000000000002023672E0, &
    -.000000000000474132E0, .000000000000112211E0,-.000000000000026804E0, &
    .000000000000006457E0, -.000000000000001568E0, .000000000000000383E0, &
    -.000000000000000094E0, .000000000000000023E0 ]
  REAL, PARAMETER :: ae14cs(26) = [ -.1892918000753017E0, -.08648117855259871E0, &
    .00722410154374659E0, -.00080975594575573E0, .00010999134432661E0, &
    -.00001717332998937E0, .00000298562751447E0,-.00000056596491457E0, &
    .00000011526808397E0, -.00000002495030440E0, .00000000569232420E0, &
    -.00000000135995766E0, .00000000033846628E0,-.00000000008737853E0, &
    .00000000002331588E0, -.00000000000641148E0, .00000000000181224E0, &
    -.00000000000052538E0, .00000000000015592E0,-.00000000000004729E0, &
    .00000000000001463E0, -.00000000000000461E0, .00000000000000148E0, &
    -.00000000000000048E0, .00000000000000016E0,-.00000000000000005E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  E1
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    ntae11 = INITS(ae11cs,39,eta)
    ntae12 = INITS(ae12cs,25,eta)
    nte11 = INITS(e11cs,19,eta)
    nte12 = INITS(e12cs,16,eta)
    ntae13 = INITS(ae13cs,25,eta)
    ntae14 = INITS(ae14cs,26,eta)
    !
    xmaxt = -LOG(R1MACH(1))
    xmax = xmaxt - LOG(xmaxt)
    first = .FALSE.
  END IF
  !
  IF ( X<=(-10.) ) THEN
    !
    ! E1(X) = -EI(-X) FOR X .LE. -10.
    !
    E1 = EXP(-X)/X*(1.+CSEVL(20./X+1.,ae11cs,ntae11))
    RETURN
    !
  ELSEIF ( X<=(-4.0) ) THEN
    !
    ! E1(X) = -EI(-X) FOR -10. .LT. X .LE. -4.
    !
    E1 = EXP(-X)/X*(1.+CSEVL((40./X+7.)/3.,ae12cs,ntae12))
    RETURN
    !
  ELSEIF ( X<=(-1.0) ) THEN
    !
    ! E1(X) = -EI(-X) FOR -4. .LT. X .LE. -1.
    !
    E1 = -LOG(ABS(X)) + CSEVL((2.*X+5.)/3.,e11cs,nte11)
    RETURN
    !
  ELSEIF ( X<=1. ) THEN
    IF ( X==0. ) CALL XERMSG('SLATEC','E1','X IS 0',2,2)
    !
    ! E1(X) = -EI(-X) FOR -1. .LT. X .LE. 1.,  X .NE. 0.
    !
    E1 = (-LOG(ABS(X))-0.6875+X) + CSEVL(X,e12cs,nte12)
    RETURN
    !
  ELSEIF ( X<=4. ) THEN
    !
    ! E1(X) = -EI(-X) FOR 1. .LT. X .LE. 4.
    !
    E1 = EXP(-X)/X*(1.+CSEVL((8./X-5.)/3.,ae13cs,ntae13))
    RETURN
    !
  ELSEIF ( X>xmax ) THEN
    !
    ! E1(X) = -EI(-X) FOR X .GT. XMAX
    !
    CALL XERMSG('SLATEC','E1','X SO BIG E1 UNDERFLOWS',1,1)
    E1 = 0.
    RETURN
  END IF
  !
  ! E1(X) = -EI(-X) FOR 4. .LT. X .LE. XMAX
  !
  E1 = EXP(-X)/X*(1.+CSEVL(8./X-1.,ae14cs,ntae14))
  RETURN
END FUNCTION E1
