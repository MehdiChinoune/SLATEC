!** E1
REAL(SP) ELEMENTAL FUNCTION E1(X)
  !> Compute the exponential integral E1(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      SINGLE PRECISION (E1-S, DE1-D)
  !***
  ! **Keywords:**  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB, SPECIAL FUNCTIONS
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
  USE service, ONLY : tiny_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  INTEGER, PARAMETER :: ntae11 = 14, ntae12 = 11, nte11 = 11, nte12 = 9, &
    ntae13 = 11, ntae14 = 11
  REAL(SP), PARAMETER :: xmaxt = -LOG(tiny_sp), xmax = xmaxt - LOG(xmaxt)
  REAL(SP), PARAMETER :: ae11cs(39) = [ .12150323971606579_SP, -.065088778513550150_SP, &
    .004897651357459670_SP, -.000649237843027216_SP, .000093840434587471_SP, &
    .000000420236380882_SP, -.000008113374735904_SP, .000002804247688663_SP, &
    .000000056487164441_SP, -.000000344809174450_SP, .000000058209273578_SP, &
    .000000038711426349_SP, -.000000012453235014_SP,-.000000005118504888_SP, &
    .000000002148771527_SP,  .000000000868459898_SP,-.000000000343650105_SP, &
    -.000000000179796603_SP, .000000000047442060_SP, .000000000040423282_SP, &
    -.000000000003543928_SP,-.000000000008853444_SP,-.000000000000960151_SP, &
    .000000000001692921_SP,  .000000000000607990_SP,-.000000000000224338_SP, &
    -.000000000000200327_SP,-.000000000000006246_SP, .000000000000045571_SP, &
    .000000000000016383_SP, -.000000000000005561_SP,-.000000000000006074_SP, &
    -.000000000000000862_SP, .000000000000001223_SP, .000000000000000716_SP, &
    -.000000000000000024_SP,-.000000000000000201_SP,-.000000000000000082_SP, &
    .000000000000000017_SP ]
  REAL(SP), PARAMETER :: ae12cs(25) = [ .58241749513472674_SP, -.15834885090578275_SP, &
    -.006764275590323141_SP, .005125843950185725_SP, .000435232492169391_SP, &
    -.000143613366305483_SP,-.000041801320556301_SP,-.000002713395758640_SP, &
    .000001151381913647_SP,  .000000420650022012_SP, .000000066581901391_SP, &
    .000000000662143777_SP, -.000000002844104870_SP,-.000000000940724197_SP, &
    -.000000000177476602_SP,-.000000000015830222_SP, .000000000002905732_SP, &
    .000000000001769356_SP,  .000000000000492735_SP, .000000000000093709_SP, &
    .000000000000010707_SP, -.000000000000000537_SP,-.000000000000000716_SP, &
    -.000000000000000244_SP,-.000000000000000058_SP ]
  REAL(SP), PARAMETER :: e11cs(19) = [ -16.113461655571494026_SP, 7.7940727787426802769_SP, &
    -1.9554058188631419507_SP, .37337293866277945612_SP,-.05692503191092901938_SP, &
    .00721107776966009185_SP, -.00078104901449841593_SP, .00007388093356262168_SP, &
    -.00000620286187580820_SP, .00000046816002303176_SP,-.00000003209288853329_SP, &
    .00000000201519974874_SP, -.00000000011673686816_SP, .00000000000627627066_SP, &
    -.00000000000031481541_SP, .00000000000001479904_SP,-.00000000000000065457_SP, &
    .00000000000000002733_SP, -.00000000000000000108_SP ]
  REAL(SP), PARAMETER :: e12cs(16) =  [ -0.037390214792202795_SP, 0.042723986062209577_SP, &
    -.1303182079849700544_SP, .01441912402469889073_SP, -.00134617078051068022_SP, &
    .00010731029253063780_SP, -.00000742999951611943_SP, .00000045377325690753_SP, &
    -.00000002476417211390_SP, .00000000122076581374_SP,-.00000000005485141480_SP, &
    .00000000000226362142_SP, -.00000000000008635897_SP, .00000000000000306291_SP, &
    -.00000000000000010148_SP, .00000000000000000315_SP ]
  REAL(SP), PARAMETER :: ae13cs(25) = [ -.60577324664060346_SP, -.11253524348366090_SP, &
    .013432266247902779_SP, -.001926845187381145_SP, .000309118337720603_SP, &
    -.000053564132129618_SP, .000009827812880247_SP,-.000001885368984916_SP, &
    .000000374943193568_SP, -.000000076823455870_SP, .000000016143270567_SP, &
    -.000000003466802211_SP, .000000000758754209_SP,-.000000000168864333_SP, &
    .000000000038145706_SP, -.000000000008733026_SP, .000000000002023672_SP, &
    -.000000000000474132_SP, .000000000000112211_SP,-.000000000000026804_SP, &
    .000000000000006457_SP, -.000000000000001568_SP, .000000000000000383_SP, &
    -.000000000000000094_SP, .000000000000000023_SP ]
  REAL(SP), PARAMETER :: ae14cs(26) = [ -.1892918000753017_SP, -.08648117855259871_SP, &
    .00722410154374659_SP, -.00080975594575573_SP, .00010999134432661_SP, &
    -.00001717332998937_SP, .00000298562751447_SP,-.00000056596491457_SP, &
    .00000011526808397_SP, -.00000002495030440_SP, .00000000569232420_SP, &
    -.00000000135995766_SP, .00000000033846628_SP,-.00000000008737853_SP, &
    .00000000002331588_SP, -.00000000000641148_SP, .00000000000181224_SP, &
    -.00000000000052538_SP, .00000000000015592_SP,-.00000000000004729_SP, &
    .00000000000001463_SP, -.00000000000000461_SP, .00000000000000148_SP, &
    -.00000000000000048_SP, .00000000000000016_SP,-.00000000000000005_SP ]
  !* FIRST EXECUTABLE STATEMENT  E1
  ! ntae11 = INITS(ae11cs,eta)
  ! ntae12 = INITS(ae12cs,eta)
  ! nte11 = INITS(e11cs,eta)
  ! nte12 = INITS(e12cs,eta)
  ! ntae13 = INITS(ae13cs,eta)
  ! ntae14 = INITS(ae14cs,eta)
  !
  IF( X==0. ) ERROR STOP 'E1 : X IS 0'

  IF( X<=(-10._SP) ) THEN
    ! E1(X) = -EI(-X) FOR X <= -10.
    E1 = EXP(-X)/X*(1._SP+CSEVL(20._SP/X+1._SP,ae11cs(1:ntae11)))
  ELSEIF( X<=(-4._SP) ) THEN
    ! E1(X) = -EI(-X) FOR -10. < X <= -4.
    E1 = EXP(-X)/X*(1._SP+CSEVL((40._SP/X+7._SP)/3._SP,ae12cs(1:ntae12)))
  ELSEIF( X<=(-1._SP) ) THEN
    ! E1(X) = -EI(-X) FOR -4. < X <= -1.
    E1 = -LOG(ABS(X)) + CSEVL((2._SP*X+5._SP)/3._SP,e11cs(1:nte11))
  ELSEIF( X<=1. ) THEN
    ! E1(X) = -EI(-X) FOR -1. < X <= 1.,  X /= 0.
    E1 = (-LOG(ABS(X))-0.6875_SP+X) + CSEVL(X,e12cs(1:nte12))
  ELSEIF( X<=4. ) THEN
    ! E1(X) = -EI(-X) FOR 1. < X <= 4.
    E1 = EXP(-X)/X*(1._SP+CSEVL((8._SP/X-5._SP)/3._SP,ae13cs(1:ntae13)))
  ELSEIF( X<=xmax ) THEN
    ! E1(X) = -EI(-X) FOR 4. < X <= XMAX
    E1 = EXP(-X)/X*(1._SP+CSEVL(8._SP/X-1._SP,ae14cs(1:ntae14)))
  ELSE
    ! E1(X) = -EI(-X) FOR X > XMAX
    ! CALL XERMSG('E1 : X SO BIG E1 UNDERFLOWS',1,1)
    E1 = 0.
  END IF

  RETURN
END FUNCTION E1