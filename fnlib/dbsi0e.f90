!*==DBSI0E.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBSI0E
REAL(8) FUNCTION DBSI0E(X)
  IMPLICIT NONE
  !*--DBSI0E5
  !*** Start of declarations inserted by SPAG
  REAL eta
  INTEGER INITDS , ntai0 , ntai02 , nti0
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBSI0E
  !***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
  !            Bessel function of the first kind of order zero.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B1
  !***TYPE      DOUBLE PRECISION (BESI0E-S, DBSI0E-D)
  !***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
  !             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
  !             ORDER ZERO, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBSI0E(X) calculates the double precision exponentially scaled
  ! modified (hyperbolic) Bessel function of the first kind of order
  ! zero for double precision argument X.  The result is the Bessel
  ! function I0(X) multiplied by EXP(-ABS(X)).
  !
  ! Series for BI0        on the interval  0.          to  9.00000E+00
  !                                        with weighted error   9.51E-34
  !                                         log weighted error  33.02
  !                               significant figures required  33.31
  !                                    decimal places required  33.65
  !
  ! Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
  !                                        with weighted error   2.74E-32
  !                                         log weighted error  31.56
  !                               significant figures required  30.15
  !                                    decimal places required  32.39
  !
  ! Series for AI02       on the interval  0.          to  1.25000E-01
  !                                        with weighted error   1.97E-32
  !                                         log weighted error  31.71
  !                               significant figures required  30.15
  !                                    decimal places required  32.63
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DBSI0E
  REAL(8) :: X , bi0cs(18) , ai0cs(46) , ai02cs(69) , xsml , y , &
    D1MACH , DCSEVL
  LOGICAL first
  SAVE bi0cs , ai0cs , ai02cs , nti0 , ntai0 , ntai02 , xsml , first
  DATA bi0cs(1)/ - .7660547252839144951081894976243285D-1/
  DATA bi0cs(2)/ + .1927337953993808269952408750881196D+1/
  DATA bi0cs(3)/ + .2282644586920301338937029292330415D+0/
  DATA bi0cs(4)/ + .1304891466707290428079334210691888D-1/
  DATA bi0cs(5)/ + .4344270900816487451378682681026107D-3/
  DATA bi0cs(6)/ + .9422657686001934663923171744118766D-5/
  DATA bi0cs(7)/ + .1434006289510691079962091878179957D-6/
  DATA bi0cs(8)/ + .1613849069661749069915419719994611D-8/
  DATA bi0cs(9)/ + .1396650044535669699495092708142522D-10/
  DATA bi0cs(10)/ + .9579451725505445344627523171893333D-13/
  DATA bi0cs(11)/ + .5333981859862502131015107744000000D-15/
  DATA bi0cs(12)/ + .2458716088437470774696785919999999D-17/
  DATA bi0cs(13)/ + .9535680890248770026944341333333333D-20/
  DATA bi0cs(14)/ + .3154382039721427336789333333333333D-22/
  DATA bi0cs(15)/ + .9004564101094637431466666666666666D-25/
  DATA bi0cs(16)/ + .2240647369123670016000000000000000D-27/
  DATA bi0cs(17)/ + .4903034603242837333333333333333333D-30/
  DATA bi0cs(18)/ + .9508172606122666666666666666666666D-33/
  DATA ai0cs(1)/ + .7575994494023795942729872037438D-1/
  DATA ai0cs(2)/ + .7591380810823345507292978733204D-2/
  DATA ai0cs(3)/ + .4153131338923750501863197491382D-3/
  DATA ai0cs(4)/ + .1070076463439073073582429702170D-4/
  DATA ai0cs(5)/ - .7901179979212894660750319485730D-5/
  DATA ai0cs(6)/ - .7826143501438752269788989806909D-6/
  DATA ai0cs(7)/ + .2783849942948870806381185389857D-6/
  DATA ai0cs(8)/ + .8252472600612027191966829133198D-8/
  DATA ai0cs(9)/ - .1204463945520199179054960891103D-7/
  DATA ai0cs(10)/ + .1559648598506076443612287527928D-8/
  DATA ai0cs(11)/ + .2292556367103316543477254802857D-9/
  DATA ai0cs(12)/ - .1191622884279064603677774234478D-9/
  DATA ai0cs(13)/ + .1757854916032409830218331247743D-10/
  DATA ai0cs(14)/ + .1128224463218900517144411356824D-11/
  DATA ai0cs(15)/ - .1146848625927298877729633876982D-11/
  DATA ai0cs(16)/ + .2715592054803662872643651921606D-12/
  DATA ai0cs(17)/ - .2415874666562687838442475720281D-13/
  DATA ai0cs(18)/ - .6084469888255125064606099639224D-14/
  DATA ai0cs(19)/ + .3145705077175477293708360267303D-14/
  DATA ai0cs(20)/ - .7172212924871187717962175059176D-15/
  DATA ai0cs(21)/ + .7874493403454103396083909603327D-16/
  DATA ai0cs(22)/ + .1004802753009462402345244571839D-16/
  DATA ai0cs(23)/ - .7566895365350534853428435888810D-17/
  DATA ai0cs(24)/ + .2150380106876119887812051287845D-17/
  DATA ai0cs(25)/ - .3754858341830874429151584452608D-18/
  DATA ai0cs(26)/ + .2354065842226992576900757105322D-19/
  DATA ai0cs(27)/ + .1114667612047928530226373355110D-19/
  DATA ai0cs(28)/ - .5398891884396990378696779322709D-20/
  DATA ai0cs(29)/ + .1439598792240752677042858404522D-20/
  DATA ai0cs(30)/ - .2591916360111093406460818401962D-21/
  DATA ai0cs(31)/ + .2238133183998583907434092298240D-22/
  DATA ai0cs(32)/ + .5250672575364771172772216831999D-23/
  DATA ai0cs(33)/ - .3249904138533230784173432285866D-23/
  DATA ai0cs(34)/ + .9924214103205037927857284710400D-24/
  DATA ai0cs(35)/ - .2164992254244669523146554299733D-24/
  DATA ai0cs(36)/ + .3233609471943594083973332991999D-25/
  DATA ai0cs(37)/ - .1184620207396742489824733866666D-26/
  DATA ai0cs(38)/ - .1281671853950498650548338687999D-26/
  DATA ai0cs(39)/ + .5827015182279390511605568853333D-27/
  DATA ai0cs(40)/ - .1668222326026109719364501503999D-27/
  DATA ai0cs(41)/ + .3625309510541569975700684800000D-28/
  DATA ai0cs(42)/ - .5733627999055713589945958399999D-29/
  DATA ai0cs(43)/ + .3736796722063098229642581333333D-30/
  DATA ai0cs(44)/ + .1602073983156851963365512533333D-30/
  DATA ai0cs(45)/ - .8700424864057229884522495999999D-31/
  DATA ai0cs(46)/ + .2741320937937481145603413333333D-31/
  DATA ai02cs(1)/ + .5449041101410883160789609622680D-1/
  DATA ai02cs(2)/ + .3369116478255694089897856629799D-2/
  DATA ai02cs(3)/ + .6889758346916823984262639143011D-4/
  DATA ai02cs(4)/ + .2891370520834756482966924023232D-5/
  DATA ai02cs(5)/ + .2048918589469063741827605340931D-6/
  DATA ai02cs(6)/ + .2266668990498178064593277431361D-7/
  DATA ai02cs(7)/ + .3396232025708386345150843969523D-8/
  DATA ai02cs(8)/ + .4940602388224969589104824497835D-9/
  DATA ai02cs(9)/ + .1188914710784643834240845251963D-10/
  DATA ai02cs(10)/ - .3149916527963241364538648629619D-10/
  DATA ai02cs(11)/ - .1321581184044771311875407399267D-10/
  DATA ai02cs(12)/ - .1794178531506806117779435740269D-11/
  DATA ai02cs(13)/ + .7180124451383666233671064293469D-12/
  DATA ai02cs(14)/ + .3852778382742142701140898017776D-12/
  DATA ai02cs(15)/ + .1540086217521409826913258233397D-13/
  DATA ai02cs(16)/ - .4150569347287222086626899720156D-13/
  DATA ai02cs(17)/ - .9554846698828307648702144943125D-14/
  DATA ai02cs(18)/ + .3811680669352622420746055355118D-14/
  DATA ai02cs(19)/ + .1772560133056526383604932666758D-14/
  DATA ai02cs(20)/ - .3425485619677219134619247903282D-15/
  DATA ai02cs(21)/ - .2827623980516583484942055937594D-15/
  DATA ai02cs(22)/ + .3461222867697461093097062508134D-16/
  DATA ai02cs(23)/ + .4465621420296759999010420542843D-16/
  DATA ai02cs(24)/ - .4830504485944182071255254037954D-17/
  DATA ai02cs(25)/ - .7233180487874753954562272409245D-17/
  DATA ai02cs(26)/ + .9921475412173698598880460939810D-18/
  DATA ai02cs(27)/ + .1193650890845982085504399499242D-17/
  DATA ai02cs(28)/ - .2488709837150807235720544916602D-18/
  DATA ai02cs(29)/ - .1938426454160905928984697811326D-18/
  DATA ai02cs(30)/ + .6444656697373443868783019493949D-19/
  DATA ai02cs(31)/ + .2886051596289224326481713830734D-19/
  DATA ai02cs(32)/ - .1601954907174971807061671562007D-19/
  DATA ai02cs(33)/ - .3270815010592314720891935674859D-20/
  DATA ai02cs(34)/ + .3686932283826409181146007239393D-20/
  DATA ai02cs(35)/ + .1268297648030950153013595297109D-22/
  DATA ai02cs(36)/ - .7549825019377273907696366644101D-21/
  DATA ai02cs(37)/ + .1502133571377835349637127890534D-21/
  DATA ai02cs(38)/ + .1265195883509648534932087992483D-21/
  DATA ai02cs(39)/ - .6100998370083680708629408916002D-22/
  DATA ai02cs(40)/ - .1268809629260128264368720959242D-22/
  DATA ai02cs(41)/ + .1661016099890741457840384874905D-22/
  DATA ai02cs(42)/ - .1585194335765885579379705048814D-23/
  DATA ai02cs(43)/ - .3302645405968217800953817667556D-23/
  DATA ai02cs(44)/ + .1313580902839239781740396231174D-23/
  DATA ai02cs(45)/ + .3689040246671156793314256372804D-24/
  DATA ai02cs(46)/ - .4210141910461689149219782472499D-24/
  DATA ai02cs(47)/ + .4791954591082865780631714013730D-25/
  DATA ai02cs(48)/ + .8459470390221821795299717074124D-25/
  DATA ai02cs(49)/ - .4039800940872832493146079371810D-25/
  DATA ai02cs(50)/ - .6434714653650431347301008504695D-26/
  DATA ai02cs(51)/ + .1225743398875665990344647369905D-25/
  DATA ai02cs(52)/ - .2934391316025708923198798211754D-26/
  DATA ai02cs(53)/ - .1961311309194982926203712057289D-26/
  DATA ai02cs(54)/ + .1503520374822193424162299003098D-26/
  DATA ai02cs(55)/ - .9588720515744826552033863882069D-28/
  DATA ai02cs(56)/ - .3483339380817045486394411085114D-27/
  DATA ai02cs(57)/ + .1690903610263043673062449607256D-27/
  DATA ai02cs(58)/ + .1982866538735603043894001157188D-28/
  DATA ai02cs(59)/ - .5317498081491816214575830025284D-28/
  DATA ai02cs(60)/ + .1803306629888392946235014503901D-28/
  DATA ai02cs(61)/ + .6213093341454893175884053112422D-29/
  DATA ai02cs(62)/ - .7692189292772161863200728066730D-29/
  DATA ai02cs(63)/ + .1858252826111702542625560165963D-29/
  DATA ai02cs(64)/ + .1237585142281395724899271545541D-29/
  DATA ai02cs(65)/ - .1102259120409223803217794787792D-29/
  DATA ai02cs(66)/ + .1886287118039704490077874479431D-30/
  DATA ai02cs(67)/ + .2160196872243658913149031414060D-30/
  DATA ai02cs(68)/ - .1605454124919743200584465949655D-30/
  DATA ai02cs(69)/ + .1965352984594290603938848073318D-31/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBSI0E
  IF ( first ) THEN
    eta = 0.1*REAL(D1MACH(3))
    nti0 = INITDS(bi0cs,18,eta)
    ntai0 = INITDS(ai0cs,46,eta)
    ntai02 = INITDS(ai02cs,69,eta)
    xsml = SQRT(4.5D0*D1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>3.0D0 ) THEN
    !
    IF ( y<=8.D0 ) DBSI0E = (0.375D0+DCSEVL((48.D0/y-11.D0)/5.D0,ai0cs,ntai0&
      ))/SQRT(y)
    IF ( y>8.D0 ) DBSI0E = (0.375D0+DCSEVL(16.D0/y-1.D0,ai02cs,ntai02))&
      /SQRT(y)
    GOTO 99999
  ENDIF
  !
  DBSI0E = 1.0D0 - X
  IF ( y>xsml ) DBSI0E = EXP(-y)*(2.75D0+DCSEVL(y*y/4.5D0-1.D0,bi0cs,nti0))
  RETURN
  !
  99999 END FUNCTION DBSI0E
