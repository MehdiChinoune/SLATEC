!*==D9LN2R.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9LN2R
DOUBLE PRECISION FUNCTION D9LN2R(X)
  IMPLICIT NONE
  !*--D9LN2R5
  !*** Start of declarations inserted by SPAG
  REAL eps , sqeps
  INTEGER INITDS , ntln21 , ntln22
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9LN2R
  !***SUBSIDIARY
  !***PURPOSE  Evaluate LOG(1+X) from second order relative accuracy so
  !            that LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X)
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4B
  !***TYPE      DOUBLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so
  ! that    LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X)
  !
  ! Series for LN21       on the interval -6.25000E-01 to  0.
  !                                        with weighted error   1.82E-32
  !                                         log weighted error  31.74
  !                               significant figures required  31.00
  !                                    decimal places required  32.59
  !
  ! Series for LN22       on the interval  0.          to  8.12500E-01
  !                                        with weighted error   6.10E-32
  !                                         log weighted error  31.21
  !                               significant figures required  30.32
  !                                    decimal places required  32.00
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  D9LN2R
  DOUBLE PRECISION X , xbig , txbig , xmax , txmax , xmin , ln21cs(50) , &
    ln22cs(37) , DCSEVL , D1MACH
  LOGICAL first
  SAVE ln21cs , ln22cs , ntln21 , ntln22 , xmin , xbig , xmax , first
  DATA ln21cs(1)/ + .18111962513478809875894953043071D+0/
  DATA ln21cs(2)/ - .15627123192872462669625155541078D+0/
  DATA ln21cs(3)/ + .28676305361557275209540627102051D-1/
  DATA ln21cs(4)/ - .55586996559481398781157725126781D-2/
  DATA ln21cs(5)/ + .11178976652299837657335666279727D-2/
  DATA ln21cs(6)/ - .23080508982327947182299279585705D-3/
  DATA ln21cs(7)/ + .48598853341100175874681558068750D-4/
  DATA ln21cs(8)/ - .10390127388903210765514242633338D-4/
  DATA ln21cs(9)/ + .22484563707390128494621804946408D-5/
  DATA ln21cs(10)/ - .49140592739266484875327802597091D-6/
  DATA ln21cs(11)/ + .10828256507077483336620152971597D-6/
  DATA ln21cs(12)/ - .24025872763420701435976675416719D-7/
  DATA ln21cs(13)/ + .53624600472708133762984443250163D-8/
  DATA ln21cs(14)/ - .12029951362138772264671646424377D-8/
  DATA ln21cs(15)/ + .27107889277591860785622551632266D-9/
  DATA ln21cs(16)/ - .61323562618319010068796728430690D-10/
  DATA ln21cs(17)/ + .13920858369159469857436908543978D-10/
  DATA ln21cs(18)/ - .31699300330223494015283057260883D-11/
  DATA ln21cs(19)/ + .72383754044307505335214326197011D-12/
  DATA ln21cs(20)/ - .16570017184764411391498805506268D-12/
  DATA ln21cs(21)/ + .38018428663117424257364422631876D-13/
  DATA ln21cs(22)/ - .87411189296972700259724429899137D-14/
  DATA ln21cs(23)/ + .20135619845055748302118751028154D-14/
  DATA ln21cs(24)/ - .46464456409033907031102008154477D-15/
  DATA ln21cs(25)/ + .10739282147018339453453338554925D-15/
  DATA ln21cs(26)/ - .24858534619937794755534021833960D-16/
  DATA ln21cs(27)/ + .57620197950800189813888142628181D-17/
  DATA ln21cs(28)/ - .13373063769804394701402199958050D-17/
  DATA ln21cs(29)/ + .31074653227331824966533807166805D-18/
  DATA ln21cs(30)/ - .72288104083040539906901957917627D-19/
  DATA ln21cs(31)/ + .16833783788037385103313258186888D-19/
  DATA ln21cs(32)/ - .39239463312069958052519372739925D-20/
  DATA ln21cs(33)/ + .91551468387536789746385528640853D-21/
  DATA ln21cs(34)/ - .21378895321320159520982095801002D-21/
  DATA ln21cs(35)/ + .49964507479047864699828564568746D-22/
  DATA ln21cs(36)/ - .11686240636080170135360806147413D-22/
  DATA ln21cs(37)/ + .27353123470391863775628686786559D-23/
  DATA ln21cs(38)/ - .64068025084792111965050345881599D-24/
  DATA ln21cs(39)/ + .15016293204334124162949071940266D-24/
  DATA ln21cs(40)/ - .35217372410398479759497145002666D-25/
  DATA ln21cs(41)/ + .82643901014814767012482733397333D-26/
  DATA ln21cs(42)/ - .19404930275943401918036617898666D-26/
  DATA ln21cs(43)/ + .45587880018841283562451588437333D-27/
  DATA ln21cs(44)/ - .10715492087545202154378625023999D-27/
  DATA ln21cs(45)/ + .25199408007927592978096674133333D-28/
  DATA ln21cs(46)/ - .59289088400120969341750476800000D-29/
  DATA ln21cs(47)/ + .13955864061057513058237153279999D-29/
  DATA ln21cs(48)/ - .32864578813478583431436697599999D-30/
  DATA ln21cs(49)/ + .77424967950478166247254698666666D-31/
  DATA ln21cs(50)/ - .18247735667260887638125226666666D-31/
  DATA ln22cs(1)/ - .2224253253502046082986015223552D+0/
  DATA ln22cs(2)/ - .6104710010807862398680104755764D-1/
  DATA ln22cs(3)/ + .7427235009750394590519629755729D-2/
  DATA ln22cs(4)/ - .9335018261636970565612779606397D-3/
  DATA ln22cs(5)/ + .1200499076872601283350731287359D-3/
  DATA ln22cs(6)/ - .1570472295282004112823352608243D-4/
  DATA ln22cs(7)/ + .2081874781051271096050783592759D-5/
  DATA ln22cs(8)/ - .2789195577646713654057213051375D-6/
  DATA ln22cs(9)/ + .3769355823760132058422895135447D-7/
  DATA ln22cs(10)/ - .5130902896527711258240589938003D-8/
  DATA ln22cs(11)/ + .7027141178150694738206218215392D-9/
  DATA ln22cs(12)/ - .9674859550134342389243972005137D-10/
  DATA ln22cs(13)/ + .1338104645924887306588496449748D-10/
  DATA ln22cs(14)/ - .1858102603534063981628453846591D-11/
  DATA ln22cs(15)/ + .2589294422527919749308600123070D-12/
  DATA ln22cs(16)/ - .3619568316141588674466025382172D-13/
  DATA ln22cs(17)/ + .5074037398016623088006858917396D-14/
  DATA ln22cs(18)/ - .7131012977031127302700938748927D-15/
  DATA ln22cs(19)/ + .1004490328554567481853386784126D-15/
  DATA ln22cs(20)/ - .1417906532184025791904405075285D-16/
  DATA ln22cs(21)/ + .2005297034743326117891086396074D-17/
  DATA ln22cs(22)/ - .2840996662339803305365396717567D-18/
  DATA ln22cs(23)/ + .4031469883969079899599878662826D-19/
  DATA ln22cs(24)/ - .5729325241832207320455498956799D-20/
  DATA ln22cs(25)/ + .8153488253890010675848928733866D-21/
  DATA ln22cs(26)/ - .1161825588549721787606027468799D-21/
  DATA ln22cs(27)/ + .1657516611662538343659339775999D-22/
  DATA ln22cs(28)/ - .2367336704710805190114017280000D-23/
  DATA ln22cs(29)/ + .3384670367975521386076569599999D-24/
  DATA ln22cs(30)/ - .4843940829215718204296396799999D-25/
  DATA ln22cs(31)/ + .6938759162514273718676138666666D-26/
  DATA ln22cs(32)/ - .9948142607031436571923797333333D-27/
  DATA ln22cs(33)/ + .1427440611211698610634752000000D-27/
  DATA ln22cs(34)/ - .2049794721898234911566506666666D-28/
  DATA ln22cs(35)/ + .2945648756401362222885546666666D-29/
  DATA ln22cs(36)/ - .4235973185184957027669333333333D-30/
  DATA ln22cs(37)/ + .6095532614003832040106666666666D-31/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9LN2R
  IF ( first ) THEN
    eps = D1MACH(3)
    ntln21 = INITDS(ln21cs,50,0.1*eps)
    ntln22 = INITDS(ln22cs,37,0.1*eps)
    !
    xmin = -1.0D0 + SQRT(D1MACH(4))
    sqeps = SQRT(eps)
    txmax = 8.0/sqeps
    xmax = txmax - (eps*txmax**2-2.D0*LOG(txmax))/(2.D0*eps*txmax)
    txbig = 6.0/SQRT(sqeps)
    xbig = txbig - (sqeps*txbig**2-2.D0*LOG(txbig))/(2.D0*sqeps*txbig)
  ENDIF
  first = .FALSE.
  !
  IF ( X<(-.625D0).OR.X>0.8125D0 ) THEN
    !
    IF ( X<xmin ) CALL XERMSG('SLATEC','D9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1'&
      ,1,1)
    IF ( X>xmax ) CALL XERMSG('SLATEC','D9LN2R',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',&
      3,2)
    IF ( X>xbig ) CALL XERMSG('SLATEC','D9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'&
      ,2,1)
    !
    D9LN2R = (LOG(1.D0+X)-X*(1.D0-0.5D0*X))/X**3
    GOTO 99999
  ENDIF
  !
  IF ( X<0.0D0 ) D9LN2R = 0.375D0 + DCSEVL(16.D0*X/5.D0+1.D0,ln21cs,ntln21)
  IF ( X>=0.0D0 ) D9LN2R = 0.375D0 + DCSEVL(32.D0*X/13.D0-1.D0,ln22cs,&
    ntln22)
  RETURN
  !
  99999 END FUNCTION D9LN2R
