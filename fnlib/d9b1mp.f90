!*==D9B1MP.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9B1MP
SUBROUTINE D9B1MP(X,Ampl,Theta)
  IMPLICIT NONE
  !*--D9B1MP5
  !*** Start of declarations inserted by SPAG
  REAL eta
  INTEGER INITDS , nbm1 , nbm12 , nbt12 , nbth1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9B1MP
  !***SUBSIDIARY
  !***PURPOSE  Evaluate the modulus and phase for the J1 and Y1 Bessel
  !            functions.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10A1
  !***TYPE      DOUBLE PRECISION (D9B1MP-D)
  !***KEYWORDS  BESSEL FUNCTION, FNLIB, MODULUS, PHASE, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate the modulus and phase for the Bessel J1 and Y1 functions.
  !
  ! Series for BM1        on the interval  1.56250E-02 to  6.25000E-02
  !                                        with weighted error   4.91E-32
  !                                         log weighted error  31.31
  !                               significant figures required  30.04
  !                                    decimal places required  32.09
  !
  ! Series for BT12       on the interval  1.56250E-02 to  6.25000E-02
  !                                        with weighted error   3.33E-32
  !                                         log weighted error  31.48
  !                               significant figures required  31.05
  !                                    decimal places required  32.27
  !
  ! Series for BM12       on the interval  0.          to  1.56250E-02
  !                                        with weighted error   5.01E-32
  !                                         log weighted error  31.30
  !                               significant figures required  29.99
  !                                    decimal places required  32.10
  !
  ! Series for BTH1       on the interval  0.          to  1.56250E-02
  !                                        with weighted error   2.82E-32
  !                                         log weighted error  31.55
  !                               significant figures required  31.12
  !                                    decimal places required  32.37
  !
  !***SEE ALSO  DBESJ1, DBESY1
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !   920618  Removed space from variable name and code restructured to
  !           use IF-THEN-ELSE.  (RWC, WRB)
  !***END PROLOGUE  D9B1MP
  REAL(8) :: X , Ampl , Theta , bm1cs(37) , bt12cs(39) , bm12cs(40) , &
    bth1cs(44) , xmax , pi4 , z , D1MACH , DCSEVL
  LOGICAL first
  SAVE bm1cs , bt12cs , bth1cs , bm12cs , pi4 , nbm1 , nbt12 , nbm12 , &
    nbth1 , xmax , first
  DATA bm1cs(1)/ + .1069845452618063014969985308538D+0/
  DATA bm1cs(2)/ + .3274915039715964900729055143445D-2/
  DATA bm1cs(3)/ - .2987783266831698592030445777938D-4/
  DATA bm1cs(4)/ + .8331237177991974531393222669023D-6/
  DATA bm1cs(5)/ - .4112665690302007304896381725498D-7/
  DATA bm1cs(6)/ + .2855344228789215220719757663161D-8/
  DATA bm1cs(7)/ - .2485408305415623878060026596055D-9/
  DATA bm1cs(8)/ + .2543393338072582442742484397174D-10/
  DATA bm1cs(9)/ - .2941045772822967523489750827909D-11/
  DATA bm1cs(10)/ + .3743392025493903309265056153626D-12/
  DATA bm1cs(11)/ - .5149118293821167218720548243527D-13/
  DATA bm1cs(12)/ + .7552535949865143908034040764199D-14/
  DATA bm1cs(13)/ - .1169409706828846444166290622464D-14/
  DATA bm1cs(14)/ + .1896562449434791571721824605060D-15/
  DATA bm1cs(15)/ - .3201955368693286420664775316394D-16/
  DATA bm1cs(16)/ + .5599548399316204114484169905493D-17/
  DATA bm1cs(17)/ - .1010215894730432443119390444544D-17/
  DATA bm1cs(18)/ + .1873844985727562983302042719573D-18/
  DATA bm1cs(19)/ - .3563537470328580219274301439999D-19/
  DATA bm1cs(20)/ + .6931283819971238330422763519999D-20/
  DATA bm1cs(21)/ - .1376059453406500152251408930133D-20/
  DATA bm1cs(22)/ + .2783430784107080220599779327999D-21/
  DATA bm1cs(23)/ - .5727595364320561689348669439999D-22/
  DATA bm1cs(24)/ + .1197361445918892672535756799999D-22/
  DATA bm1cs(25)/ - .2539928509891871976641440426666D-23/
  DATA bm1cs(26)/ + .5461378289657295973069619199999D-24/
  DATA bm1cs(27)/ - .1189211341773320288986289493333D-24/
  DATA bm1cs(28)/ + .2620150977340081594957824000000D-25/
  DATA bm1cs(29)/ - .5836810774255685901920938666666D-26/
  DATA bm1cs(30)/ + .1313743500080595773423615999999D-26/
  DATA bm1cs(31)/ - .2985814622510380355332778666666D-27/
  DATA bm1cs(32)/ + .6848390471334604937625599999999D-28/
  DATA bm1cs(33)/ - .1584401568222476721192960000000D-28/
  DATA bm1cs(34)/ + .3695641006570938054301013333333D-29/
  DATA bm1cs(35)/ - .8687115921144668243012266666666D-30/
  DATA bm1cs(36)/ + .2057080846158763462929066666666D-30/
  DATA bm1cs(37)/ - .4905225761116225518523733333333D-31/
  DATA bt12cs(1)/ + .73823860128742974662620839792764D+0/
  DATA bt12cs(2)/ - .33361113174483906384470147681189D-2/
  DATA bt12cs(3)/ + .61463454888046964698514899420186D-4/
  DATA bt12cs(4)/ - .24024585161602374264977635469568D-5/
  DATA bt12cs(5)/ + .14663555577509746153210591997204D-6/
  DATA bt12cs(6)/ - .11841917305589180567005147504983D-7/
  DATA bt12cs(7)/ + .11574198963919197052125466303055D-8/
  DATA bt12cs(8)/ - .13001161129439187449366007794571D-9/
  DATA bt12cs(9)/ + .16245391141361731937742166273667D-10/
  DATA bt12cs(10)/ - .22089636821403188752155441770128D-11/
  DATA bt12cs(11)/ + .32180304258553177090474358653778D-12/
  DATA bt12cs(12)/ - .49653147932768480785552021135381D-13/
  DATA bt12cs(13)/ + .80438900432847825985558882639317D-14/
  DATA bt12cs(14)/ - .13589121310161291384694712682282D-14/
  DATA bt12cs(15)/ + .23810504397147214869676529605973D-15/
  DATA bt12cs(16)/ - .43081466363849106724471241420799D-16/
  DATA bt12cs(17)/ + .80202544032771002434993512550400D-17/
  DATA bt12cs(18)/ - .15316310642462311864230027468799D-17/
  DATA bt12cs(19)/ + .29928606352715568924073040554666D-18/
  DATA bt12cs(20)/ - .59709964658085443393815636650666D-19/
  DATA bt12cs(21)/ + .12140289669415185024160852650666D-19/
  DATA bt12cs(22)/ - .25115114696612948901006977706666D-20/
  DATA bt12cs(23)/ + .52790567170328744850738380799999D-21/
  DATA bt12cs(24)/ - .11260509227550498324361161386666D-21/
  DATA bt12cs(25)/ + .24348277359576326659663462400000D-22/
  DATA bt12cs(26)/ - .53317261236931800130038442666666D-23/
  DATA bt12cs(27)/ + .11813615059707121039205990399999D-23/
  DATA bt12cs(28)/ - .26465368283353523514856789333333D-24/
  DATA bt12cs(29)/ + .59903394041361503945577813333333D-25/
  DATA bt12cs(30)/ - .13690854630829503109136383999999D-25/
  DATA bt12cs(31)/ + .31576790154380228326413653333333D-26/
  DATA bt12cs(32)/ - .73457915082084356491400533333333D-27/
  DATA bt12cs(33)/ + .17228081480722747930705920000000D-27/
  DATA bt12cs(34)/ - .40716907961286507941068800000000D-28/
  DATA bt12cs(35)/ + .96934745136779622700373333333333D-29/
  DATA bt12cs(36)/ - .23237636337765716765354666666666D-29/
  DATA bt12cs(37)/ + .56074510673522029406890666666666D-30/
  DATA bt12cs(38)/ - .13616465391539005860522666666666D-30/
  DATA bt12cs(39)/ + .33263109233894654388906666666666D-31/
  DATA bm12cs(1)/ + .9807979156233050027272093546937D-1/
  DATA bm12cs(2)/ + .1150961189504685306175483484602D-2/
  DATA bm12cs(3)/ - .4312482164338205409889358097732D-5/
  DATA bm12cs(4)/ + .5951839610088816307813029801832D-7/
  DATA bm12cs(5)/ - .1704844019826909857400701586478D-8/
  DATA bm12cs(6)/ + .7798265413611109508658173827401D-10/
  DATA bm12cs(7)/ - .4958986126766415809491754951865D-11/
  DATA bm12cs(8)/ + .4038432416421141516838202265144D-12/
  DATA bm12cs(9)/ - .3993046163725175445765483846645D-13/
  DATA bm12cs(10)/ + .4619886183118966494313342432775D-14/
  DATA bm12cs(11)/ - .6089208019095383301345472619333D-15/
  DATA bm12cs(12)/ + .8960930916433876482157048041249D-16/
  DATA bm12cs(13)/ - .1449629423942023122916518918925D-16/
  DATA bm12cs(14)/ + .2546463158537776056165149648068D-17/
  DATA bm12cs(15)/ - .4809472874647836444259263718620D-18/
  DATA bm12cs(16)/ + .9687684668292599049087275839124D-19/
  DATA bm12cs(17)/ - .2067213372277966023245038117551D-19/
  DATA bm12cs(18)/ + .4646651559150384731802767809590D-20/
  DATA bm12cs(19)/ - .1094966128848334138241351328339D-20/
  DATA bm12cs(20)/ + .2693892797288682860905707612785D-21/
  DATA bm12cs(21)/ - .6894992910930374477818970026857D-22/
  DATA bm12cs(22)/ + .1830268262752062909890668554740D-22/
  DATA bm12cs(23)/ - .5025064246351916428156113553224D-23/
  DATA bm12cs(24)/ + .1423545194454806039631693634194D-23/
  DATA bm12cs(25)/ - .4152191203616450388068886769801D-24/
  DATA bm12cs(26)/ + .1244609201503979325882330076547D-24/
  DATA bm12cs(27)/ - .3827336370569304299431918661286D-25/
  DATA bm12cs(28)/ + .1205591357815617535374723981835D-25/
  DATA bm12cs(29)/ - .3884536246376488076431859361124D-26/
  DATA bm12cs(30)/ + .1278689528720409721904895283461D-26/
  DATA bm12cs(31)/ - .4295146689447946272061936915912D-27/
  DATA bm12cs(32)/ + .1470689117829070886456802707983D-27/
  DATA bm12cs(33)/ - .5128315665106073128180374017796D-28/
  DATA bm12cs(34)/ + .1819509585471169385481437373286D-28/
  DATA bm12cs(35)/ - .6563031314841980867618635050373D-29/
  DATA bm12cs(36)/ + .2404898976919960653198914875834D-29/
  DATA bm12cs(37)/ - .8945966744690612473234958242979D-30/
  DATA bm12cs(38)/ + .3376085160657231026637148978240D-30/
  DATA bm12cs(39)/ - .1291791454620656360913099916966D-30/
  DATA bm12cs(40)/ + .5008634462958810520684951501254D-31/
  DATA bth1cs(1)/ + .74749957203587276055443483969695D+0/
  DATA bth1cs(2)/ - .12400777144651711252545777541384D-2/
  DATA bth1cs(3)/ + .99252442404424527376641497689592D-5/
  DATA bth1cs(4)/ - .20303690737159711052419375375608D-6/
  DATA bth1cs(5)/ + .75359617705690885712184017583629D-8/
  DATA bth1cs(6)/ - .41661612715343550107630023856228D-9/
  DATA bth1cs(7)/ + .30701618070834890481245102091216D-10/
  DATA bth1cs(8)/ - .28178499637605213992324008883924D-11/
  DATA bth1cs(9)/ + .30790696739040295476028146821647D-12/
  DATA bth1cs(10)/ - .38803300262803434112787347554781D-13/
  DATA bth1cs(11)/ + .55096039608630904934561726208562D-14/
  DATA bth1cs(12)/ - .86590060768383779940103398953994D-15/
  DATA bth1cs(13)/ + .14856049141536749003423689060683D-15/
  DATA bth1cs(14)/ - .27519529815904085805371212125009D-16/
  DATA bth1cs(15)/ + .54550796090481089625036223640923D-17/
  DATA bth1cs(16)/ - .11486534501983642749543631027177D-17/
  DATA bth1cs(17)/ + .25535213377973900223199052533522D-18/
  DATA bth1cs(18)/ - .59621490197413450395768287907849D-19/
  DATA bth1cs(19)/ + .14556622902372718620288302005833D-19/
  DATA bth1cs(20)/ - .37022185422450538201579776019593D-20/
  DATA bth1cs(21)/ + .97763074125345357664168434517924D-21/
  DATA bth1cs(22)/ - .26726821639668488468723775393052D-21/
  DATA bth1cs(23)/ + .75453300384983271794038190655764D-22/
  DATA bth1cs(24)/ - .21947899919802744897892383371647D-22/
  DATA bth1cs(25)/ + .65648394623955262178906999817493D-23/
  DATA bth1cs(26)/ - .20155604298370207570784076869519D-23/
  DATA bth1cs(27)/ + .63417768556776143492144667185670D-24/
  DATA bth1cs(28)/ - .20419277885337895634813769955591D-24/
  DATA bth1cs(29)/ + .67191464220720567486658980018551D-25/
  DATA bth1cs(30)/ - .22569079110207573595709003687336D-25/
  DATA bth1cs(31)/ + .77297719892989706370926959871929D-26/
  DATA bth1cs(32)/ - .26967444512294640913211424080920D-26/
  DATA bth1cs(33)/ + .95749344518502698072295521933627D-27/
  DATA bth1cs(34)/ - .34569168448890113000175680827627D-27/
  DATA bth1cs(35)/ + .12681234817398436504211986238374D-27/
  DATA bth1cs(36)/ - .47232536630722639860464993713445D-28/
  DATA bth1cs(37)/ + .17850008478186376177858619796417D-28/
  DATA bth1cs(38)/ - .68404361004510395406215223566746D-29/
  DATA bth1cs(39)/ + .26566028671720419358293422672212D-29/
  DATA bth1cs(40)/ - .10450402527914452917714161484670D-29/
  DATA bth1cs(41)/ + .41618290825377144306861917197064D-30/
  DATA bth1cs(42)/ - .16771639203643714856501347882887D-30/
  DATA bth1cs(43)/ + .68361997776664389173535928028528D-31/
  DATA bth1cs(44)/ - .28172247861233641166739574622810D-31/
  DATA pi4/0.785398163397448309615660845819876D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9B1MP
  IF ( first ) THEN
    eta = 0.1*REAL(D1MACH(3))
    nbm1 = INITDS(bm1cs,37,eta)
    nbt12 = INITDS(bt12cs,39,eta)
    nbm12 = INITDS(bm12cs,40,eta)
    nbth1 = INITDS(bth1cs,44,eta)
    !
    xmax = 1.0D0/D1MACH(4)
  ENDIF
  first = .FALSE.
  !
  IF ( X<4.0D0 ) THEN
    CALL XERMSG('SLATEC','D9B1MP','X must be .GE. 4',1,2)
    Ampl = 0.0D0
    Theta = 0.0D0
  ELSEIF ( X<=8.0D0 ) THEN
    z = (128.0D0/(X*X)-5.0D0)/3.0D0
    Ampl = (0.75D0+DCSEVL(z,bm1cs,nbm1))/SQRT(X)
    Theta = X - 3.0D0*pi4 + DCSEVL(z,bt12cs,nbt12)/X
  ELSE
    IF ( X>xmax ) CALL XERMSG('SLATEC','D9B1MP',&
      'No precision because X is too big',2,2)
    !
    z = 128.0D0/(X*X) - 1.0D0
    Ampl = (0.75D0+DCSEVL(z,bm12cs,nbm12))/SQRT(X)
    Theta = X - 3.0D0*pi4 + DCSEVL(z,bth1cs,nbth1)/X
  ENDIF
END SUBROUTINE D9B1MP
