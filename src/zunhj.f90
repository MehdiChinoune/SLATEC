!DECK ZUNHJ
SUBROUTINE ZUNHJ(Zr,Zi,Fnu,Ipmtr,Tol,Phir,Phii,Argr,Argi,Zeta1r,Zeta1i,&
    Zeta2r,Zeta2i,Asumr,Asumi,Bsumr,Bsumi)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ZUNHJ
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUNHJ-A, ZUNHJ-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     REFERENCES
  !         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
  !         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
  !
  !         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
  !         PRESS, N.Y., 1974, PAGE 420
  !
  !     ABSTRACT
  !         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
  !         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
  !         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
  !
  !         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
  !
  !         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
  !         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
  !
  !               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
  !
  !         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
  !         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
  !
  !         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
  !         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
  !         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
  !
  !***SEE ALSO  ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZLOG, ZSQRT
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
  !***END PROLOGUE  ZUNHJ
  !     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
  !    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
  !    *ZETA2,ZTH
  REAL(8) :: alfa, ang, ap, ar, Argi, Argr, Asumi, Asumr, &
    atol, aw2, azth, beta, br, Bsumi, Bsumr, btol, &
    c, conei, coner, cri, crr, dri, drr, ex1, ex2, &
    Fnu, fn13, fn23, gama, gpi, hpi, Phii, Phir, pi, &
    pp, pr, przthi, przthr, ptfni, ptfnr, raw, raw2, &
    razth, rfnu, rfnu2, rfn13, rtzti, rtztr, rzthi, &
    rzthr, sti, str, sumai, sumar, sumbi, sumbr, &
    test, tfni, tfnr, thpi, Tol, tzai, tzar, t2i, &
    t2r, upi, upr, wi, wr, w2i, w2r, zai, zar, zbi, &
    zbr, zci, zcr, zeroi, zeror, zetai, zetar, &
    Zeta1i, Zeta1r, Zeta2i, Zeta2r, Zi, Zr, zthi, &
    zthr, ZABS, ac, D1MACH
  INTEGER ias, ibs, Ipmtr, is, j, jr, ju, k, kmax, kp1, ks, l, &
    lr, lrp1, l1, l2, m, idum
  DIMENSION ar(14), br(14), c(105), alfa(180), beta(210), gama(30), &
    ap(30), pr(30), pi(30), upr(14), upi(14), crr(14), cri(14)&
    , drr(14), dri(14)
  EXTERNAL ZABS, ZLOG, ZSQRT
  DATA ar(1), ar(2), ar(3), ar(4), ar(5), ar(6), ar(7), ar(8), &
    ar(9), ar(10), ar(11), ar(12), ar(13), ar(14)&
    /1.00000000000000000D+00, 1.04166666666666667D-01, &
    8.35503472222222222D-02, 1.28226574556327160D-01, &
    2.91849026464140464D-01, 8.81627267443757652D-01, &
    3.32140828186276754D+00, 1.49957629868625547D+01, &
    7.89230130115865181D+01, 4.74451538868264323D+02, &
    3.20749009089066193D+03, 2.40865496408740049D+04, &
    1.98923119169509794D+05, 1.79190200777534383D+06/
  DATA br(1), br(2), br(3), br(4), br(5), br(6), br(7), br(8), &
    br(9), br(10), br(11), br(12), br(13), br(14)&
    /1.00000000000000000D+00, -1.45833333333333333D-01, &
    -9.87413194444444444D-02, -1.43312053915895062D-01, &
    -3.17227202678413548D-01, -9.42429147957120249D-01, &
    -3.51120304082635426D+00, -1.57272636203680451D+01, &
    -8.22814390971859444D+01, -4.92355370523670524D+02, &
    -3.31621856854797251D+03, -2.48276742452085896D+04, &
    -2.04526587315129788D+05, -1.83844491706820990D+06/
  DATA c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), &
    c(10), c(11), c(12), c(13), c(14), c(15), c(16), c(17), &
    c(18), c(19), c(20), c(21), c(22), c(23), &
    c(24)/1.00000000000000000D+00, -2.08333333333333333D-01, &
    1.25000000000000000D-01, 3.34201388888888889D-01, &
    -4.01041666666666667D-01, 7.03125000000000000D-02, &
    -1.02581259645061728D+00, 1.84646267361111111D+00, &
    -8.91210937500000000D-01, 7.32421875000000000D-02, &
    4.66958442342624743D+00, -1.12070026162229938D+01, &
    8.78912353515625000D+00, -2.36408691406250000D+00, &
    1.12152099609375000D-01, -2.82120725582002449D+01, &
    8.46362176746007346D+01, -9.18182415432400174D+01, &
    4.25349987453884549D+01, -7.36879435947963170D+00, &
    2.27108001708984375D-01, 2.12570130039217123D+02, &
    -7.65252468141181642D+02, 1.05999045252799988D+03/
  DATA c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32), &
    c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40), &
    c(41), c(42), c(43), c(44), c(45), c(46), c(47), &
    c(48)/ - 6.99579627376132541D+02, 2.18190511744211590D+02, &
    -2.64914304869515555D+01, 5.72501420974731445D-01, &
    -1.91945766231840700D+03, 8.06172218173730938D+03, &
    -1.35865500064341374D+04, 1.16553933368645332D+04, &
    -5.30564697861340311D+03, 1.20090291321635246D+03, &
    -1.08090919788394656D+02, 1.72772750258445740D+00, &
    2.02042913309661486D+04, -9.69805983886375135D+04, &
    1.92547001232531532D+05, -2.03400177280415534D+05, &
    1.22200464983017460D+05, -4.11926549688975513D+04, &
    7.10951430248936372D+03, -4.93915304773088012D+02, &
    6.07404200127348304D+00, -2.42919187900551333D+05, &
    1.31176361466297720D+06, -2.99801591853810675D+06/
  DATA c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56), &
    c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64), &
    c(65), c(66), c(67), c(68), c(69), c(70), c(71), &
    c(72)/3.76327129765640400D+06, -2.81356322658653411D+06, &
    1.26836527332162478D+06, -3.31645172484563578D+05, &
    4.52187689813627263D+04, -2.49983048181120962D+03, &
    2.43805296995560639D+01, 3.28446985307203782D+06, &
    -1.97068191184322269D+07, 5.09526024926646422D+07, &
    -7.41051482115326577D+07, 6.63445122747290267D+07, &
    -3.75671766607633513D+07, 1.32887671664218183D+07, &
    -2.78561812808645469D+06, 3.08186404612662398D+05, &
    -1.38860897537170405D+04, 1.10017140269246738D+02, &
    -4.93292536645099620D+07, 3.25573074185765749D+08, &
    -9.39462359681578403D+08, 1.55359689957058006D+09, &
    -1.62108055210833708D+09, 1.10684281682301447D+09/
  DATA c(73), c(74), c(75), c(76), c(77), c(78), c(79), c(80), &
    c(81), c(82), c(83), c(84), c(85), c(86), c(87), c(88), &
    c(89), c(90), c(91), c(92), c(93), c(94), c(95), &
    c(96)/ - 4.95889784275030309D+08, 1.42062907797533095D+08, &
    -2.44740627257387285D+07, 2.24376817792244943D+06, &
    -8.40054336030240853D+04, 5.51335896122020586D+02, &
    8.14789096118312115D+08, -5.86648149205184723D+09, &
    1.86882075092958249D+10, -3.46320433881587779D+10, &
    4.12801855797539740D+10, -3.30265997498007231D+10, &
    1.79542137311556001D+10, -6.56329379261928433D+09, &
    1.55927986487925751D+09, -2.25105661889415278D+08, &
    1.73951075539781645D+07, -5.49842327572288687D+05, &
    3.03809051092238427D+03, -1.46792612476956167D+10, &
    1.14498237732025810D+11, -3.99096175224466498D+11, &
    8.19218669548577329D+11, -1.09837515608122331D+12/
  DATA c(97), c(98), c(99), c(100), c(101), c(102), c(103), c(104), &
    c(105)/1.00815810686538209D+12, -6.45364869245376503D+11, &
    2.87900649906150589D+11, -8.78670721780232657D+10, &
    1.76347306068349694D+10, -2.16716498322379509D+09, &
    1.43157876718888981D+08, -3.87183344257261262D+06, &
    1.82577554742931747D+04/
  DATA alfa(1), alfa(2), alfa(3), alfa(4), alfa(5), alfa(6), alfa(7), &
    alfa(8), alfa(9), alfa(10), alfa(11), alfa(12), alfa(13), &
    alfa(14), alfa(15), alfa(16), alfa(17), alfa(18), alfa(19), &
    alfa(20), alfa(21), alfa(22)/ - 4.44444444444444444D-03, &
    -9.22077922077922078D-04, -8.84892884892884893D-05, &
    1.65927687832449737D-04, 2.46691372741792910D-04, &
    2.65995589346254780D-04, 2.61824297061500945D-04, &
    2.48730437344655609D-04, 2.32721040083232098D-04, &
    2.16362485712365082D-04, 2.00738858762752355D-04, &
    1.86267636637545172D-04, 1.73060775917876493D-04, &
    1.61091705929015752D-04, 1.50274774160908134D-04, &
    1.40503497391269794D-04, 1.31668816545922806D-04, &
    1.23667445598253261D-04, 1.16405271474737902D-04, &
    1.09798298372713369D-04, 1.03772410422992823D-04, &
    9.82626078369363448D-05/
  DATA alfa(23), alfa(24), alfa(25), alfa(26), alfa(27), alfa(28), &
    alfa(29), alfa(30), alfa(31), alfa(32), alfa(33), alfa(34), &
    alfa(35), alfa(36), alfa(37), alfa(38), alfa(39), alfa(40), &
    alfa(41), alfa(42), alfa(43), alfa(44)/9.32120517249503256D-05, &
    8.85710852478711718D-05, 8.42963105715700223D-05, &
    8.03497548407791151D-05, 7.66981345359207388D-05, &
    7.33122157481777809D-05, 7.01662625163141333D-05, &
    6.72375633790160292D-05, 6.93735541354588974D-04, &
    2.32241745182921654D-04, -1.41986273556691197D-05, &
    -1.16444931672048640D-04, -1.50803558053048762D-04, &
    -1.55121924918096223D-04, -1.46809756646465549D-04, &
    -1.33815503867491367D-04, -1.19744975684254051D-04, &
    -1.06184319207974020D-04, -9.37699549891194492D-05, &
    -8.26923045588193274D-05, -7.29374348155221211D-05, &
    -6.44042357721016283D-05/
  DATA alfa(45), alfa(46), alfa(47), alfa(48), alfa(49), alfa(50), &
    alfa(51), alfa(52), alfa(53), alfa(54), alfa(55), alfa(56), &
    alfa(57), alfa(58), alfa(59), alfa(60), alfa(61), alfa(62), &
    alfa(63), alfa(64), alfa(65), alfa(66)&
    / - 5.69611566009369048D-05, -5.04731044303561628D-05, &
    -4.48134868008882786D-05, -3.98688727717598864D-05, &
    -3.55400532972042498D-05, -3.17414256609022480D-05, &
    -2.83996793904174811D-05, -2.54522720634870566D-05, &
    -2.28459297164724555D-05, -2.05352753106480604D-05, &
    -1.84816217627666085D-05, -1.66519330021393806D-05, &
    -1.50179412980119482D-05, -1.35554031379040526D-05, &
    -1.22434746473858131D-05, -1.10641884811308169D-05, &
    -3.54211971457743841D-04, -1.56161263945159416D-04, &
    3.04465503594936410D-05, 1.30198655773242693D-04, &
    1.67471106699712269D-04, 1.70222587683592569D-04/
  DATA alfa(67), alfa(68), alfa(69), alfa(70), alfa(71), alfa(72), &
    alfa(73), alfa(74), alfa(75), alfa(76), alfa(77), alfa(78), &
    alfa(79), alfa(80), alfa(81), alfa(82), alfa(83), alfa(84), &
    alfa(85), alfa(86), alfa(87), alfa(88)/1.56501427608594704D-04, &
    1.36339170977445120D-04, 1.14886692029825128D-04, &
    9.45869093034688111D-05, 7.64498419250898258D-05, &
    6.07570334965197354D-05, 4.74394299290508799D-05, &
    3.62757512005344297D-05, 2.69939714979224901D-05, &
    1.93210938247939253D-05, 1.30056674793963203D-05, &
    7.82620866744496661D-06, 3.59257485819351583D-06, &
    1.44040049814251817D-07, -2.65396769697939116D-06, &
    -4.91346867098485910D-06, -6.72739296091248287D-06, &
    -8.17269379678657923D-06, -9.31304715093561232D-06, &
    -1.02011418798016441D-05, -1.08805962510592880D-05, &
    -1.13875481509603555D-05/
  DATA alfa(89), alfa(90), alfa(91), alfa(92), alfa(93), alfa(94), &
    alfa(95), alfa(96), alfa(97), alfa(98), alfa(99), alfa(100), &
    alfa(101), alfa(102), alfa(103), alfa(104), alfa(105), alfa(106)&
    , alfa(107), alfa(108), alfa(109), alfa(110)&
    / - 1.17519675674556414D-05, -1.19987364870944141D-05, &
    3.78194199201772914D-04, 2.02471952761816167D-04, &
    -6.37938506318862408D-05, -2.38598230603005903D-04, &
    -3.10916256027361568D-04, -3.13680115247576316D-04, &
    -2.78950273791323387D-04, -2.28564082619141374D-04, &
    -1.75245280340846749D-04, -1.25544063060690348D-04, &
    -8.22982872820208365D-05, -4.62860730588116458D-05, &
    -1.72334302366962267D-05, 5.60690482304602267D-06, &
    2.31395443148286800D-05, 3.62642745856793957D-05, &
    4.58006124490188752D-05, 5.24595294959114050D-05, &
    5.68396208545815266D-05, 5.94349820393104052D-05/
  DATA alfa(111), alfa(112), alfa(113), alfa(114), alfa(115), alfa(116)&
    , alfa(117), alfa(118), alfa(119), alfa(120), alfa(121), &
    alfa(122), alfa(123), alfa(124), alfa(125), alfa(126), alfa(127)&
    , alfa(128), alfa(129), alfa(130)/6.06478527578421742D-05, &
    6.08023907788436497D-05, 6.01577894539460388D-05, &
    5.89199657344698500D-05, 5.72515823777593053D-05, &
    5.52804375585852577D-05, 5.31063773802880170D-05, &
    5.08069302012325706D-05, 4.84418647620094842D-05, &
    4.60568581607475370D-05, -6.91141397288294174D-04, &
    -4.29976633058871912D-04, 1.83067735980039018D-04, &
    6.60088147542014144D-04, 8.75964969951185931D-04, &
    8.77335235958235514D-04, 7.49369585378990637D-04, &
    5.63832329756980918D-04, 3.68059319971443156D-04, &
    1.88464535514455599D-04/
  DATA alfa(131), alfa(132), alfa(133), alfa(134), alfa(135), alfa(136)&
    , alfa(137), alfa(138), alfa(139), alfa(140), alfa(141), &
    alfa(142), alfa(143), alfa(144), alfa(145), alfa(146), alfa(147)&
    , alfa(148), alfa(149), alfa(150)/3.70663057664904149D-05, &
    -8.28520220232137023D-05, -1.72751952869172998D-04, &
    -2.36314873605872983D-04, -2.77966150694906658D-04, &
    -3.02079514155456919D-04, -3.12594712643820127D-04, &
    -3.12872558758067163D-04, -3.05678038466324377D-04, &
    -2.93226470614557331D-04, -2.77255655582934777D-04, &
    -2.59103928467031709D-04, -2.39784014396480342D-04, &
    -2.20048260045422848D-04, -2.00443911094971498D-04, &
    -1.81358692210970687D-04, -1.63057674478657464D-04, &
    -1.45712672175205844D-04, -1.29425421983924587D-04, &
    -1.14245691942445952D-04/
  DATA alfa(151), alfa(152), alfa(153), alfa(154), alfa(155), alfa(156)&
    , alfa(157), alfa(158), alfa(159), alfa(160), alfa(161), &
    alfa(162), alfa(163), alfa(164), alfa(165), alfa(166), alfa(167)&
    , alfa(168), alfa(169), alfa(170)/1.92821964248775885D-03, &
    1.35592576302022234D-03, -7.17858090421302995D-04, &
    -2.58084802575270346D-03, -3.49271130826168475D-03, &
    -3.46986299340960628D-03, -2.82285233351310182D-03, &
    -1.88103076404891354D-03, -8.89531718383947600D-04, &
    3.87912102631035228D-06, 7.28688540119691412D-04, &
    1.26566373053457758D-03, 1.62518158372674427D-03, &
    1.83203153216373172D-03, 1.91588388990527909D-03, &
    1.90588846755546138D-03, 1.82798982421825727D-03, &
    1.70389506421121530D-03, 1.55097127171097686D-03, &
    1.38261421852276159D-03/
  DATA alfa(171), alfa(172), alfa(173), alfa(174), alfa(175), alfa(176)&
    , alfa(177), alfa(178), alfa(179), alfa(180)&
    /1.20881424230064774D-03, 1.03676532638344962D-03, &
    8.71437918068619115D-04, 7.16080155297701002D-04, &
    5.72637002558129372D-04, 4.42089819465802277D-04, &
    3.24724948503090564D-04, 2.20342042730246599D-04, &
    1.28412898401353882D-04, 4.82005924552095464D-05/
  DATA beta(1), beta(2), beta(3), beta(4), beta(5), beta(6), beta(7), &
    beta(8), beta(9), beta(10), beta(11), beta(12), beta(13), &
    beta(14), beta(15), beta(16), beta(17), beta(18), beta(19), &
    beta(20), beta(21), beta(22)/1.79988721413553309D-02, &
    5.59964911064388073D-03, 2.88501402231132779D-03, &
    1.80096606761053941D-03, 1.24753110589199202D-03, &
    9.22878876572938311D-04, 7.14430421727287357D-04, &
    5.71787281789704872D-04, 4.69431007606481533D-04, &
    3.93232835462916638D-04, 3.34818889318297664D-04, &
    2.88952148495751517D-04, 2.52211615549573284D-04, &
    2.22280580798883327D-04, 1.97541838033062524D-04, &
    1.76836855019718004D-04, 1.59316899661821081D-04, &
    1.44347930197333986D-04, 1.31448068119965379D-04, &
    1.20245444949302884D-04, 1.10449144504599392D-04, &
    1.01828770740567258D-04/
  DATA beta(23), beta(24), beta(25), beta(26), beta(27), beta(28), &
    beta(29), beta(30), beta(31), beta(32), beta(33), beta(34), &
    beta(35), beta(36), beta(37), beta(38), beta(39), beta(40), &
    beta(41), beta(42), beta(43), beta(44)/9.41998224204237509D-05, &
    8.74130545753834437D-05, 8.13466262162801467D-05, &
    7.59002269646219339D-05, 7.09906300634153481D-05, &
    6.65482874842468183D-05, 6.25146958969275078D-05, &
    5.88403394426251749D-05, -1.49282953213429172D-03, &
    -8.78204709546389328D-04, -5.02916549572034614D-04, &
    -2.94822138512746025D-04, -1.75463996970782828D-04, &
    -1.04008550460816434D-04, -5.96141953046457895D-05, &
    -3.12038929076098340D-05, -1.26089735980230047D-05, &
    -2.42892608575730389D-07, 8.05996165414273571D-06, &
    1.36507009262147391D-05, 1.73964125472926261D-05, &
    1.98672978842133780D-05/
  DATA beta(45), beta(46), beta(47), beta(48), beta(49), beta(50), &
    beta(51), beta(52), beta(53), beta(54), beta(55), beta(56), &
    beta(57), beta(58), beta(59), beta(60), beta(61), beta(62), &
    beta(63), beta(64), beta(65), beta(66)/2.14463263790822639D-05, &
    2.23954659232456514D-05, 2.28967783814712629D-05, &
    2.30785389811177817D-05, 2.30321976080909144D-05, &
    2.28236073720348722D-05, 2.25005881105292418D-05, &
    2.20981015361991429D-05, 2.16418427448103905D-05, &
    2.11507649256220843D-05, 2.06388749782170737D-05, &
    2.01165241997081666D-05, 1.95913450141179244D-05, &
    1.90689367910436740D-05, 1.85533719641636667D-05, &
    1.80475722259674218D-05, 5.52213076721292790D-04, &
    4.47932581552384646D-04, 2.79520653992020589D-04, &
    1.52468156198446602D-04, 6.93271105657043598D-05, &
    1.76258683069991397D-05/
  DATA beta(67), beta(68), beta(69), beta(70), beta(71), beta(72), &
    beta(73), beta(74), beta(75), beta(76), beta(77), beta(78), &
    beta(79), beta(80), beta(81), beta(82), beta(83), beta(84), &
    beta(85), beta(86), beta(87), beta(88)&
    / - 1.35744996343269136D-05, -3.17972413350427135D-05, &
    -4.18861861696693365D-05, -4.69004889379141029D-05, &
    -4.87665447413787352D-05, -4.87010031186735069D-05, &
    -4.74755620890086638D-05, -4.55813058138628452D-05, &
    -4.33309644511266036D-05, -4.09230193157750364D-05, &
    -3.84822638603221274D-05, -3.60857167535410501D-05, &
    -3.37793306123367417D-05, -3.15888560772109621D-05, &
    -2.95269561750807315D-05, -2.75978914828335759D-05, &
    -2.58006174666883713D-05, -2.41308356761280200D-05, &
    -2.25823509518346033D-05, -2.11479656768912971D-05, &
    -1.98200638885294927D-05, -1.85909870801065077D-05/
  DATA beta(89), beta(90), beta(91), beta(92), beta(93), beta(94), &
    beta(95), beta(96), beta(97), beta(98), beta(99), beta(100), &
    beta(101), beta(102), beta(103), beta(104), beta(105), beta(106)&
    , beta(107), beta(108), beta(109), beta(110)&
    / - 1.74532699844210224D-05, -1.63997823854497997D-05, &
    -4.74617796559959808D-04, -4.77864567147321487D-04, &
    -3.20390228067037603D-04, -1.61105016119962282D-04, &
    -4.25778101285435204D-05, 3.44571294294967503D-05, &
    7.97092684075674924D-05, 1.03138236708272200D-04, &
    1.12466775262204158D-04, 1.13103642108481389D-04, &
    1.08651634848774268D-04, 1.01437951597661973D-04, &
    9.29298396593363896D-05, 8.40293133016089978D-05, &
    7.52727991349134062D-05, 6.69632521975730872D-05, &
    5.92564547323194704D-05, 5.22169308826975567D-05, &
    4.58539485165360646D-05, 4.01445513891486808D-05/
  DATA beta(111), beta(112), beta(113), beta(114), beta(115), beta(116)&
    , beta(117), beta(118), beta(119), beta(120), beta(121), &
    beta(122), beta(123), beta(124), beta(125), beta(126), beta(127)&
    , beta(128), beta(129), beta(130)/3.50481730031328081D-05, &
    3.05157995034346659D-05, 2.64956119950516039D-05, &
    2.29363633690998152D-05, 1.97893056664021636D-05, &
    1.70091984636412623D-05, 1.45547428261524004D-05, &
    1.23886640995878413D-05, 1.04775876076583236D-05, &
    8.79179954978479373D-06, 7.36465810572578444D-04, &
    8.72790805146193976D-04, 6.22614862573135066D-04, &
    2.85998154194304147D-04, 3.84737672879366102D-06, &
    -1.87906003636971558D-04, -2.97603646594554535D-04, &
    -3.45998126832656348D-04, -3.53382470916037712D-04, &
    -3.35715635775048757D-04/
  DATA beta(131), beta(132), beta(133), beta(134), beta(135), beta(136)&
    , beta(137), beta(138), beta(139), beta(140), beta(141), &
    beta(142), beta(143), beta(144), beta(145), beta(146), beta(147)&
    , beta(148), beta(149), beta(150)/ - 3.04321124789039809D-04, &
    -2.66722723047612821D-04, -2.27654214122819527D-04, &
    -1.89922611854562356D-04, -1.55058918599093870D-04, &
    -1.23778240761873630D-04, -9.62926147717644187D-05, &
    -7.25178327714425337D-05, -5.22070028895633801D-05, &
    -3.50347750511900522D-05, -2.06489761035551757D-05, &
    -8.70106096849767054D-06, 1.13698686675100290D-06, &
    9.16426474122778849D-06, 1.56477785428872620D-05, &
    2.08223629482466847D-05, 2.48923381004595156D-05, &
    2.80340509574146325D-05, 3.03987774629861915D-05, &
    3.21156731406700616D-05/
  DATA beta(151), beta(152), beta(153), beta(154), beta(155), beta(156)&
    , beta(157), beta(158), beta(159), beta(160), beta(161), &
    beta(162), beta(163), beta(164), beta(165), beta(166), beta(167)&
    , beta(168), beta(169), beta(170)/ - 1.80182191963885708D-03, &
    -2.43402962938042533D-03, -1.83422663549856802D-03, &
    -7.62204596354009765D-04, 2.39079475256927218D-04, &
    9.49266117176881141D-04, 1.34467449701540359D-03, &
    1.48457495259449178D-03, 1.44732339830617591D-03, &
    1.30268261285657186D-03, 1.10351597375642682D-03, &
    8.86047440419791759D-04, 6.73073208165665473D-04, &
    4.77603872856582378D-04, 3.05991926358789362D-04, &
    1.60315694594721630D-04, 4.00749555270613286D-05, &
    -5.66607461635251611D-05, -1.32506186772982638D-04, &
    -1.90296187989614057D-04/
  DATA beta(171), beta(172), beta(173), beta(174), beta(175), beta(176)&
    , beta(177), beta(178), beta(179), beta(180), beta(181), &
    beta(182), beta(183), beta(184), beta(185), beta(186), beta(187)&
    , beta(188), beta(189), beta(190)/ - 2.32811450376937408D-04, &
    -2.62628811464668841D-04, -2.82050469867598672D-04, &
    -2.93081563192861167D-04, -2.97435962176316616D-04, &
    -2.96557334239348078D-04, -2.91647363312090861D-04, &
    -2.83696203837734166D-04, -2.73512317095673346D-04, &
    -2.61750155806768580D-04, 6.38585891212050914D-03, &
    9.62374215806377941D-03, 7.61878061207001043D-03, &
    2.83219055545628054D-03, -2.09841352012720090D-03, &
    -5.73826764216626498D-03, -7.70804244495414620D-03, &
    -8.21011692264844401D-03, -7.65824520346905413D-03, &
    -6.47209729391045177D-03/
  DATA beta(191), beta(192), beta(193), beta(194), beta(195), beta(196)&
    , beta(197), beta(198), beta(199), beta(200), beta(201), &
    beta(202), beta(203), beta(204), beta(205), beta(206), beta(207)&
    , beta(208), beta(209), beta(210)/ - 4.99132412004966473D-03, &
    -3.45612289713133280D-03, -2.01785580014170775D-03, &
    -7.59430686781961401D-04, 2.84173631523859138D-04, &
    1.10891667586337403D-03, 1.72901493872728771D-03, &
    2.16812590802684701D-03, 2.45357710494539735D-03, &
    2.61281821058334862D-03, 2.67141039656276912D-03, &
    2.65203073395980430D-03, 2.57411652877287315D-03, &
    2.45389126236094427D-03, 2.30460058071795494D-03, &
    2.13684837686712662D-03, 1.95896528478870911D-03, &
    1.77737008679454412D-03, 1.59690280765839059D-03, &
    1.42111975664438546D-03/
  DATA gama(1), gama(2), gama(3), gama(4), gama(5), gama(6), gama(7), &
    gama(8), gama(9), gama(10), gama(11), gama(12), gama(13), &
    gama(14), gama(15), gama(16), gama(17), gama(18), gama(19), &
    gama(20), gama(21), gama(22)/6.29960524947436582D-01, &
    2.51984209978974633D-01, 1.54790300415655846D-01, &
    1.10713062416159013D-01, 8.57309395527394825D-02, &
    6.97161316958684292D-02, 5.86085671893713576D-02, &
    5.04698873536310685D-02, 4.42600580689154809D-02, &
    3.93720661543509966D-02, 3.54283195924455368D-02, &
    3.21818857502098231D-02, 2.94646240791157679D-02, &
    2.71581677112934479D-02, 2.51768272973861779D-02, &
    2.34570755306078891D-02, 2.19508390134907203D-02, &
    2.06210828235646240D-02, 1.94388240897880846D-02, &
    1.83810633800683158D-02, 1.74293213231963172D-02, &
    1.65685837786612353D-02/
  DATA gama(23), gama(24), gama(25), gama(26), gama(27), gama(28), &
    gama(29), gama(30)/1.57865285987918445D-02, &
    1.50729501494095594D-02, 1.44193250839954639D-02, &
    1.38184805735341786D-02, 1.32643378994276568D-02, &
    1.27517121970498651D-02, 1.22761545318762767D-02, &
    1.18338262398482403D-02/
  DATA ex1, ex2, hpi, gpi, thpi/3.33333333333333333D-01, &
    6.66666666666666667D-01, 1.57079632679489662D+00, &
    3.14159265358979324D+00, 4.71238898038468986D+00/
  DATA zeror, zeroi, coner, conei/0.0D0, 0.0D0, 1.0D0, 0.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZUNHJ
  rfnu = 1.0D0/Fnu
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST (Z/FNU TOO SMALL)
  !-----------------------------------------------------------------------
  test = D1MACH(1)*1.0D+3
  ac = Fnu*test
  IF ( ABS(Zr)>ac.OR.ABS(Zi)>ac ) THEN
    zbr = Zr*rfnu
    zbi = Zi*rfnu
    rfnu2 = rfnu*rfnu
    !-----------------------------------------------------------------------
    !     COMPUTE IN THE FOURTH QUADRANT
    !-----------------------------------------------------------------------
    fn13 = Fnu**ex1
    fn23 = fn13*fn13
    rfn13 = 1.0D0/fn13
    w2r = coner - zbr*zbr + zbi*zbi
    w2i = conei - zbr*zbi - zbr*zbi
    aw2 = ZABS(w2r,w2i)
    IF ( aw2>0.25D0 ) THEN
      !-----------------------------------------------------------------------
      !     ABS(W2).GT.0.25D0
      !-----------------------------------------------------------------------
      CALL ZSQRT(w2r,w2i,wr,wi)
      IF ( wr<0.0D0 ) wr = 0.0D0
      IF ( wi<0.0D0 ) wi = 0.0D0
      str = coner + wr
      sti = wi
      CALL ZDIV(str,sti,zbr,zbi,zar,zai)
      CALL ZLOG(zar,zai,zcr,zci,idum)
      IF ( zci<0.0D0 ) zci = 0.0D0
      IF ( zci>hpi ) zci = hpi
      IF ( zcr<0.0D0 ) zcr = 0.0D0
      zthr = (zcr-wr)*1.5D0
      zthi = (zci-wi)*1.5D0
      Zeta1r = zcr*Fnu
      Zeta1i = zci*Fnu
      Zeta2r = wr*Fnu
      Zeta2i = wi*Fnu
      azth = ZABS(zthr,zthi)
      ang = thpi
      IF ( zthr<0.0D0.OR.zthi>=0.0D0 ) THEN
        ang = hpi
        IF ( zthr/=0.0D0 ) THEN
          ang = DATAN(zthi/zthr)
          IF ( zthr<0.0D0 ) ang = ang + gpi
        ENDIF
      ENDIF
      pp = azth**ex2
      ang = ang*ex2
      zetar = pp*COS(ang)
      zetai = pp*SIN(ang)
      IF ( zetai<0.0D0 ) zetai = 0.0D0
      Argr = zetar*fn23
      Argi = zetai*fn23
      CALL ZDIV(zthr,zthi,zetar,zetai,rtztr,rtzti)
      CALL ZDIV(rtztr,rtzti,wr,wi,zar,zai)
      tzar = zar + zar
      tzai = zai + zai
      CALL ZSQRT(tzar,tzai,str,sti)
      Phir = str*rfn13
      Phii = sti*rfn13
      IF ( Ipmtr/=1 ) THEN
        raw = 1.0D0/SQRT(aw2)
        str = wr*raw
        sti = -wi*raw
        tfnr = str*rfnu*raw
        tfni = sti*rfnu*raw
        razth = 1.0D0/azth
        str = zthr*razth
        sti = -zthi*razth
        rzthr = str*razth*rfnu
        rzthi = sti*razth*rfnu
        zcr = rzthr*ar(2)
        zci = rzthi*ar(2)
        raw2 = 1.0D0/aw2
        str = w2r*raw2
        sti = -w2i*raw2
        t2r = str*raw2
        t2i = sti*raw2
        str = t2r*c(2) + c(3)
        sti = t2i*c(2)
        upr(2) = str*tfnr - sti*tfni
        upi(2) = str*tfni + sti*tfnr
        Bsumr = upr(2) + zcr
        Bsumi = upi(2) + zci
        Asumr = zeror
        Asumi = zeroi
        IF ( rfnu>=Tol ) THEN
          przthr = rzthr
          przthi = rzthi
          ptfnr = tfnr
          ptfni = tfni
          upr(1) = coner
          upi(1) = conei
          pp = 1.0D0
          btol = Tol*(ABS(Bsumr)+ABS(Bsumi))
          ks = 0
          kp1 = 2
          l = 3
          ias = 0
          ibs = 0
          DO lr = 2, 12, 2
            lrp1 = lr + 1
            !-----------------------------------------------------------------------
            !     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
            !     NEXT SUMA AND SUMB
            !-----------------------------------------------------------------------
            DO k = lr, lrp1
              ks = ks + 1
              kp1 = kp1 + 1
              l = l + 1
              zar = c(l)
              zai = zeroi
              DO j = 2, kp1
                l = l + 1
                str = zar*t2r - t2i*zai + c(l)
                zai = zar*t2i + zai*t2r
                zar = str
              ENDDO
              str = ptfnr*tfnr - ptfni*tfni
              ptfni = ptfnr*tfni + ptfni*tfnr
              ptfnr = str
              upr(kp1) = ptfnr*zar - ptfni*zai
              upi(kp1) = ptfni*zar + ptfnr*zai
              crr(ks) = przthr*br(ks+1)
              cri(ks) = przthi*br(ks+1)
              str = przthr*rzthr - przthi*rzthi
              przthi = przthr*rzthi + przthi*rzthr
              przthr = str
              drr(ks) = przthr*ar(ks+2)
              dri(ks) = przthi*ar(ks+2)
            ENDDO
            pp = pp*rfnu2
            IF ( ias/=1 ) THEN
              sumar = upr(lrp1)
              sumai = upi(lrp1)
              ju = lrp1
              DO jr = 1, lr
                ju = ju - 1
                sumar = sumar + crr(jr)*upr(ju) - cri(jr)*upi(ju)
                sumai = sumai + crr(jr)*upi(ju) + cri(jr)*upr(ju)
              ENDDO
              Asumr = Asumr + sumar
              Asumi = Asumi + sumai
              test = ABS(sumar) + ABS(sumai)
              IF ( pp<Tol.AND.test<Tol ) ias = 1
            ENDIF
            IF ( ibs/=1 ) THEN
              sumbr = upr(lr+2) + upr(lrp1)*zcr - upi(lrp1)*zci
              sumbi = upi(lr+2) + upr(lrp1)*zci + upi(lrp1)*zcr
              ju = lrp1
              DO jr = 1, lr
                ju = ju - 1
                sumbr = sumbr + drr(jr)*upr(ju) - dri(jr)*upi(ju)
                sumbi = sumbi + drr(jr)*upi(ju) + dri(jr)*upr(ju)
              ENDDO
              Bsumr = Bsumr + sumbr
              Bsumi = Bsumi + sumbi
              test = ABS(sumbr) + ABS(sumbi)
              IF ( pp<btol.AND.test<btol ) ibs = 1
            ENDIF
            IF ( ias==1.AND.ibs==1 ) EXIT
          ENDDO
        ENDIF
        Asumr = Asumr + coner
        str = -Bsumr*rfn13
        sti = -Bsumi*rfn13
        CALL ZDIV(str,sti,rtztr,rtzti,Bsumr,Bsumi)
      ENDIF
    ELSE
      !-----------------------------------------------------------------------
      !     POWER SERIES FOR ABS(W2).LE.0.25D0
      !-----------------------------------------------------------------------
      k = 1
      pr(1) = coner
      pi(1) = conei
      sumar = gama(1)
      sumai = zeroi
      ap(1) = 1.0D0
      IF ( aw2>=Tol ) THEN
        DO k = 2, 30
          pr(k) = pr(k-1)*w2r - pi(k-1)*w2i
          pi(k) = pr(k-1)*w2i + pi(k-1)*w2r
          sumar = sumar + pr(k)*gama(k)
          sumai = sumai + pi(k)*gama(k)
          ap(k) = ap(k-1)*aw2
          IF ( ap(k)<Tol ) GOTO 20
        ENDDO
        k = 30
      ENDIF
      20       kmax = k
      zetar = w2r*sumar - w2i*sumai
      zetai = w2r*sumai + w2i*sumar
      Argr = zetar*fn23
      Argi = zetai*fn23
      CALL ZSQRT(sumar,sumai,zar,zai)
      CALL ZSQRT(w2r,w2i,str,sti)
      Zeta2r = str*Fnu
      Zeta2i = sti*Fnu
      str = coner + ex2*(zetar*zar-zetai*zai)
      sti = conei + ex2*(zetar*zai+zetai*zar)
      Zeta1r = str*Zeta2r - sti*Zeta2i
      Zeta1i = str*Zeta2i + sti*Zeta2r
      zar = zar + zar
      zai = zai + zai
      CALL ZSQRT(zar,zai,str,sti)
      Phir = str*rfn13
      Phii = sti*rfn13
      IF ( Ipmtr/=1 ) THEN
        !-----------------------------------------------------------------------
        !     SUM SERIES FOR ASUM AND BSUM
        !-----------------------------------------------------------------------
        sumbr = zeror
        sumbi = zeroi
        DO k = 1, kmax
          sumbr = sumbr + pr(k)*beta(k)
          sumbi = sumbi + pi(k)*beta(k)
        ENDDO
        Asumr = zeror
        Asumi = zeroi
        Bsumr = sumbr
        Bsumi = sumbi
        l1 = 0
        l2 = 30
        btol = Tol*(ABS(Bsumr)+ABS(Bsumi))
        atol = Tol
        pp = 1.0D0
        ias = 0
        ibs = 0
        IF ( rfnu2>=Tol ) THEN
          DO is = 2, 7
            atol = atol/rfnu2
            pp = pp*rfnu2
            IF ( ias/=1 ) THEN
              sumar = zeror
              sumai = zeroi
              DO k = 1, kmax
                m = l1 + k
                sumar = sumar + pr(k)*alfa(m)
                sumai = sumai + pi(k)*alfa(m)
                IF ( ap(k)<atol ) EXIT
              ENDDO
              Asumr = Asumr + sumar*pp
              Asumi = Asumi + sumai*pp
              IF ( pp<Tol ) ias = 1
            ENDIF
            IF ( ibs/=1 ) THEN
              sumbr = zeror
              sumbi = zeroi
              DO k = 1, kmax
                m = l2 + k
                sumbr = sumbr + pr(k)*beta(m)
                sumbi = sumbi + pi(k)*beta(m)
                IF ( ap(k)<atol ) EXIT
              ENDDO
              Bsumr = Bsumr + sumbr*pp
              Bsumi = Bsumi + sumbi*pp
              IF ( pp<btol ) ibs = 1
            ENDIF
            IF ( ias==1.AND.ibs==1 ) EXIT
            l1 = l1 + 30
            l2 = l2 + 30
          ENDDO
        ENDIF
        Asumr = Asumr + coner
        pp = rfnu*rfn13
        Bsumr = Bsumr*pp
        Bsumi = Bsumi*pp
      ENDIF
    ENDIF
  ELSE
    Zeta1r = 2.0D0*ABS(LOG(test)) + Fnu
    Zeta1i = 0.0D0
    Zeta2r = Fnu
    Zeta2i = 0.0D0
    Phir = 1.0D0
    Phii = 0.0D0
    Argr = 1.0D0
    Argi = 0.0D0
    RETURN
  ENDIF
  RETURN
END SUBROUTINE ZUNHJ
