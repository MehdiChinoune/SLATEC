!** DJAIRY
SUBROUTINE DJAIRY(X,Rx,C,Ai,Dai)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DBESJ and DBESY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (JAIRY-S, DJAIRY-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !           Daniel, S. L., (SNLA)
  !           Weston, M. K., (SNLA)
  !***
  ! **Description:**
  !
  !                  DJAIRY computes the Airy function AI(X)
  !                   and its derivative DAI(X) for DASYJY
  !
  !                                   INPUT
  !
  !         X - Argument, computed by DASYJY, X unrestricted
  !        RX - RX=SQRT(ABS(X)), computed by DASYJY
  !         C - C=2.*(ABS(X)**1.5)/3., computed by DASYJY
  !
  !                                  OUTPUT
  !
  !        AI - Value of function AI(X)
  !       DAI - Value of the derivative DAI(X)
  !
  !***
  ! **See also:**  DBESJ, DBESY
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  
  !
  INTEGER i, j, m1, m1d, m2, m2d, m3, m3d, m4, m4d, n1, n1d, &
    n2, n2d, n3, n3d, n4, n4d
  REAL(8) :: a, Ai, ajn, ajp, ak1, ak2, ak3, b, C, ccv, &
    con2, con3, con4, con5, cv, da, Dai, dajn, dajp, &
    dak1, dak2, dak3, db, ec, e1, e2, fpi12, f1, &
    f2, rtrx, Rx, scv, t, temp1, temp2, tt, X
  DIMENSION ajp(19), ajn(19), a(15), b(15)
  DIMENSION ak1(14), ak2(23), ak3(14)
  DIMENSION dajp(19), dajn(19), da(15), db(15)
  DIMENSION dak1(14), dak2(24), dak3(14)
  SAVE n1, n2, n3, n4, m1, m2, m3, m4, fpi12, con2, con3, con4, &
    con5, ak1, ak2, ak3, ajp, ajn, a, b, n1d, n2d, n3d, n4d, &
    m1d, m2d, m3d, m4d, dak1, dak2, dak3, dajp, dajn, da, db
  DATA n1, n2, n3, n4/14, 23, 19, 15/
  DATA m1, m2, m3, m4/12, 21, 17, 13/
  DATA fpi12, con2, con3, con4, con5/1.30899693899575D+00, &
    5.03154716196777D+00, 3.80004589867293D-01, 8.33333333333333D-01, &
    8.66025403784439D-01/
  DATA ak1(1), ak1(2), ak1(3), ak1(4), ak1(5), ak1(6), ak1(7), ak1(8)&
    , ak1(9), ak1(10), ak1(11), ak1(12), ak1(13), ak1(14)&
    /2.20423090987793D-01, -1.25290242787700D-01, &
    1.03881163359194D-02, 8.22844152006343D-04, -2.34614345891226D-04, &
    1.63824280172116D-05, 3.06902589573189D-07, -1.29621999359332D-07, &
    8.22908158823668D-09, 1.53963968623298D-11, -3.39165465615682D-11, &
    2.03253257423626D-12, -1.10679546097884D-14, -5.16169497785080D-15/
  DATA ak2(1), ak2(2), ak2(3), ak2(4), ak2(5), ak2(6), ak2(7), ak2(8)&
    , ak2(9), ak2(10), ak2(11), ak2(12), ak2(13), ak2(14), ak2(15)&
    , ak2(16), ak2(17), ak2(18), ak2(19), ak2(20), ak2(21), ak2(22)&
    , ak2(23)/2.74366150869598D-01, 5.39790969736903D-03, &
    -1.57339220621190D-03, 4.27427528248750D-04, &
    -1.12124917399925D-04, 2.88763171318904D-05, &
    -7.36804225370554D-06, 1.87290209741024D-06, &
    -4.75892793962291D-07, 1.21130416955909D-07, &
    -3.09245374270614D-08, 7.92454705282654D-09, &
    -2.03902447167914D-09, 5.26863056595742D-10, &
    -1.36704767639569D-10, 3.56141039013708D-11, &
    -9.31388296548430D-12, 2.44464450473635D-12, &
    -6.43840261990955D-13, 1.70106030559349D-13, &
    -4.50760104503281D-14, 1.19774799164811D-14, -3.19077040865066D-15/
  DATA ak3(1), ak3(2), ak3(3), ak3(4), ak3(5), ak3(6), ak3(7), ak3(8)&
    , ak3(9), ak3(10), ak3(11), ak3(12), ak3(13), ak3(14)&
    /2.80271447340791D-01, -1.78127042844379D-03, &
    4.03422579628999D-05, -1.63249965269003D-06, 9.21181482476768D-08, &
    -6.52294330229155D-09, 5.47138404576546D-10, &
    -5.24408251800260D-11, 5.60477904117209D-12, &
    -6.56375244639313D-13, 8.31285761966247D-14, &
    -1.12705134691063D-14, 1.62267976598129D-15, -2.46480324312426D-16/
  DATA ajp(1), ajp(2), ajp(3), ajp(4), ajp(5), ajp(6), ajp(7), ajp(8)&
    , ajp(9), ajp(10), ajp(11), ajp(12), ajp(13), ajp(14), ajp(15)&
    , ajp(16), ajp(17), ajp(18), ajp(19)/7.78952966437581D-02, &
    -1.84356363456801D-01, 3.01412605216174D-02, 3.05342724277608D-02, &
    -4.95424702513079D-03, -1.72749552563952D-03, &
    2.43137637839190D-04, 5.04564777517082D-05, -6.16316582695208D-06, &
    -9.03986745510768D-07, 9.70243778355884D-08, 1.09639453305205D-08, &
    -1.04716330588766D-09, -9.60359441344646D-11, &
    8.25358789454134D-12, 6.36123439018768D-13, -4.96629614116015D-14, &
    -3.29810288929615D-15, 2.35798252031104D-16/
  DATA ajn(1), ajn(2), ajn(3), ajn(4), ajn(5), ajn(6), ajn(7), ajn(8)&
    , ajn(9), ajn(10), ajn(11), ajn(12), ajn(13), ajn(14), ajn(15)&
    , ajn(16), ajn(17), ajn(18), ajn(19)/3.80497887617242D-02, &
    -2.45319541845546D-01, 1.65820623702696D-01, 7.49330045818789D-02, &
    -2.63476288106641D-02, -5.92535597304981D-03, &
    1.44744409589804D-03, 2.18311831322215D-04, -4.10662077680304D-05, &
    -4.66874994171766D-06, 7.15218807277160D-07, 6.52964770854633D-08, &
    -8.44284027565946D-09, -6.44186158976978D-10, &
    7.20802286505285D-11, 4.72465431717846D-12, -4.66022632547045D-13, &
    -2.67762710389189D-14, 2.36161316570019D-15/
  DATA a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), &
    a(10), a(11), a(12), a(13), a(14), a(15)/4.90275424742791D-01, &
    1.57647277946204D-03, -9.66195963140306D-05, 1.35916080268815D-07, &
    2.98157342654859D-07, -1.86824767559979D-08, &
    -1.03685737667141D-09, 3.28660818434328D-10, &
    -2.57091410632780D-11, -2.32357655300677D-12, &
    9.57523279048255D-13, -1.20340828049719D-13, &
    -2.90907716770715D-15, 4.55656454580149D-15, -9.99003874810259D-16/
  DATA b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), &
    b(10), b(11), b(12), b(13), b(14), b(15)/2.78593552803079D-01, &
    -3.52915691882584D-03, -2.31149677384994D-05, &
    4.71317842263560D-06, -1.12415907931333D-07, &
    -2.00100301184339D-08, 2.60948075302193D-09, &
    -3.55098136101216D-11, -3.50849978423875D-11, &
    5.83007187954202D-12, -2.04644828753326D-13, &
    -1.10529179476742D-13, 2.87724778038775D-14, &
    -2.88205111009939D-15, -3.32656311696166D-16/
  DATA n1d, n2d, n3d, n4d/14, 24, 19, 15/
  DATA m1d, m2d, m3d, m4d/12, 22, 17, 13/
  DATA dak1(1), dak1(2), dak1(3), dak1(4), dak1(5), dak1(6), dak1(7), &
    dak1(8), dak1(9), dak1(10), dak1(11), dak1(12), dak1(13), &
    dak1(14)/2.04567842307887D-01, -6.61322739905664D-02, &
    -8.49845800989287D-03, 3.12183491556289D-03, &
    -2.70016489829432D-04, -6.35636298679387D-06, &
    3.02397712409509D-06, -2.18311195330088D-07, &
    -5.36194289332826D-10, 1.13098035622310D-09, &
    -7.43023834629073D-11, 4.28804170826891D-13, 2.23810925754539D-13, &
    -1.39140135641182D-14/
  DATA dak2(1), dak2(2), dak2(3), dak2(4), dak2(5), dak2(6), dak2(7), &
    dak2(8), dak2(9), dak2(10), dak2(11), dak2(12), dak2(13), &
    dak2(14), dak2(15), dak2(16), dak2(17), dak2(18), dak2(19), &
    dak2(20), dak2(21), dak2(22), dak2(23), dak2(24)&
    /2.93332343883230D-01, -8.06196784743112D-03, &
    2.42540172333140D-03, -6.82297548850235D-04, 1.85786427751181D-04, &
    -4.97457447684059D-05, 1.32090681239497D-05, &
    -3.49528240444943D-06, 9.24362451078835D-07, &
    -2.44732671521867D-07, 6.49307837648910D-08, &
    -1.72717621501538D-08, 4.60725763604656D-09, &
    -1.23249055291550D-09, 3.30620409488102D-10, &
    -8.89252099772401D-11, 2.39773319878298D-11, &
    -6.48013921153450D-12, 1.75510132023731D-12, &
    -4.76303829833637D-13, 1.29498241100810D-13, &
    -3.52679622210430D-14, 9.62005151585923D-15, -2.62786914342292D-15/
  DATA dak3(1), dak3(2), dak3(3), dak3(4), dak3(5), dak3(6), dak3(7), &
    dak3(8), dak3(9), dak3(10), dak3(11), dak3(12), dak3(13), &
    dak3(14)/2.84675828811349D-01, 2.53073072619080D-03, &
    -4.83481130337976D-05, 1.84907283946343D-06, &
    -1.01418491178576D-07, 7.05925634457153D-09, &
    -5.85325291400382D-10, 5.56357688831339D-11, &
    -5.90889094779500D-12, 6.88574353784436D-13, &
    -8.68588256452194D-14, 1.17374762617213D-14, &
    -1.68523146510923D-15, 2.55374773097056D-16/
  DATA dajp(1), dajp(2), dajp(3), dajp(4), dajp(5), dajp(6), dajp(7), &
    dajp(8), dajp(9), dajp(10), dajp(11), dajp(12), dajp(13), &
    dajp(14), dajp(15), dajp(16), dajp(17), dajp(18), dajp(19)&
    /6.53219131311457D-02, -1.20262933688823D-01, &
    9.78010236263823D-03, 1.67948429230505D-02, -1.97146140182132D-03, &
    -8.45560295098867D-04, 9.42889620701976D-05, 2.25827860945475D-05, &
    -2.29067870915987D-06, -3.76343991136919D-07, &
    3.45663933559565D-08, 4.29611332003007D-09, -3.58673691214989D-10, &
    -3.57245881361895D-11, 2.72696091066336D-12, 2.26120653095771D-13, &
    -1.58763205238303D-14, -1.12604374485125D-15, 7.31327529515367D-17/
  DATA dajn(1), dajn(2), dajn(3), dajn(4), dajn(5), dajn(6), dajn(7), &
    dajn(8), dajn(9), dajn(10), dajn(11), dajn(12), dajn(13), &
    dajn(14), dajn(15), dajn(16), dajn(17), dajn(18), dajn(19)&
    /1.08594539632967D-02, 8.53313194857091D-02, &
    -3.15277068113058D-01, -8.78420725294257D-02, &
    5.53251906976048D-02, 9.41674060503241D-03, -3.32187026018996D-03, &
    -4.11157343156826D-04, 1.01297326891346D-04, 9.87633682208396D-06, &
    -1.87312969812393D-06, -1.50798500131468D-07, &
    2.32687669525394D-08, 1.59599917419225D-09, -2.07665922668385D-10, &
    -1.24103350500302D-11, 1.39631765331043D-12, 7.39400971155740D-14, &
    -7.32887475627500D-15/
  DATA da(1), da(2), da(3), da(4), da(5), da(6), da(7), da(8), &
    da(9), da(10), da(11), da(12), da(13), da(14), da(15)&
    /4.91627321104601D-01, 3.11164930427489D-03, 8.23140762854081D-05, &
    -4.61769776172142D-06, -6.13158880534626D-08, &
    2.87295804656520D-08, -1.81959715372117D-09, &
    -1.44752826642035D-10, 4.53724043420422D-11, &
    -3.99655065847223D-12, -3.24089119830323D-13, &
    1.62098952568741D-13, -2.40765247974057D-14, 1.69384811284491D-16, &
    8.17900786477396D-16/
  DATA db(1), db(2), db(3), db(4), db(5), db(6), db(7), db(8), &
    db(9), db(10), db(11), db(12), db(13), db(14), db(15)&
    / - 2.77571356944231D-01, 4.44212833419920D-03, &
    -8.42328522190089D-05, -2.58040318418710D-06, &
    3.42389720217621D-07, -6.24286894709776D-09, &
    -2.36377836844577D-09, 3.16991042656673D-10, &
    -4.40995691658191D-12, -5.18674221093575D-12, &
    9.64874015137022D-13, -4.90190576608710D-14, &
    -1.77253430678112D-14, 5.55950610442662D-15, -7.11793337579530D-16/
  !* FIRST EXECUTABLE STATEMENT  DJAIRY
  IF ( X<0.0D0 ) THEN
    !
    IF ( C>5.0D0 ) THEN
      !
      t = 10.0D0/C - 1.0D0
      tt = t + t
      j = n4
      f1 = a(j)
      e1 = b(j)
      f2 = 0.0D0
      e2 = 0.0D0
      DO i = 1, m4
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + a(j)
        e1 = tt*e1 - e2 + b(j)
        f2 = temp1
        e2 = temp2
      ENDDO
      temp1 = t*f1 - f2 + a(1)
      temp2 = t*e1 - e2 + b(1)
      rtrx = SQRT(Rx)
      cv = C - fpi12
      ccv = COS(cv)
      scv = SIN(cv)
      Ai = (temp1*ccv-temp2*scv)/rtrx
      j = n4d
      f1 = da(j)
      e1 = db(j)
      f2 = 0.0D0
      e2 = 0.0D0
      DO i = 1, m4d
        j = j - 1
        temp1 = f1
        temp2 = e1
        f1 = tt*f1 - f2 + da(j)
        e1 = tt*e1 - e2 + db(j)
        f2 = temp1
        e2 = temp2
      ENDDO
      temp1 = t*f1 - f2 + da(1)
      temp2 = t*e1 - e2 + db(1)
      e1 = ccv*con5 + 0.5D0*scv
      e2 = scv*con5 - 0.5D0*ccv
      Dai = (temp1*e1-temp2*e2)*rtrx
      RETURN
    ENDIF
  ELSEIF ( C>5.0D0 ) THEN
    !
    t = 10.0D0/C - 1.0D0
    tt = t + t
    j = n1
    f1 = ak3(j)
    f2 = 0.0D0
    DO i = 1, m1
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + ak3(j)
      f2 = temp1
    ENDDO
    rtrx = SQRT(Rx)
    ec = EXP(-C)
    Ai = ec*(t*f1-f2+ak3(1))/rtrx
    j = n1d
    f1 = dak3(j)
    f2 = 0.0D0
    DO i = 1, m1d
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + dak3(j)
      f2 = temp1
    ENDDO
    Dai = -rtrx*ec*(t*f1-f2+dak3(1))
    RETURN
  ELSEIF ( X>1.20D0 ) THEN
    !
    t = (X+X-con2)*con3
    tt = t + t
    j = n2
    f1 = ak2(j)
    f2 = 0.0D0
    DO i = 1, m2
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + ak2(j)
      f2 = temp1
    ENDDO
    rtrx = SQRT(Rx)
    ec = EXP(-C)
    Ai = ec*(t*f1-f2+ak2(1))/rtrx
    j = n2d
    f1 = dak2(j)
    f2 = 0.0D0
    DO i = 1, m2d
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + dak2(j)
      f2 = temp1
    ENDDO
    Dai = -ec*(t*f1-f2+dak2(1))*rtrx
    RETURN
  ELSE
    t = (X+X-1.2D0)*con4
    tt = t + t
    j = n1
    f1 = ak1(j)
    f2 = 0.0D0
    DO i = 1, m1
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + ak1(j)
      f2 = temp1
    ENDDO
    Ai = t*f1 - f2 + ak1(1)
    !
    j = n1d
    f1 = dak1(j)
    f2 = 0.0D0
    DO i = 1, m1d
      j = j - 1
      temp1 = f1
      f1 = tt*f1 - f2 + dak1(j)
      f2 = temp1
    ENDDO
    Dai = -(t*f1-f2+dak1(1))
    RETURN
  ENDIF
  t = 0.4D0*C - 1.0D0
  tt = t + t
  j = n3
  f1 = ajp(j)
  e1 = ajn(j)
  f2 = 0.0D0
  e2 = 0.0D0
  DO i = 1, m3
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + ajp(j)
    e1 = tt*e1 - e2 + ajn(j)
    f2 = temp1
    e2 = temp2
  ENDDO
  Ai = (t*e1-e2+ajn(1)) - X*(t*f1-f2+ajp(1))
  j = n3d
  f1 = dajp(j)
  e1 = dajn(j)
  f2 = 0.0D0
  e2 = 0.0D0
  DO i = 1, m3d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + dajp(j)
    e1 = tt*e1 - e2 + dajn(j)
    f2 = temp1
    e2 = temp2
  ENDDO
  Dai = X*X*(t*f1-f2+dajp(1)) + (t*e1-e2+dajn(1))
  RETURN
END SUBROUTINE DJAIRY
