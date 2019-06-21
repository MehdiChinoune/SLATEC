!** ZUNIK
SUBROUTINE ZUNIK(Zrr,Zri,Fnu,Ikflg,Ipmtr,Tol,Init,Phir,Phii,Zeta1r,Zeta1i,&
    Zeta2r,Zeta2i,Sumr,Sumi,Cwrkr,Cwrki)
  !> Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUNIK-A, ZUNIK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
  !        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
  !        RESPECTIVELY BY
  !
  !        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
  !
  !        WHERE       ZETA=-ZETA1 + ZETA2       OR
  !                          ZETA1 - ZETA2
  !
  !        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
  !        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
  !        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
  !        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
  !        ZETA1,ZETA2.
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZDIV, ZLOG, ZSQRT

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added EXTERNAL statement with ZLOG and ZSQRT.  (RWC)
  USE service, ONLY : D1MACH
  !     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
  !    *ZETA2,ZN,ZR
  REAL(DP) :: ac, crfni, crfnr, Cwrki(16), Cwrkr(16), Fnu, Phii, Phir, rfn, si, sr, &
    sri, srr, sti, str, Sumi, Sumr, test, ti, Tol, tr, t2i, t2r, Zeta1i, Zeta1r, &
    Zeta2i, Zeta2r, zni, znr, Zri, Zrr
  INTEGER :: i, idum, Ikflg, Init, Ipmtr, j, k, l
  REAL(DP), PARAMETER :: zeror = 0._DP, zeroi = 0._DP, coner = 1._DP, conei = 0._DP
  REAL(DP), PARAMETER :: con(2) = [ 3.98942280401432678E-01_DP, 1.25331413731550025E+00_DP ]
  REAL(DP), PARAMETER :: c(120) = [ 1.00000000000000000E+00_DP, -2.08333333333333333E-01_DP, &
    1.25000000000000000E-01_DP,  3.34201388888888889E-01_DP,-4.01041666666666667E-01_DP, &
    7.03125000000000000E-02_DP, -1.02581259645061728E+00_DP, 1.84646267361111111E+00_DP, &
    -8.91210937500000000E-01_DP, 7.32421875000000000E-02_DP, 4.66958442342624743E+00_DP, &
    -1.12070026162229938E+01_DP, 8.78912353515625000E+00_DP,-2.36408691406250000E+00_DP, &
    1.12152099609375000E-01_DP, -2.82120725582002449E+01_DP, 8.46362176746007346E+01_DP, &
    -9.18182415432400174E+01_DP, 4.25349987453884549E+01_DP,-7.36879435947963170E+00_DP, &
    2.27108001708984375E-01_DP,  2.12570130039217123E+02_DP,-7.65252468141181642E+02_DP, &
    1.05999045252799988E+03_DP, -6.99579627376132541E+02_DP, 2.18190511744211590E+02_DP, &
    -2.64914304869515555E+01_DP, 5.72501420974731445E-01_DP,-1.91945766231840700E+03_DP, &
    8.06172218173730938E+03_DP, -1.35865500064341374E+04_DP, 1.16553933368645332E+04_DP, &
    -5.30564697861340311E+03_DP, 1.20090291321635246E+03_DP,-1.08090919788394656E+02_DP, &
    1.72772750258445740E+00_DP,  2.02042913309661486E+04_DP,-9.69805983886375135E+04_DP, &
    1.92547001232531532E+05_DP, -2.03400177280415534E+05_DP, 1.22200464983017460E+05_DP, &
    -4.11926549688975513E+04_DP, 7.10951430248936372E+03_DP,-4.93915304773088012E+02_DP, &
    6.07404200127348304E+00_DP, -2.42919187900551333E+05_DP, 1.31176361466297720E+06_DP, &
    -2.99801591853810675E+06_DP, 3.76327129765640400E+06_DP,-2.81356322658653411E+06_DP, &
    1.26836527332162478E+06_DP, -3.31645172484563578E+05_DP, 4.52187689813627263E+04_DP, &
    -2.49983048181120962E+03_DP, 2.43805296995560639E+01_DP, 3.28446985307203782E+06_DP, &
    -1.97068191184322269E+07_DP, 5.09526024926646422E+07_DP,-7.41051482115326577E+07_DP, &
    6.63445122747290267E+07_DP, -3.75671766607633513E+07_DP, 1.32887671664218183E+07_DP, &
    -2.78561812808645469E+06_DP, 3.08186404612662398E+05_DP,-1.38860897537170405E+04_DP, &
    1.10017140269246738E+02_DP, -4.93292536645099620E+07_DP, 3.25573074185765749E+08_DP, &
    -9.39462359681578403E+08_DP, 1.55359689957058006E+09_DP,-1.62108055210833708E+09_DP, &
    1.10684281682301447E+09_DP, -4.95889784275030309E+08_DP, 1.42062907797533095E+08_DP, &
    -2.44740627257387285E+07_DP, 2.24376817792244943E+06_DP,-8.40054336030240853E+04_DP, &
    5.51335896122020586E+02_DP,  8.14789096118312115E+08_DP,-5.86648149205184723E+09_DP, &
    1.86882075092958249E+10_DP, -3.46320433881587779E+10_DP, 4.12801855797539740E+10_DP, &
    -3.30265997498007231E+10_DP, 1.79542137311556001E+10_DP,-6.56329379261928433E+09_DP, &
    1.55927986487925751E+09_DP, -2.25105661889415278E+08_DP, 1.73951075539781645E+07_DP, &
    -5.49842327572288687E+05_DP, 3.03809051092238427E+03_DP,-1.46792612476956167E+10_DP, &
    1.14498237732025810E+11_DP, -3.99096175224466498E+11_DP, 8.19218669548577329E+11_DP, &
    -1.09837515608122331E+12_DP, 1.00815810686538209E+12_DP,-6.45364869245376503E+11_DP, &
    2.87900649906150589E+11_DP, -8.78670721780232657E+10_DP, 1.76347306068349694E+10_DP, &
    -2.16716498322379509E+09_DP, 1.43157876718888981E+08_DP,-3.87183344257261262E+06_DP, &
    1.82577554742931747E+04_DP,  2.86464035717679043E+11_DP,-2.40629790002850396E+12_DP, &
    9.10934118523989896E+12_DP, -2.05168994109344374E+13_DP, 3.05651255199353206E+13_DP, &
    -3.16670885847851584E+13_DP, 2.33483640445818409E+13_DP,-1.23204913055982872E+13_DP, &
    4.61272578084913197E+12_DP, -1.19655288019618160E+12_DP, 2.05914503232410016E+11_DP, &
    -2.18229277575292237E+10_DP, 1.24700929351271032E+09_DP,-2.91883881222208134E+07_DP, &
    1.18838426256783253E+05_DP ]
  !* FIRST EXECUTABLE STATEMENT  ZUNIK
  IF( Init==0 ) THEN
    !-----------------------------------------------------------------------
    !     INITIALIZE ALL VARIABLES
    !-----------------------------------------------------------------------
    rfn = 1._DP/Fnu
    !-----------------------------------------------------------------------
    !     OVERFLOW TEST (ZR/FNU TOO SMALL)
    !-----------------------------------------------------------------------
    test = D1MACH(1)*1.E+3_DP
    ac = Fnu*test
    IF( ABS(Zrr)>ac .OR. ABS(Zri)>ac ) THEN
      tr = Zrr*rfn
      ti = Zri*rfn
      sr = coner + (tr*tr-ti*ti)
      si = conei + (tr*ti+ti*tr)
      CALL ZSQRT(sr,si,srr,sri)
      str = coner + srr
      sti = conei + sri
      CALL ZDIV(str,sti,tr,ti,znr,zni)
      CALL ZLOG(znr,zni,str,sti,idum)
      Zeta1r = Fnu*str
      Zeta1i = Fnu*sti
      Zeta2r = Fnu*srr
      Zeta2i = Fnu*sri
      CALL ZDIV(coner,conei,srr,sri,tr,ti)
      srr = tr*rfn
      sri = ti*rfn
      CALL ZSQRT(srr,sri,Cwrkr(16),Cwrki(16))
      Phir = Cwrkr(16)*con(Ikflg)
      Phii = Cwrki(16)*con(Ikflg)
      IF( Ipmtr/=0 ) RETURN
      CALL ZDIV(coner,conei,sr,si,t2r,t2i)
      Cwrkr(1) = coner
      Cwrki(1) = conei
      crfnr = coner
      crfni = conei
      ac = 1._DP
      l = 1
      DO k = 2, 15
        sr = zeror
        si = zeroi
        DO j = 1, k
          l = l + 1
          str = sr*t2r - si*t2i + c(l)
          si = sr*t2i + si*t2r
          sr = str
        END DO
        str = crfnr*srr - crfni*sri
        crfni = crfnr*sri + crfni*srr
        crfnr = str
        Cwrkr(k) = crfnr*sr - crfni*si
        Cwrki(k) = crfnr*si + crfni*sr
        ac = ac*rfn
        test = ABS(Cwrkr(k)) + ABS(Cwrki(k))
        IF( ac<Tol .AND. test<Tol ) GOTO 20
      END DO
      k = 15
      20  Init = k
    ELSE
      Zeta1r = 2._DP*ABS(LOG(test)) + Fnu
      Zeta1i = 0._DP
      Zeta2r = Fnu
      Zeta2i = 0._DP
      Phir = 1._DP
      Phii = 0._DP
      RETURN
    END IF
  END IF
  IF( Ikflg==2 ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE SUM FOR THE K FUNCTION
    !-----------------------------------------------------------------------
    sr = zeror
    si = zeroi
    tr = coner
    DO i = 1, Init
      sr = sr + tr*Cwrkr(i)
      si = si + tr*Cwrki(i)
      tr = -tr
    END DO
    Sumr = sr
    Sumi = si
    Phir = Cwrkr(16)*con(2)
    Phii = Cwrki(16)*con(2)
    RETURN
  END IF
  !-----------------------------------------------------------------------
  !     COMPUTE SUM FOR THE I FUNCTION
  !-----------------------------------------------------------------------
  sr = zeror
  si = zeroi
  DO i = 1, Init
    sr = sr + Cwrkr(i)
    si = si + Cwrki(i)
  END DO
  Sumr = sr
  Sumi = si
  Phir = Cwrkr(16)*con(1)
  Phii = Cwrki(16)*con(1)
  RETURN
END SUBROUTINE ZUNIK
