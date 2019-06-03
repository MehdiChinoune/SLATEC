!** ZUNIK
SUBROUTINE ZUNIK(Zrr,Zri,Fnu,Ikflg,Ipmtr,Tol,Init,Phir,Phii,Zeta1r,Zeta1i,&
    Zeta2r,Zeta2i,Sumr,Sumi,Cwrkr,Cwrki)
  !>
  !  Subsidiary to ZBESI and ZBESK
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
  INTEGER i, idum, Ikflg, Init, Ipmtr, j, k, l
  REAL(DP), PARAMETER :: zeror = 0.0D0, zeroi = 0.0D0, coner = 1.0D0, conei = 0.0D0
  REAL(DP), PARAMETER :: con(2) = [ 3.98942280401432678D-01, 1.25331413731550025D+00 ]
  REAL(DP), PARAMETER :: c(120) = [ 1.00000000000000000D+00, -2.08333333333333333D-01, &
    1.25000000000000000D-01,  3.34201388888888889D-01,-4.01041666666666667D-01, &
    7.03125000000000000D-02, -1.02581259645061728D+00, 1.84646267361111111D+00, &
    -8.91210937500000000D-01, 7.32421875000000000D-02, 4.66958442342624743D+00, &
    -1.12070026162229938D+01, 8.78912353515625000D+00,-2.36408691406250000D+00, &
    1.12152099609375000D-01, -2.82120725582002449D+01, 8.46362176746007346D+01, &
    -9.18182415432400174D+01, 4.25349987453884549D+01,-7.36879435947963170D+00, &
    2.27108001708984375D-01,  2.12570130039217123D+02,-7.65252468141181642D+02, &
    1.05999045252799988D+03, -6.99579627376132541D+02, 2.18190511744211590D+02, &
    -2.64914304869515555D+01, 5.72501420974731445D-01,-1.91945766231840700D+03, &
    8.06172218173730938D+03, -1.35865500064341374D+04, 1.16553933368645332D+04, &
    -5.30564697861340311D+03, 1.20090291321635246D+03,-1.08090919788394656D+02, &
    1.72772750258445740D+00,  2.02042913309661486D+04,-9.69805983886375135D+04, &
    1.92547001232531532D+05, -2.03400177280415534D+05, 1.22200464983017460D+05, &
    -4.11926549688975513D+04, 7.10951430248936372D+03,-4.93915304773088012D+02, &
    6.07404200127348304D+00, -2.42919187900551333D+05, 1.31176361466297720D+06, &
    -2.99801591853810675D+06, 3.76327129765640400D+06,-2.81356322658653411D+06, &
    1.26836527332162478D+06, -3.31645172484563578D+05, 4.52187689813627263D+04, &
    -2.49983048181120962D+03, 2.43805296995560639D+01, 3.28446985307203782D+06, &
    -1.97068191184322269D+07, 5.09526024926646422D+07,-7.41051482115326577D+07, &
    6.63445122747290267D+07, -3.75671766607633513D+07, 1.32887671664218183D+07, &
    -2.78561812808645469D+06, 3.08186404612662398D+05,-1.38860897537170405D+04, &
    1.10017140269246738D+02, -4.93292536645099620D+07, 3.25573074185765749D+08, &
    -9.39462359681578403D+08, 1.55359689957058006D+09,-1.62108055210833708D+09, &
    1.10684281682301447D+09, -4.95889784275030309D+08, 1.42062907797533095D+08, &
    -2.44740627257387285D+07, 2.24376817792244943D+06,-8.40054336030240853D+04, &
    5.51335896122020586D+02,  8.14789096118312115D+08,-5.86648149205184723D+09, &
    1.86882075092958249D+10, -3.46320433881587779D+10, 4.12801855797539740D+10, &
    -3.30265997498007231D+10, 1.79542137311556001D+10,-6.56329379261928433D+09, &
    1.55927986487925751D+09, -2.25105661889415278D+08, 1.73951075539781645D+07, &
    -5.49842327572288687D+05, 3.03809051092238427D+03,-1.46792612476956167D+10, &
    1.14498237732025810D+11, -3.99096175224466498D+11, 8.19218669548577329D+11, &
    -1.09837515608122331D+12, 1.00815810686538209D+12,-6.45364869245376503D+11, &
    2.87900649906150589D+11, -8.78670721780232657D+10, 1.76347306068349694D+10, &
    -2.16716498322379509D+09, 1.43157876718888981D+08,-3.87183344257261262D+06, &
    1.82577554742931747D+04,  2.86464035717679043D+11,-2.40629790002850396D+12, &
    9.10934118523989896D+12, -2.05168994109344374D+13, 3.05651255199353206D+13, &
    -3.16670885847851584D+13, 2.33483640445818409D+13,-1.23204913055982872D+13, &
    4.61272578084913197D+12, -1.19655288019618160D+12, 2.05914503232410016D+11, &
    -2.18229277575292237D+10, 1.24700929351271032D+09,-2.91883881222208134D+07, &
    1.18838426256783253D+05 ]
  !* FIRST EXECUTABLE STATEMENT  ZUNIK
  IF ( Init==0 ) THEN
    !-----------------------------------------------------------------------
    !     INITIALIZE ALL VARIABLES
    !-----------------------------------------------------------------------
    rfn = 1.0D0/Fnu
    !-----------------------------------------------------------------------
    !     OVERFLOW TEST (ZR/FNU TOO SMALL)
    !-----------------------------------------------------------------------
    test = D1MACH(1)*1.0D+3
    ac = Fnu*test
    IF ( ABS(Zrr)>ac.OR.ABS(Zri)>ac ) THEN
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
      IF ( Ipmtr/=0 ) RETURN
      CALL ZDIV(coner,conei,sr,si,t2r,t2i)
      Cwrkr(1) = coner
      Cwrki(1) = conei
      crfnr = coner
      crfni = conei
      ac = 1.0D0
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
        IF ( ac<Tol.AND.test<Tol ) GOTO 20
      END DO
      k = 15
      20  Init = k
    ELSE
      Zeta1r = 2.0D0*ABS(LOG(test)) + Fnu
      Zeta1i = 0.0D0
      Zeta2r = Fnu
      Zeta2i = 0.0D0
      Phir = 1.0D0
      Phii = 0.0D0
      RETURN
    END IF
  END IF
  IF ( Ikflg==2 ) THEN
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
