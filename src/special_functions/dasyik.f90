!** DASYIK
PURE SUBROUTINE DASYIK(X,Fnu,Kode,Flgik,Ra,Arg,In,Y)
  !> Subsidiary to DBESI and DBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (ASYIK-S, DASYIK-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !    DASYIK computes Bessel functions I and K for arguments X>0.0 and orders
  !    FNU>=35 on FLGIK = 1 and FLGIK = -1 respectively.
  !
  !                                    INPUT
  !
  !      X    - Argument, X>0.0D0
  !      FNU  - Order of first Bessel function
  !      KODE - A parameter to indicate the scaling option
  !             KODE=1 returns Y(I)=        I_{FNU+I-1}(X), I=1,IN
  !                    or      Y(I)=        K_{FNU+I-1}(X), I=1,IN
  !                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
  !             KODE=2 returns Y(I)=EXP(-X)*I_{FNU+I-1}(X), I=1,IN
  !                    or      Y(I)=EXP( X)*K_{FNU+I-1}(X), I=1,IN
  !                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
  !     FLGIK - Selection parameter for I or K FUNCTION
  !             FLGIK =  1.0D0 gives the I function
  !             FLGIK = -1.0D0 gives the K function
  !        RA - SQRT(1.+Z*Z), Z=X/FNU
  !       ARG - Argument of the leading exponential
  !        IN - Number of functions desired, IN=1 or 2
  !
  !                                    OUTPUT
  !
  !         Y - A vector whose first IN components contain the sequence
  !
  !     Abstract  **** A double precision routine ****
  !         DASYIK implements the uniform asymptotic expansion of
  !         the I and K Bessel functions for FNU>=35 and real
  !         X>0.0D0. The forms are identical except for a change
  !         in sign of some of the terms. This change in sign is
  !         accomplished by means of the FLAG FLGIK = 1 or -1.
  !
  !***
  ! **See also:**  DBESI, DBESK
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  USE service, ONLY : D1MACH
  !
  INTEGER, INTENT(IN) :: In, Kode
  REAL(DP), INTENT(IN) :: Flgik, Fnu, X
  REAL(DP), INTENT(OUT) :: Arg, Ra, Y(In)
  INTEGER :: j, jn, k, kk, l
  REAL(DP) :: ak, ap, coef, etx, fn, gln, s1, s2, t, tol, t2, z
  REAL(DP), PARAMETER :: con(2) = [ 3.98942280401432678E-01_DP, 1.25331413731550025E+00_DP ]
  REAL(DP), PARAMETER :: c(65) = [ -2.08333333333333E-01_DP, 1.25000000000000E-01_DP, &
    3.34201388888889E-01_DP, -4.01041666666667E-01_DP, 7.03125000000000E-02_DP, &
    -1.02581259645062E+00_DP, 1.84646267361111E+00_DP, -8.91210937500000E-01_DP, &
    7.32421875000000E-02_DP, 4.66958442342625E+00_DP, -1.12070026162230E+01_DP, &
    8.78912353515625E+00_DP, -2.36408691406250E+00_DP, 1.12152099609375E-01_DP, &
    -2.82120725582002E+01_DP, 8.46362176746007E+01_DP, -9.18182415432400E+01_DP, &
    4.25349987453885E+01_DP, -7.36879435947963E+00_DP, 2.27108001708984E-01_DP, &
    2.12570130039217E+02_DP, -7.65252468141182E+02_DP, 1.05999045252800E+03_DP, &
    -6.99579627376133E+02_DP, 2.18190511744212E+02_DP, -2.64914304869516E+01_DP, &
    5.72501420974731E-01_DP, -1.91945766231841E+03_DP, 8.06172218173731E+03_DP, &
    -1.35865500064341E+04_DP, 1.16553933368645E+04_DP, -5.30564697861340E+03_DP, &
    1.20090291321635E+03_DP, -1.08090919788395E+02_DP, 1.72772750258446E+00_DP, &
    2.02042913309661E+04_DP, -9.69805983886375E+04_DP, 1.92547001232532E+05_DP, &
    -2.03400177280416E+05_DP, 1.22200464983017E+05_DP, -4.11926549688976E+04_DP, &
    7.10951430248936E+03_DP, -4.93915304773088E+02_DP, 6.07404200127348E+00_DP, &
    -2.42919187900551E+05_DP, 1.31176361466298E+06_DP, -2.99801591853811E+06_DP, &
    3.76327129765640E+06_DP, -2.81356322658653E+06_DP, 1.26836527332162E+06_DP, &
    -3.31645172484564E+05_DP, 4.52187689813627E+04_DP, -2.49983048181121E+03_DP, &
    2.43805296995561E+01_DP, 3.28446985307204E+06_DP, -1.97068191184322E+07_DP, &
    5.09526024926646E+07_DP, -7.41051482115327E+07_DP, 6.63445122747290E+07_DP, &
    -3.75671766607634E+07_DP, 1.32887671664218E+07_DP, -2.78561812808645E+06_DP, &
    3.08186404612662E+05_DP, -1.38860897537170E+04_DP, 1.10017140269247E+02_DP ]
  !* FIRST EXECUTABLE STATEMENT  DASYIK
  tol = D1MACH(3)
  tol = MAX(tol,1.E-15_DP)
  fn = Fnu
  z = (3._DP-Flgik)/2._DP
  kk = INT(z)
  DO jn = 1, In
    IF( jn/=1 ) THEN
      fn = fn - Flgik
      z = X/fn
      Ra = SQRT(1._DP+z*z)
      gln = LOG((1._DP+Ra)/z)
      etx = Kode - 1
      t = Ra*(1._DP-etx) + etx/(z+Ra)
      Arg = fn*(t-gln)*Flgik
    END IF
    coef = EXP(Arg)
    t = 1._DP/Ra
    t2 = t*t
    t = t/fn
    t = SIGN(t,Flgik)
    s2 = 1._DP
    ap = 1._DP
    l = 0
    DO k = 2, 11
      l = l + 1
      s1 = c(l)
      DO j = 2, k
        l = l + 1
        s1 = s1*t2 + c(l)
      END DO
      ap = ap*t
      ak = ap*s1
      s2 = s2 + ak
      IF( MAX(ABS(ak),ABS(ap))<tol ) EXIT
    END DO
    t = ABS(t)
    Y(jn) = s2*coef*SQRT(t)*con(kk)
  END DO

END SUBROUTINE DASYIK