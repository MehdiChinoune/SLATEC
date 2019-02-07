!*==DASYIK.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DASYIK
SUBROUTINE DASYIK(X,Fnu,Kode,Flgik,Ra,Arg,In,Y)
  IMPLICIT NONE
  !*--DASYIK5
  !***BEGIN PROLOGUE  DASYIK
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBESI and DBESK
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (ASYIK-S, DASYIK-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !                    DASYIK computes Bessel functions I and K
  !                  for arguments X.GT.0.0 and orders FNU.GE.35
  !                  on FLGIK = 1 and FLGIK = -1 respectively.
  !
  !                                    INPUT
  !
  !      X    - Argument, X.GT.0.0D0
  !      FNU  - Order of first Bessel function
  !      KODE - A parameter to indicate the scaling option
  !             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
  !                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
  !                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
  !             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
  !                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
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
  !         the I and K Bessel functions for FNU.GE.35 and real
  !         X.GT.0.0D0. The forms are identical except for a change
  !         in sign of some of the terms. This change in sign is
  !         accomplished by means of the FLAG FLGIK = 1 or -1.
  !
  !***SEE ALSO  DBESI, DBESK
  !***ROUTINES CALLED  D1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  !***END PROLOGUE  DASYIK
  !
  INTEGER In , j , jn , k , kk , Kode , l
  DOUBLE PRECISION ak , ap , Arg , c , coef , con , etx , Flgik , fn , Fnu , &
    gln , Ra , s1 , s2 , t , tol , t2 , X , Y , z
  DOUBLE PRECISION D1MACH
  DIMENSION Y(*) , c(65) , con(2)
  SAVE con , c
  DATA con(1) , con(2)/3.98942280401432678D-01 , 1.25331413731550025D+00/
  DATA c(1) , c(2) , c(3) , c(4) , c(5) , c(6) , c(7) , c(8) , c(9) , &
    c(10) , c(11) , c(12) , c(13) , c(14) , c(15) , c(16) , c(17) , &
    c(18) , c(19) , c(20) , c(21) , c(22) , c(23) , &
    c(24)/ - 2.08333333333333D-01 , 1.25000000000000D-01 , &
    3.34201388888889D-01 , -4.01041666666667D-01 , 7.03125000000000D-02 , &
    -1.02581259645062D+00 , 1.84646267361111D+00 , &
    -8.91210937500000D-01 , 7.32421875000000D-02 , 4.66958442342625D+00 , &
    -1.12070026162230D+01 , 8.78912353515625D+00 , &
    -2.36408691406250D+00 , 1.12152099609375D-01 , &
    -2.82120725582002D+01 , 8.46362176746007D+01 , &
    -9.18182415432400D+01 , 4.25349987453885D+01 , &
    -7.36879435947963D+00 , 2.27108001708984D-01 , 2.12570130039217D+02 , &
    -7.65252468141182D+02 , 1.05999045252800D+03 , -6.99579627376133D+02/
  DATA c(25) , c(26) , c(27) , c(28) , c(29) , c(30) , c(31) , c(32) , &
    c(33) , c(34) , c(35) , c(36) , c(37) , c(38) , c(39) , c(40) , &
    c(41) , c(42) , c(43) , c(44) , c(45) , c(46) , c(47) , &
    c(48)/2.18190511744212D+02 , -2.64914304869516D+01 , &
    5.72501420974731D-01 , -1.91945766231841D+03 , 8.06172218173731D+03 , &
    -1.35865500064341D+04 , 1.16553933368645D+04 , &
    -5.30564697861340D+03 , 1.20090291321635D+03 , &
    -1.08090919788395D+02 , 1.72772750258446D+00 , 2.02042913309661D+04 , &
    -9.69805983886375D+04 , 1.92547001232532D+05 , &
    -2.03400177280416D+05 , 1.22200464983017D+05 , &
    -4.11926549688976D+04 , 7.10951430248936D+03 , &
    -4.93915304773088D+02 , 6.07404200127348D+00 , &
    -2.42919187900551D+05 , 1.31176361466298D+06 , &
    -2.99801591853811D+06 , 3.76327129765640D+06/
  DATA c(49) , c(50) , c(51) , c(52) , c(53) , c(54) , c(55) , c(56) , &
    c(57) , c(58) , c(59) , c(60) , c(61) , c(62) , c(63) , c(64) , &
    c(65)/ - 2.81356322658653D+06 , 1.26836527332162D+06 , &
    -3.31645172484564D+05 , 4.52187689813627D+04 , &
    -2.49983048181121D+03 , 2.43805296995561D+01 , 3.28446985307204D+06 , &
    -1.97068191184322D+07 , 5.09526024926646D+07 , &
    -7.41051482115327D+07 , 6.63445122747290D+07 , &
    -3.75671766607634D+07 , 1.32887671664218D+07 , &
    -2.78561812808645D+06 , 3.08186404612662D+05 , &
    -1.38860897537170D+04 , 1.10017140269247D+02/
  !***FIRST EXECUTABLE STATEMENT  DASYIK
  tol = D1MACH(3)
  tol = MAX(tol,1.0D-15)
  fn = Fnu
  z = (3.0D0-Flgik)/2.0D0
  kk = INT(z)
  DO jn = 1 , In
    IF ( jn/=1 ) THEN
      fn = fn - Flgik
      z = X/fn
      Ra = SQRT(1.0D0+z*z)
      gln = LOG((1.0D0+Ra)/z)
      etx = Kode - 1
      t = Ra*(1.0D0-etx) + etx/(z+Ra)
      Arg = fn*(t-gln)*Flgik
    ENDIF
    coef = EXP(Arg)
    t = 1.0D0/Ra
    t2 = t*t
    t = t/fn
    t = SIGN(t,Flgik)
    s2 = 1.0D0
    ap = 1.0D0
    l = 0
    DO k = 2 , 11
      l = l + 1
      s1 = c(l)
      DO j = 2 , k
        l = l + 1
        s1 = s1*t2 + c(l)
      ENDDO
      ap = ap*t
      ak = ap*s1
      s2 = s2 + ak
      IF ( MAX(ABS(ak),ABS(ap))<tol ) EXIT
    ENDDO
    t = ABS(t)
    Y(jn) = s2*coef*SQRT(t)*con(kk)
  ENDDO
END SUBROUTINE DASYIK
