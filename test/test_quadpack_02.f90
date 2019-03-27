MODULE TEST40_MOD
  IMPLICIT NONE

CONTAINS
  !** CDQAG
  SUBROUTINE CDQAG(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAG.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAG-S, CDQAG-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF1G, DF2G, DF3G, DPRIN, DQAG

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      error, exact1, exact2, exact3, pi, result, uflow, work(400)
    INTEGER ier, ip, Ipass, iwork(100), key, Kprint, last, lenw, limit, neval
    DATA pi/0.31415926535897932D+01/
    DATA exact1/0.1154700538379252D+01/
    DATA exact2/0.11780972450996172D+00/
    DATA exact3/0.1855802D+02/
    !* FIRST EXECUTABLE STATEMENT  CDQAG
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 100
    lenw = limit*4
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    key = 6
    epsrel = MAX(SQRT(epmach),0.1D-07)
    a = 0.0D+00
    b = 0.1D+01
    CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    limit = 1
    lenw = limit*4
    b = pi*0.2D+01
    CALL DQAG(DF2G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = D1MACH(1)
    limit = 100
    lenw = limit*4
    CALL DQAG(DF2G,a,b,uflow,0.0D+00,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    b = 0.1D+01
    CALL DQAG(DF3G,a,b,epsabs,epsrel,1,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    lenw = 1
    CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAG FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAG PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAG
  !** CDQAGI
  SUBROUTINE CDQAGI(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAGI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAGI-S, CDQAGI-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DPRIN, DQAGI, DT0, DT1, DT2, DT3, DT4, DT5

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891009  Removed unreferenced variables.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    
    !
    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(4), inf
    REAL(8) :: abserr, bound, D1MACH, epmach, epsabs, epsrel, &
      error, exact0, exact1, exact2, exact3, exact4, oflow, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, Lun, neval
    DATA exact0/2.0D+00/, exact1/0.115470066904D1/
    DATA exact2/0.909864525656D-02/
    DATA exact3/0.31415926535897932D+01/
    DATA exact4/0.19984914554328673D+04/
    !* FIRST EXECUTABLE STATEMENT  CDQAGI
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAGI QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    bound = 0.0D+00
    inf = 1
    CALL DQAGI(DT0,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL DQAGI(DT1,bound,inf,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = D1MACH(1)
    CALL DQAGI(DT2,bound,inf,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    CALL DQAGI(DT3,bound,inf,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 4 OR 3 OR 1 OR 0
    !
    CALL DQAGI(DT4,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 3
    ierv(3) = 1
    ierv(4) = 0
    ip = 0
    IF ( ier==4.OR.ier==3.OR.ier==1.OR.ier==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
    oflow = D1MACH(2)
    CALL DQAGI(DT5,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    CALL DQAGI(DT1,bound,inf,epsabs,0.0D+00,result,abserr,neval,ier,limit,&
      lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGI FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGI PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAGI
  !** CDQAGP
  SUBROUTINE CDQAGP(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAGP.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAGP-S, CDQAGP-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF1P, DF2P, DF3P, DF4P, DPRIN, DQAGP

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(4)
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      error, exact1, exact2, exact3, oflow, points(5), p1, p2, result, uflow, work(405)
    INTEGER ier, ip, Ipass, iwork(205), Kprint, last, leniw, lenw, limit, &
      Lun, neval, npts2
    DATA exact1/0.4285277667368085D+01/
    DATA exact2/0.909864525656D-2/
    DATA exact3/0.31415926535897932D+01/
    DATA p1/0.1428571428571428D+00/
    DATA p2/0.6666666666666667D+00/
    !* FIRST EXECUTABLE STATEMENT  CDQAGP
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAGP QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    npts2 = 4
    limit = 100
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    a = 0.0D+00
    b = 0.1D+01
    points(1) = p1
    points(2) = p2
    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    error = ABS(result-exact1)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    leniw = 10
    lenw = leniw*2 - npts2
    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2, 4, 1 OR 3
    !
    npts2 = 3
    points(1) = 0.1D+00
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    uflow = D1MACH(1)
    a = 0.1D+00
    CALL DQAGP(DF2P,a,b,npts2,points,uflow,0.0D+00,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 3
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1.OR.ier==3 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    npts2 = 2
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    a = 0.0D+00
    b = 0.5D+01
    CALL DQAGP(DF3P,a,b,npts2,points,uflow,0.0D+00,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
    b = 0.1D+01
    CALL DQAGP(DF4P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    oflow = D1MACH(2)
    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    npts2 = 5
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    points(1) = p1
    points(2) = p2
    points(3) = 0.3D+01
    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGP FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGP PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAGP
  !** CDQAGS
  SUBROUTINE CDQAGS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAGS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAGS-S, CDQAGS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF0S, DF1S, DF2S, DF3S, DF4S, DF5S, DPRIN,
    !                    DQAGS

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   911114  Modified test on IER=4 to allow IER=5.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(5), Lun
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      error, exact0, exact1, exact2, exact3, exact4, oflow, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    DATA exact0/0.2D+01/
    DATA exact1/0.115470066904D+01/
    DATA exact2/0.909864525656D-02/
    DATA exact3/0.31415926535897932D+01/
    DATA exact4/0.19984914554328673D+04/
    !* FIRST EXECUTABLE STATEMENT  CDQAGS
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAGS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    a = 0.0D+00
    b = 0.1D+01
    CALL DQAGS(DF0S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL DQAGS(DF1S,a,b,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = D1MACH(1)
    a = 0.1D+00
    CALL DQAGS(DF2S,a,b,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    a = 0.0D+00
    b = 0.5D+01
    CALL DQAGS(DF3S,a,b,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 4, OR 5 OR 3 OR 1 OR 0
    !
    b = 0.1D+01
    CALL DQAGS(DF4S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 5
    ierv(3) = 3
    ierv(4) = 1
    ierv(5) = 0
    ip = 0
    IF ( ier==5.OR.ier==4.OR.ier==3.OR.ier==1.OR.ier==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,5)
    !
    ! TEST ON IER = 5
    !
    oflow = D1MACH(2)
    CALL DQAGS(DF5S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    CALL DQAGS(DF1S,a,b,epsabs,0.0D+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGS FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGS PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAGS
  !** CDQAWC
  SUBROUTINE CDQAWC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAWC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAWC-S, CDQAWC-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF0C, DF1C, DPRIN, DQAWC

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      error, exact0, exact1, c, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    DATA exact0/-0.6284617285065624D+03/
    DATA exact1/0.1855802D+01/
    !* FIRST EXECUTABLE STATEMENT  CDQAWC
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWC QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    c = 0.5D+00
    a = -1.0D+00
    b = 1.0D+00
    limit = 200
    lenw = limit*4
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = D1MACH(1)
    CALL DQAWC(DF0C,a,b,c,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    CALL DQAWC(DF1C,0.0D+00,b,c,uflow,0.0D+00,result,abserr,neval,ier,limit,&
      lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    epsabs = 0.0D+00
    epsrel = 0.0D+00
    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWC FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWC PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAWC
  !** CDQAWF
  SUBROUTINE CDQAWF(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAWF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAWF-S, CDQAWF-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF0F, DF1F, DPRIN, DQAWF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(4), integr, iwork(450), leniw, Lun, maxp1
    REAL(8) :: a, abserr, D1MACH, epsabs, epmach, error, exact0, &
      omega, pi, result, uflow, work(1425)
    INTEGER ier, ip, Ipass, Kprint, lenw, limit, limlst, lst, neval
    DATA exact0/0.1422552162575912D+01/
    DATA pi/0.31415926535897932D+01/
    !* FIRST EXECUTABLE STATEMENT  CDQAWF
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWF QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    limlst = 50
    limit = 200
    leniw = limit*2 + limlst
    lenw = leniw*2 + maxp1*25
    epmach = D1MACH(4)
    epsabs = MAX(SQRT(epmach),0.1D-02)
    a = 0.0D+00
    omega = 0.8D+01
    integr = 2
    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsabs ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    limlst = 3
    leniw = 403
    lenw = leniw*2 + maxp1*25
    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    limlst = 50
    leniw = limit*2 + limlst
    lenw = leniw*2 + maxp1*25
    uflow = D1MACH(1)
    CALL DQAWF(DF1F,a,0.0D+00,1,uflow,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 6
    !
    limlst = 50
    leniw = 20
    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 7
    !
    limlst = 50
    leniw = 52
    lenw = leniw*2 + maxp1*25
    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==7 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,7,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWF FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWF PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAWF
  !** CDQAWO
  SUBROUTINE CDQAWO(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAWO.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAWO-S, CDQAWO-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF0O, DF1O, DF2O, DPRIN, DQAWO

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER leniw
    REAL(8) :: a, abserr, b, epmach, epsabs, epsrel, error, &
      exact0, oflow, omega, pi, result, D1MACH, uflow, work(1325)
    INTEGER ier, ierv(4), integr, ip, Ipass, iwork(400), Kprint, last, lenw, &
      Lun, maxp1, neval
    DATA exact0/0.1042872789432789D+05/
    DATA pi/0.31415926535897932D+01/
    !* FIRST EXECUTABLE STATEMENT  CDQAWO
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWO QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    leniw = 400
    lenw = leniw*2 + maxp1*25
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    a = 0.0D+00
    b = pi
    omega = 0.1D+01
    integr = 2
    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    leniw = 2
    lenw = leniw*2 + maxp1*25
    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = D1MACH(1)
    leniw = 400
    lenw = leniw*2 + maxp1*25
    CALL DQAWO(DF0O,a,b,omega,integr,uflow,0.0D+00,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    b = 0.5D+01
    omega = 0.0D+00
    integr = 1
    CALL DQAWO(DF1O,a,b,omega,integr,uflow,0.0D+00,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
    b = 0.1D+01
    oflow = D1MACH(2)
    CALL DQAWO(DF2O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    integr = 3
    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWO FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWO PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAWO
  !** CDQAWS
  SUBROUTINE CDQAWS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQAWS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQAWS-S, CDQAWS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF0WS, DF1WS, DPRIN, DQAWS

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      error, exact0, exact1, alfa, beta, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval, integr
    DATA exact0/0.5350190569223644D+00/
    DATA exact1/0.1998491554328673D+04/
    !* FIRST EXECUTABLE STATEMENT  CDQAWS
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    alfa = -0.5D+00
    beta = -0.5D+00
    integr = 1
    a = 0.0D+00
    b = 0.1D+01
    limit = 200
    lenw = limit*4
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
      ier,limit,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
      ier,2,8,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = D1MACH(1)
    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,uflow,0.0D+00,result,abserr,neval,&
      ier,limit,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
      ier,limit,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    integr = 0
    CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
      ier,limit,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWS FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWS PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQAWS
  !** CDQNG
  SUBROUTINE CDQNG(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQNG.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (CQNG-S, CDQNG-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DF1N, DF2N, DPRIN, DQNG

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER Lun
    REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
      exact1, error, exact2, result, uflow
    INTEGER ier, ierv(1), ip, Ipass, Kprint, neval
    DATA exact1/0.7281029132255818D+00/
    DATA exact2/0.1D+02/
    !* FIRST EXECUTABLE STATEMENT  CDQNG
    IF ( Kprint>=2 ) WRITE (Lun,'(''1DQNG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    epsabs = 0.0D+00
    epmach = D1MACH(4)
    uflow = D1MACH(1)
    epsrel = MAX(SQRT(epmach),0.1D-07)
    a = 0.0D+00
    b = 0.1D+01
    CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
    CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL DQNG(DF2N,a,b,uflow,0.0D+00,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    epsabs = 0.0D+00
    epsrel = 0.0D+00
    CALL DQNG(DF1N,a,b,epsabs,0.0D+00,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQNG FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQNG PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CDQNG
  !** DF0C
  REAL(8) FUNCTION DF0C(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0C
    DF0C = 1.D0/(X*X+1.D-4)
  END FUNCTION DF0C
  !** DF0F
  REAL(8) FUNCTION DF0F(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0F
    DF0F = 0.0D+00
    IF ( X/=0.0D+00 ) DF0F = SIN(0.5D+02*X)/(X*SQRT(X))
  END FUNCTION DF0F
  !** DF0O
  REAL(8) FUNCTION DF0O(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0O
    DF0O = (0.2D+01*SIN(X))**14
  END FUNCTION DF0O
  !** DF0S
  REAL(8) FUNCTION DF0S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0S
    DF0S = 0.0D+00
    IF ( X/=0.0D+00 ) DF0S = 0.1D+01/SQRT(X)
  END FUNCTION DF0S
  !** DF0WS
  REAL(8) FUNCTION DF0WS(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0WS
    DF0WS = SIN(0.1D+02*X)
  END FUNCTION DF0WS
  !** DF1C
  REAL(8) FUNCTION DF1C(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1C
    DF1C = 0.0D+00
    IF ( X/=0.33D+00 ) DF1C = (X-0.5D+00)*ABS(X-0.33D+00)**(-0.9D+00)
  END FUNCTION DF1C
  !** DF1F
  REAL(8) FUNCTION DF1F(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DF1F
    x1 = X + 0.1D+01
    DF1F = 0.5D+01/x1/x1
    y = 0.5D+01/x1
    IF ( y>3.1415926535897932D0 ) DF1F = 0.0D0
  END FUNCTION DF1F
  !** DF1G
  REAL(8) FUNCTION DF1G(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: pi, X
    DATA pi/3.1415926535897932D0/
    !* FIRST EXECUTABLE STATEMENT  DF1G
    DF1G = 2.0D0/(2.0D0+SIN(10.0D0*pi*X))
  END FUNCTION DF1G
  !** DF1N
  REAL(8) FUNCTION DF1N(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1N
    DF1N = 1.0D0/(X**4+X**2+1.0D0)
  END FUNCTION DF1N
  !** DF1O
  REAL(8) FUNCTION DF1O(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1O
    DF1O = 0.1D+01
    IF ( X>0.31415926535897932D+01 ) DF1O = 0.0D+00
  END FUNCTION DF1O
  !** DF1P
  REAL(8) FUNCTION DF1P(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: alfa1, alfa2, p1, p2, X, d1, d2
    !* FIRST EXECUTABLE STATEMENT  DF1P
    !  P1 = 1/7, P2 = 2/3
    DATA p1/0.1428571428571428D+00/
    DATA p2/0.6666666666666667D+00/
    alfa1 = -0.25D0
    alfa2 = -0.5D0
    d1 = ABS(X-p1)
    d2 = ABS(X-p2)
    DF1P = 0.0D+00
    IF ( d1/=0.0D+00.AND.d2/=0.0D+00 ) DF1P = d1**alfa1 + d2**alfa2
  END FUNCTION DF1P
  !** DF1S
  REAL(8) FUNCTION DF1S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1S
    DF1S = 0.2D+01/(0.2D+01+SIN(0.314159D+02*X))
  END FUNCTION DF1S
  !** DF1WS
  REAL(8) FUNCTION DF1WS(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1WS
    DF1WS = 0.00D+00
    IF ( X-0.33D+00/=0.00D+00 ) DF1WS = ABS(X-0.33D+00)**(-0.999D+00)
  END FUNCTION DF1WS
  !** DF2G
  REAL(8) FUNCTION DF2G(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2G
    DF2G = X*SIN(0.3D+02*X)*COS(0.5D+02*X)
  END FUNCTION DF2G
  !** DF2N
  REAL(8) FUNCTION DF2N(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2N
    DF2N = X**(-0.9D+00)
  END FUNCTION DF2N
  !** DF2O
  REAL(8) FUNCTION DF2O(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2O
    DF2O = 0.0D+00
    IF ( X/=0.0D+00 ) DF2O = 0.1D+01/(X*X*SQRT(X))
  END FUNCTION DF2O
  !** DF2P
  REAL(8) FUNCTION DF2P(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2P
    DF2P = SIN(0.314159D+03*X)/(0.314159D+01*X)
  END FUNCTION DF2P
  !** DF2S
  REAL(8) FUNCTION DF2S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2S
    DF2S = 0.1D+03
    IF ( X/=0.0D+00 ) DF2S = SIN(0.314159D+03*X)/(0.314159D+01*X)
  END FUNCTION DF2S
  !** DF3G
  REAL(8) FUNCTION DF3G(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3G
    DF3G = ABS(X-0.33D+00)**(-.90D+00)
  END FUNCTION DF3G
  !** DF3P
  REAL(8) FUNCTION DF3P(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3P
    DF3P = 0.1D+01
    IF ( X>0.31415926535897932D+01 ) DF3P = 0.0D+00
  END FUNCTION DF3P
  !** DF3S
  REAL(8) FUNCTION DF3S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3S
    DF3S = 0.1D+01
    IF ( X>0.31415926535897932D+01 ) DF3S = 0.0D+00
  END FUNCTION DF3S
  !** DF4P
  REAL(8) FUNCTION DF4P(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF4P
    DF4P = 0.0D+00
    IF ( X>0.0D+00 ) DF4P = 0.1D+01/(X*SQRT(X))
  END FUNCTION DF4P
  !** DF4S
  REAL(8) FUNCTION DF4S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF4S
    DF4S = 0.00D+00
    IF ( X-0.33D+00/=0.00D+00 ) DF4S = ABS(X-0.33D+00)**(-0.999D+00)
  END FUNCTION DF4S
  !** DF5S
  REAL(8) FUNCTION DF5S(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: X
    !* FIRST EXECUTABLE STATEMENT  DF5S
    DF5S = 0.0D+00
    IF ( X/=0.0D+00 ) DF5S = 1.0D+00/(X*SQRT(X))
  END FUNCTION DF5S
  !** DT0
  REAL(8) FUNCTION DT0(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF0S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT0
    a = 0.0D+00
    b = 0.1D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT0 = (b-a)*DF0S(y)/x1/x1
  END FUNCTION DT0
  !** DT1
  REAL(8) FUNCTION DT1(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF1S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT1
    a = 0.0D+00
    b = 0.1D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT1 = (b-a)*DF1S(y)/x1/x1
  END FUNCTION DT1
  !** DT2
  REAL(8) FUNCTION DT2(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF2S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT2
    a = 0.1D+00
    b = 0.1D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT2 = (b-a)*DF2S(y)/x1/x1
  END FUNCTION DT2
  !** DT3
  REAL(8) FUNCTION DT3(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF3S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT3
    a = 0.0D+00
    b = 0.5D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT3 = (b-a)*DF3S(y)/x1/x1
  END FUNCTION DT3
  !** DT4
  REAL(8) FUNCTION DT4(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF4S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT4
    a = 0.0D+00
    b = 0.1D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT4 = (b-a)*DF4S(y)/x1/x1
  END FUNCTION DT4
  !** DT5
  REAL(8) FUNCTION DT5(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF5S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL(8) :: a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT5
    a = 0.0D+00
    b = 0.1D+01
    x1 = X + 0.1D+01
    y = (b-a)/x1 + a
    DT5 = (b-a)*DF5S(y)/x1/x1
  END FUNCTION DT5
  !** DPRIN
  SUBROUTINE DPRIN(Lun,Num1,Kprint,Ip,Exact,Result,Abserr,Neval,Ierv,Lierv)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to CDQAG, CDQAG, CDQAGI, CDQAGP, CDQAGS, CDQAWC,
    !            CDQAWF, CDQAWO, CDQAWS, and CDQNG.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  Piessens, Robert
    !             Applied Mathematics and Programming Division
    !             K. U. Leuven
    !           de Doncker, Elise
    !             Applied Mathematics and Programming Division
    !             K. U. Leuven
    !***
    ! **Description:**
    !
    !   This program is called by the (double precision) Quadpack quick
    !   check routines for printing out their messages.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   811027  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910627  Code completely rewritten.  (WRB)
    
    !     .. Scalar Arguments ..
    REAL(8) :: Abserr, Exact, Result
    INTEGER Ip, Kprint, Lierv, Lun, Neval, Num1
    !     .. Array Arguments ..
    INTEGER Ierv(*)
    !     .. Local Scalars ..
    REAL(8) :: error
    INTEGER ier, k
    !     .. Intrinsic Functions ..
    INTRINSIC ABS
    !* FIRST EXECUTABLE STATEMENT  DPRIN
    ier = Ierv(1)
    error = ABS(Exact-Result)
    !
    IF ( Kprint>=2 ) THEN
      IF ( Ip/=1 ) THEN
        !
        !         Write failure messages.
        !
        WRITE (UNIT=Lun,FMT=99002) Num1
        IF ( Num1==0 ) WRITE (UNIT=Lun,FMT=99003)
        IF ( Num1>0 ) WRITE (UNIT=Lun,FMT=99004) Num1
        IF ( Lierv>1 ) WRITE (UNIT=Lun,FMT=99005) (Ierv(k),k=2,Lierv)
        IF ( Num1==6 ) WRITE (UNIT=Lun,FMT=99006)
        WRITE (UNIT=Lun,FMT=99007)
        WRITE (UNIT=Lun,FMT=99008)
        IF ( Num1/=5 ) THEN
          WRITE (UNIT=Lun,FMT=99009) Exact, Result, error, Abserr, ier, Neval
        ELSE
          WRITE (Lun,FMT=99010) Result, Abserr, ier, Neval
        ENDIF
      ELSEIF ( Kprint>=3 ) THEN
        !
        !           Write PASS message.
        !
        WRITE (UNIT=Lun,FMT=99001) Num1
      ENDIF
    ENDIF
    !
    RETURN
    !
    99001 FORMAT (' TEST ON IER = ',I2,' PASSED')
    99002 FORMAT (' TEST ON IER = ',I1,' FAILED.')
    99003 FORMAT (' WE MUST HAVE IER = 0, ERROR.LE.ABSERR AND ABSERR.LE',&
      '.MAX(EPSABS,EPSREL*ABS(EXACT))')
    99004 FORMAT (' WE MUST HAVE IER = ',I1)
    99005 FORMAT (' OR IER =     ',8(I1,2X))
    99006 FORMAT (' RESULT, ABSERR, NEVAL AND EVENTUALLY LAST SHOULD BE',' ZERO')
    99007 FORMAT (' WE HAVE   ')
    99008 FORMAT (7X,'EXACT',11X,'RESULT',6X,'ERROR',4X,'ABSERR',4X,'IER     NEVAL',&
      /,' ',42X,'(EST.ERR.)(FLAG)(NO F-EVAL)')
    99009 FORMAT (' ',2(D15.7,1X),2(D9.2,1X),I4,4X,I6)
    99010 FORMAT (5X,'INFINITY',4X,D15.7,11X,D9.2,I5,4X,I6)
  END SUBROUTINE DPRIN
END MODULE TEST40_MOD
!** TEST40
PROGRAM TEST40
  USE TEST40_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2
  !***
  ! **Type:**      DOUBLE PRECISION (TEST39-S, TEST40-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        DQAG     DQAGI    DQAGP    DQAGS    DQAWC
  !        DQAWF    DQAWO    DQAWS    DQNG
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CDQAG, CDQAGI, CDQAGP, CDQAGS, CDQAWC, CDQAWF,
  !                    CDQAWO, CDQAWS, CDQNG, I1MACH, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST40
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test double precision QUADPACK routines
  !
  !     Test DQAG.
  !
  CALL CDQAG(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGS.
  !
  CALL CDQAGS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGP.
  !
  CALL CDQAGP(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGI.
  !
  CALL CDQAGI(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWO.
  !
  CALL CDQAWO(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWF.
  !
  CALL CDQAWF(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWS.
  !
  CALL CDQAWS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWC.
  !
  CALL CDQAWC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQNG.
  !
  CALL CDQNG(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST40 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST40 *************')
  ENDIF
  STOP
END PROGRAM TEST40
