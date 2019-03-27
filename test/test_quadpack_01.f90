MODULE TEST39_MOD
  IMPLICIT NONE

CONTAINS
  !** CQAG
  SUBROUTINE CQAG(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAG.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAG-S, CDQAG-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F1G, F2G, F3G, QAG, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact1, &
      exact2, exact3, pi, result, uflow, work(400)
    INTEGER ier, ip, Ipass, iwork(100), key, Kprint, last, lenw, limit, neval
    DATA pi/0.31415926535897932E+01/
    DATA exact1/0.1154700538379252E+01/
    DATA exact2/0.11780972450996172E+00/
    DATA exact3/0.1855802E+02/
    !* FIRST EXECUTABLE STATEMENT  CQAG
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 100
    lenw = limit*4
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    key = 6
    epsrel = MAX(SQRT(epmach),0.1E-07)
    a = 0.0E+00
    b = 0.1E+01
    CALL QAG(F1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    limit = 1
    lenw = limit*4
    b = pi*0.2E+01
    CALL QAG(F2G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = R1MACH(1)
    limit = 100
    lenw = limit*4
    CALL QAG(F2G,a,b,uflow,0.0E+00,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    b = 0.1E+01
    CALL QAG(F3G,a,b,epsabs,epsrel,1,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    lenw = 1
    CALL QAG(F1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAG FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAG PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAG
  !** CQAGI
  SUBROUTINE CQAGI(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAGI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAGI-S, CDQAGI-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, QAGI, R1MACH, T0, T1, T2, T3, T4, T5

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891009  Removed unreferenced variables.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(4), inf
    REAL abserr, bound, R1MACH, epmach, epsabs, epsrel, error, exact0, &
      exact1, exact2, exact3, exact4, oflow, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, Lun, neval
    DATA exact0/2.0E+00/, exact1/0.115470066904E1/
    DATA exact2/0.909864525656E-02/
    DATA exact3/0.31415926535897932E+01/
    DATA exact4/0.19984914554328673E+04/
    !* FIRST EXECUTABLE STATEMENT  CQAGI
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAGI QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    bound = 0.0E+00
    inf = 1
    CALL QAGI(T0,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL QAGI(T1,bound,inf,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = R1MACH(1)
    CALL QAGI(T2,bound,inf,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1
    !
    CALL QAGI(T3,bound,inf,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 4 OR 3 OR 1
    !
    CALL QAGI(T4,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ierv(2) = 3
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==4.OR.ier==3.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
    oflow = R1MACH(2)
    CALL QAGI(T5,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    CALL QAGI(T1,bound,inf,epsabs,0.0E+00,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAGI FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAGI PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAGI
  !** CQAGP
  SUBROUTINE CQAGP(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAGP.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAGP-S, CDQAGP-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F1P, F2P, F3P, F4P, QAGP, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(4)
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact1, &
      exact2, exact3, oflow, points(5), p1, p2, result, uflow, work(405)
    INTEGER ier, ip, Ipass, iwork(205), Kprint, last, leniw, lenw, limit, &
      Lun, neval, npts2
    DATA exact1/0.4285277667368085E+01/
    DATA exact2/0.909864525656E-2/
    DATA exact3/0.31415926535897932E+01/
    DATA p1/0.1428571428571428E+00/
    DATA p2/0.6666666666666667E+00/
    !* FIRST EXECUTABLE STATEMENT  CQAGP
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAGP QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    npts2 = 4
    limit = 100
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    a = 0.0E+00
    b = 0.1E+01
    points(1) = p1
    points(2) = p2
    CALL QAGP(F1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    error = ABS(result-exact1)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    leniw = 10
    lenw = leniw*2 - npts2
    CALL QAGP(F1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2, 4, 1 OR 3
    !
    npts2 = 3
    points(1) = 0.1E+00
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    uflow = R1MACH(1)
    a = 0.1E+00
    CALL QAGP(F2P,a,b,npts2,points,uflow,0.0E+00,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 3
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1.OR.ier==3 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    npts2 = 2
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    a = 0.0E+00
    b = 0.5E+01
    CALL QAGP(F3P,a,b,npts2,points,uflow,0.0E+00,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
    b = 0.1E+01
    CALL QAGP(F4P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    oflow = R1MACH(2)
    CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    npts2 = 5
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    points(1) = p1
    points(2) = p2
    points(3) = 0.3E+01
    CALL QAGP(F1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAGP FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAGP PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAGP
  !** CQAGS
  SUBROUTINE CQAGS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAGS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAGS-S, CDQAGS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F0S, F1S, F2S, F3S, F4S, F5S, QAGS, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   911114  Modified test on IER=4 to allow IER=5.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(5), Lun
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact0, &
      exact1, exact2, exact3, exact4, oflow, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    DATA exact0/0.2E+01/
    DATA exact1/0.115470066904E+01/
    DATA exact2/0.909864525656E-02/
    DATA exact3/0.31415926535897932E+01/
    DATA exact4/0.19984914554328673E+04/
    !* FIRST EXECUTABLE STATEMENT  CQAGS
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAGS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    a = 0.0E+00
    b = 0.1E+01
    CALL QAGS(F0S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL QAGS(F1S,a,b,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = R1MACH(1)
    a = 0.1E+00
    CALL QAGS(F2S,a,b,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
    a = 0.0E+00
    b = 0.5E+01
    CALL QAGS(F3S,a,b,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ierv(4) = 2
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 4, OR 5 OR 3 OR 1 OR 0
    !
    b = 0.1E+01
    epsrel = 1.E-4
    CALL QAGS(F4S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    !      IER=4
    ierv(1) = ier
    ierv(2) = 5
    ierv(3) = 3
    ierv(4) = 1
    ierv(5) = 0
    ip = 0
    IF ( ier==5.OR.ier==4.OR.ier==3.OR.ier==1.OR.ier==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,5)
    !
    ! TEST ON IER = 5
    !
    oflow = R1MACH(2)
    CALL QAGS(F5S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    CALL QAGS(F1S,a,b,epsabs,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAGS FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAGS PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAGS
  !** CQAWC
  SUBROUTINE CQAWC(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAWC.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAWC-S, CDQAWC-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F0C, F1C, QAWC, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    
    !
    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact0, &
      exact1, c, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    DATA exact0/-0.6284617285065624E+03/
    DATA exact1/0.1855802E+01/
    !* FIRST EXECUTABLE STATEMENT  CQAWC
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWC QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    c = 0.5E+00
    a = -1.0E+00
    b = 1.0E+00
    limit = 200
    lenw = limit*4
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = R1MACH(1)
    CALL QAWC(F0C,a,b,c,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    CALL QAWC(F1C,0.0E+00,b,c,uflow,0.0E+00,result,abserr,neval,ier,limit,&
      lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    epsabs = 0.0E+00
    epsrel = 0.0E+00
    CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAWC FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAWC PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAWC
  !** CQAWF
  SUBROUTINE CQAWF(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAWF.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAWF-S, CDQAWF-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F0F, F1F, QAWF, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    
    !
    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(3), integr, iwork(450), leniw, Lun, maxp1
    REAL a, abserr, R1MACH, epsabs, epmach, error, exact0, &
      omega, pi, result, uflow, work(1425)
    INTEGER ier, ip, Ipass, Kprint, lenw, limit, limlst, lst, neval
    DATA exact0/0.1422552162575912E+01/
    DATA pi/0.31415926535897932E+01/
    !* FIRST EXECUTABLE STATEMENT  CQAWF
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWF QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    limlst = 50
    limit = 200
    leniw = limit*2 + limlst
    lenw = leniw*2 + maxp1*25
    epmach = R1MACH(4)
    epsabs = MAX(SQRT(epmach),0.1E-02)
    a = 0.0E+00
    omega = 0.8E+01
    integr = 2
    CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsabs ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    limlst = 3
    leniw = 403
    lenw = leniw*2 + maxp1*25
    CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 3 OR 4 OR 1
    !
    limlst = 50
    leniw = limit*2 + limlst
    lenw = leniw*2 + maxp1*25
    uflow = R1MACH(1)
    CALL QAWF(F1F,a,0.0E+00,1,uflow,result,abserr,neval,ier,limlst,lst,leniw,&
      maxp1,lenw,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 6
    !
    limlst = 50
    leniw = 20
    CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 7
    !
    limlst = 50
    leniw = 52
    lenw = leniw*2 + maxp1*25
    CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==7 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,7,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAWF FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAWF PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAWF
  !** CQAWO
  SUBROUTINE CQAWO(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAWO.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAWO-S, CDQAWO-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F0O, F1O, F2O, QAWO, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    
    !
    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER leniw
    REAL a, abserr, b, epmach, epsabs, epsrel, error, exact0, &
      oflow, omega, pi, result, R1MACH, uflow, work(1325)
    INTEGER ier, ierv(4), integr, ip, Ipass, iwork(400), Kprint, last, lenw, &
      Lun, maxp1, neval
    DATA exact0/0.1042872789432789E+05/
    DATA pi/0.31415926535897932E+01/
    !* FIRST EXECUTABLE STATEMENT  CQAWO
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWO QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    leniw = 400
    lenw = leniw*2 + maxp1*25
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    a = 0.0E+00
    b = pi
    omega = 0.1E+01
    integr = 2
    CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    leniw = 2
    lenw = leniw*2 + maxp1*25
    CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
    uflow = R1MACH(1)
    leniw = 400
    lenw = leniw*2 + maxp1*25
    CALL QAWO(F0O,a,b,omega,integr,uflow,0.0E+00,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1
    !
    b = 0.5E+01
    omega = 0.0E+00
    integr = 1
    CALL QAWO(F1O,a,b,omega,integr,uflow,0.0E+00,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 4
    ierv(3) = 1
    ip = 0
    IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 5
    !
    b = 0.1E+01
    oflow = R1MACH(2)
    CALL QAWO(F2O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==5 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    integr = 3
    CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
      last==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAWO FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAWO PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAWO
  !** CQAWS
  SUBROUTINE CQAWS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QAWS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQAWS-S, CDQAWS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F0WS, F1WS, QAWS, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER ierv(2), Lun
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact0, &
      exact1, alfa, beta, result, uflow, work(800)
    INTEGER ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval, integr
    DATA exact0/0.5350190569223644E+00/
    DATA exact1/0.1998491554328673E+04/
    !* FIRST EXECUTABLE STATEMENT  CQAWS
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    alfa = -0.5E+00
    beta = -0.5E+00
    integr = 1
    a = 0.0E+00
    b = 0.1E+01
    limit = 200
    lenw = limit*4
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    CALL QAWS(F0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
      limit,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL QAWS(F0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
      2,8,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
    uflow = R1MACH(1)
    CALL QAWS(F0WS,a,b,alfa,beta,integr,uflow,0.0E+00,result,abserr,neval,ier,&
      limit,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==2.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
    CALL QAWS(F1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
      limit,lenw,last,iwork,work)
    ierv(1) = ier
    ierv(2) = 1
    ip = 0
    IF ( ier==3.OR.ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
    integr = 0
    CALL QAWS(F1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
      limit,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    IF ( ier==6 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQAWS FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQAWS PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQAWS
  !** CQNG
  SUBROUTINE CQNG(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for QNG.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (CQNG-S, CDQNG-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  CPRIN, F1N, F2N, QNG, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Added PASS/FAIL message and changed the name of the first
    !           argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER Lun
    REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, exact1, error, &
      exact2, result, uflow
    INTEGER ier, ierv(1), ip, Ipass, Kprint, neval
    DATA exact1/0.7281029132255818E+00/
    DATA exact2/0.1E+02/
    !* FIRST EXECUTABLE STATEMENT  CQNG
    IF ( Kprint>=2 ) WRITE (Lun,'(''1QNG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    epsabs = 0.0E+00
    epmach = R1MACH(4)
    uflow = R1MACH(1)
    epsrel = MAX(SQRT(epmach),0.1E-07)
    a = 0.0E+00
    b = 0.1E+01
    CALL QNG(F1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL CPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
    CALL QNG(F2N,a,b,uflow,0.0E+00,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    IF ( ier==1 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL CPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
    epsabs = 0.0E+00
    epsrel = 0.0E+00
    CALL QNG(F1N,a,b,epsabs,0.0E+00,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0 ) ip = 1
    IF ( ip==0 ) Ipass = 0
    IF ( Kprint/=0 ) CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF ( Kprint>=1 ) THEN
      IF ( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CQNG FAILED''/)')
      ELSEIF ( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CQNG PASSED''/)')
      ENDIF
    ENDIF
  END SUBROUTINE CQNG
  !** F0C
  REAL FUNCTION F0C(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F0C
    F0C = 1.E0/(X*X+1.E-4)
  END FUNCTION F0C
  !** F0F
  REAL FUNCTION F0F(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F0F
    F0F = 0.0
    IF ( X/=0.0 ) F0F = SIN(0.5E+02*X)/(X*SQRT(X))
  END FUNCTION F0F
  !** F0O
  REAL FUNCTION F0O(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F0O
    F0O = (2.0E0*SIN(X))**14
  END FUNCTION F0O
  !** F0S
  REAL FUNCTION F0S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F0S
    F0S = 0.0
    IF ( X/=0.0 ) F0S = 1.0/SQRT(X)
  END FUNCTION F0S
  !** F0WS
  REAL FUNCTION F0WS(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F0WS
    F0WS = SIN(10.0*X)
  END FUNCTION F0WS
  !** F1C
  REAL FUNCTION F1C(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F1C
    F1C = 0.0
    IF ( X/=0.33 ) F1C = (X-0.5)*ABS(X-0.33)**(-0.9)
  END FUNCTION F1C
  !** F1F
  REAL FUNCTION F1F(X)
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
    
    REAL X, x1, y
    !* FIRST EXECUTABLE STATEMENT  F1F
    x1 = X + 1.0
    F1F = 5.0/x1/x1
    y = 5.0/x1
    IF ( y>3.1415926535897932 ) F1F = 0.0
  END FUNCTION F1F
  !** F1G
  REAL FUNCTION F1G(X)
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
    
    REAL pi, X
    DATA pi/3.1415926535897932/
    !* FIRST EXECUTABLE STATEMENT  F1G
    F1G = 2.0/(2.0+SIN(10.0*pi*X))
  END FUNCTION F1G
  !** F1N
  REAL FUNCTION F1N(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F1N
    F1N = 1.0E0/(X**4+X**2+1.0E0)
  END FUNCTION F1N
  !** F1O
  REAL FUNCTION F1O(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F1O
    F1O = 1.0
    IF ( X>3.1415926535897932 ) F1O = 0.0
  END FUNCTION F1O
  !** F1P
  REAL FUNCTION F1P(X)
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
    
    REAL alfa1, alfa2, p1, p2, X, d1, d2
    !  P1 = 1/7, P2 = 2/3
    DATA p1/0.1428571428571428E+00/
    DATA p2/0.6666666666666667E+00/
    !* FIRST EXECUTABLE STATEMENT  F1P
    alfa1 = -0.25E0
    alfa2 = -0.5E0
    d1 = ABS(X-p1)
    d2 = ABS(X-p2)
    F1P = 0.0E+00
    IF ( d1/=0.0E+00.AND.d2/=0.0E+00 ) F1P = d1**alfa1 + d2**alfa2
  END FUNCTION F1P
  !** F1S
  REAL FUNCTION F1S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F1S
    F1S = 0.2E+01/(0.2E+01+SIN(0.314159E+02*X))
  END FUNCTION F1S
  !** F1WS
  REAL FUNCTION F1WS(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F1WS
    F1WS = ABS(X-0.33E+00)**(-0.999E+00)
  END FUNCTION F1WS
  !** F2G
  REAL FUNCTION F2G(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F2G
    F2G = X*SIN(0.3E+02*X)*COS(0.5E+02*X)
  END FUNCTION F2G
  !** F2N
  REAL FUNCTION F2N(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F2N
    F2N = X**(-0.9E+00)
  END FUNCTION F2N
  !** F2O
  REAL FUNCTION F2O(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F2O
    F2O = 0.0E+00
    IF ( X/=0.0E+00 ) F2O = 1.0/(X*X*SQRT(X))
  END FUNCTION F2O
  !** F2P
  REAL FUNCTION F2P(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F2P
    F2P = SIN(0.314159E+03*X)/(0.314159E+01*X)
  END FUNCTION F2P
  !** F2S
  REAL FUNCTION F2S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F2S
    F2S = 100.0
    IF ( X/=0.0 ) F2S = SIN(0.314159E+03*X)/(0.314159E+01*X)
  END FUNCTION F2S
  !** F3G
  REAL FUNCTION F3G(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F3G
    F3G = ABS(X-0.33E+00)**(-0.9E+00)
  END FUNCTION F3G
  !** F3P
  REAL FUNCTION F3P(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F3P
    F3P = 1.0
    IF ( X>3.1415926535897932 ) F3P = 0.0
  END FUNCTION F3P
  !** F3S
  REAL FUNCTION F3S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F3S
    F3S = 0.1E+01
    IF ( X>3.1415926535897932 ) F3S = 0.0
  END FUNCTION F3S
  !** F4P
  REAL FUNCTION F4P(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F4P
    F4P = 0.0
    IF ( X>0.0 ) F4P = 1.0/(X*SQRT(X))
  END FUNCTION F4P
  !** F4S
  REAL FUNCTION F4S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F4S
    IF ( X==.33E+00 ) THEN
      F4S = 0.0
      RETURN
    ENDIF
    F4S = ABS(X-0.33E+00)**(-0.999E+00)
    RETURN
    RETURN
  END FUNCTION F4S
  !** F5S
  REAL FUNCTION F5S(X)
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
    
    REAL X
    !* FIRST EXECUTABLE STATEMENT  F5S
    F5S = 0.0
    IF ( X/=0.0 ) F5S = 1.0/(X*SQRT(X))
  END FUNCTION F5S
  !** T0
  REAL FUNCTION T0(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F0S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T0
    a = 0.0E+00
    b = 0.1E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T0 = (b-a)*F0S(y)/x1/x1
  END FUNCTION T0
  !** T1
  REAL FUNCTION T1(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F1S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T1
    a = 0.0E+00
    b = 0.1E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T1 = (b-a)*F1S(y)/x1/x1
  END FUNCTION T1
  !** T2
  REAL FUNCTION T2(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F2S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T2
    a = 0.1E+00
    b = 0.1E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T2 = (b-a)*F2S(y)/x1/x1
  END FUNCTION T2
  !** T3
  REAL FUNCTION T3(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F3S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T3
    a = 0.0E+00
    b = 0.5E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T3 = (b-a)*F3S(y)/x1/x1
  END FUNCTION T3
  !** T4
  REAL FUNCTION T4(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F4S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T4
    a = 0.0E+00
    b = 0.1E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T4 = (b-a)*F4S(y)/x1/x1
  END FUNCTION T4
  !** T5
  REAL FUNCTION T5(X)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  F5S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    
    REAL a, b, X, x1, y
    !* FIRST EXECUTABLE STATEMENT  T5
    a = 0.0E+00
    b = 0.1E+01
    x1 = X + 0.1E+01
    y = (b-a)/x1 + a
    T5 = (b-a)*F5S(y)/x1/x1
  END FUNCTION T5
  !** CPRIN
  SUBROUTINE CPRIN(Lun,Num1,Kprint,Ip,Exact,Result,Abserr,Neval,Ierv,Lierv)
    IMPLICIT NONE
    !>
    !***
    !  Subsidiary to CQAG, CQAG, CQAGI, CQAGP, CQAGS, CQAWC,
    !            CQAWF, CQAWO, CQAWS, and CQNG.
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
    !   This program is called by the (single precision) Quadpack quick
    !   check routines for printing out their messages.
    !
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   810401  DATE WRITTEN
    !   890831  Modified array declarations.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   910627  Code completely rewritten.  (WRB)
    
    !     .. Scalar Arguments ..
    REAL Abserr, Exact, Result
    INTEGER Ip, Kprint, Lierv, Lun, Neval, Num1
    !     .. Array Arguments ..
    INTEGER Ierv(*)
    !     .. Local Scalars ..
    REAL error
    INTEGER ier, k
    !     .. Intrinsic Functions ..
    INTRINSIC ABS
    !* FIRST EXECUTABLE STATEMENT  CPRIN
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
    99009 FORMAT (' ',2(E15.7,1X),2(E9.2,1X),I4,4X,I6)
    99010 FORMAT (5X,'INFINITY',4X,E15.7,11X,E9.2,I5,4X,I6)
  END SUBROUTINE CPRIN
END MODULE TEST39_MOD
!** TEST39
PROGRAM TEST39
  USE TEST39_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2
  !***
  ! **Type:**      SINGLE PRECISION (TEST39-S, TEST40-D)
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
  !        QAG      QAGI     QAGP     QAGS     QAWC
  !        QAWF     QAWO     QAWS     QNG
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CQAG, CQAGI, CQAGP, CQAGS, CQAWC, CQAWF, CQAWO,
  !                    CQAWS, CQNG, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST39
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
  !     Test single precision QUADPACK routines
  !
  !     Test QAG.
  !
  CALL CQAG(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAGS.
  !
  CALL CQAGS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAGP.
  !
  CALL CQAGP(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAGI.
  !
  CALL CQAGI(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAWO.
  !
  CALL CQAWO(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAWF.
  !
  CALL CQAWF(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAWS.
  !
  CALL CQAWS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QAWC.
  !
  CALL CQAWC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test QNG.
  !
  CALL CQNG(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST39 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST39 *************')
  ENDIF
  STOP
END PROGRAM TEST39
