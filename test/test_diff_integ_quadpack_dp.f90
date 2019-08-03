MODULE TEST40_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** CDQAG
  SUBROUTINE CDQAG(Lun,Kprint,Ipass)
    !> Quick check for DQAG.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAG

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(2), Lun
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, result, work(400)
    INTEGER :: ier, ip, Ipass, iwork(100), key, Kprint, last, lenw, limit, neval
    REAL(DP), PARAMETER :: pi = 0.31415926535897932E+01_DP
    REAL(DP), PARAMETER :: exact1 = 0.1154700538379252E+01_DP
    REAL(DP), PARAMETER :: exact2 = 0.11780972450996172_DP
    REAL(DP), PARAMETER :: exact3 = 0.1855802E+02_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAG
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 100
    lenw = limit*4
    epsabs = 0._DP
    epmach = eps_dp
    key = 6
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    a = 0._DP
    b = 1._DP
    CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF( ier==0 .AND. error<=abserr .AND. abserr<=epsrel*ABS(exact1) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    limit = 1
!    lenw = limit*4
!    b = pi*2._DP
!    CALL DQAG(DF2G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
!    uflow = tiny_dp
!    limit = 100
!    lenw = limit*4
!    CALL DQAG(DF2G,a,b,uflow,0._DP,key,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
!    b = 1._DP
!    CALL DQAG(DF3G,a,b,epsabs,epsrel,1,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==3 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
!    lenw = 1
!    CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 .AND. &
!      last==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAG FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAG PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAG
  !** CDQAGI
  SUBROUTINE CDQAGI(Lun,Kprint,Ipass)
    !> Quick check for DQAGI.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAGI
    !
    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(4), inf
    REAL(DP) :: abserr, bound, epmach, epsabs, epsrel, error, result, work(800)
    INTEGER :: ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, Lun, neval
    REAL(DP), PARAMETER :: exact0 = 2._DP, exact1 = 0.115470066904E1_DP
    REAL(DP), PARAMETER :: exact2 = 0.909864525656E-02_DP
    REAL(DP), PARAMETER :: exact3 = 0.31415926535897932E+01_DP
    REAL(DP), PARAMETER :: exact4 = 0.19984914554328673E+04_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAGI
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAGI QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    bound = 0._DP
    inf = 1
    CALL DQAGI(DT0,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF( ier==0 .AND. error<=epsrel*ABS(exact0) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    CALL DQAGI(DT1,bound,inf,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
!      iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
!    uflow = tiny_dp
!    CALL DQAGI(DT2,bound,inf,uflow,0._DP,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==4 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
!    CALL DQAGI(DT3,bound,inf,uflow,0._DP,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 2
!    ip = 0
!    IF( ier==3 .OR. ier==4 .OR. ier==1 .OR. ier==2 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 4 OR 3 OR 1 OR 0
    !
!    CALL DQAGI(DT4,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 3
!    ierv(3) = 1
!    ierv(4) = 0
!    ip = 0
!    IF( ier==4 .OR. ier==3 .OR. ier==1 .OR. ier==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
!    oflow = huge_dp
!    CALL DQAGI(DT5,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==5 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
!    CALL DQAGI(DT1,bound,inf,epsabs,0._DP,result,abserr,neval,ier,limit,&
!      lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 .AND. &
!      last==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGI FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGI PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAGI
  !** CDQAGP
  SUBROUTINE CDQAGP(Lun,Kprint,Ipass)
    !> Quick check for DQAGP.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAGP

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(4)
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, &
      points(5), result, work(405)
    INTEGER :: ier, ip, Ipass, iwork(205), Kprint, last, leniw, lenw, limit, &
      Lun, neval, npts2
    REAL(DP), PARAMETER :: exact1 = 0.4285277667368085E+01_DP
    REAL(DP), PARAMETER :: exact2 = 0.909864525656E-2_DP
    REAL(DP), PARAMETER :: exact3 = 0.31415926535897932E+01_DP
    REAL(DP), PARAMETER :: p1 = 0.1428571428571428_DP
    REAL(DP), PARAMETER :: p2 = 0.6666666666666667_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAGP
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAGP QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    npts2 = 4
    limit = 100
    leniw = limit*2 + npts2
    lenw = limit*4 + npts2
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    a = 0._DP
    b = 1._DP
    points(1) = p1
    points(2) = p2
    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,lenw,last,iwork,work)
    error = ABS(result-exact1)
    ierv(1) = ier
    ip = 0
    IF( ier==0 .AND. error<=epsrel*ABS(exact1) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    leniw = 10
!    lenw = leniw*2 - npts2
!    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2, 4, 1 OR 3
    !
!    npts2 = 3
!    points(1) = 0.1_DP
!    leniw = limit*2 + npts2
!    lenw = limit*4 + npts2
!    uflow = tiny_dp
!    a = 0.1_DP
!    CALL DQAGP(DF2P,a,b,npts2,points,uflow,0._DP,result,abserr,neval,ier,&
!      leniw,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 3
!    ip = 0
!    IF( ier==2 .OR. ier==4 .OR. ier==1 .OR. ier==3 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
!    npts2 = 2
!    leniw = limit*2 + npts2
!    lenw = limit*4 + npts2
!    a = 0._DP
!    b = 5._DP
!    CALL DQAGP(DF3P,a,b,npts2,points,uflow,0._DP,result,abserr,neval,ier,&
!      leniw,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 2
!    ip = 0
!    IF( ier==3 .OR. ier==4 .OR. ier==1 .OR. ier==2 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
!    b = 1._DP
!    CALL DQAGP(DF4P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==5 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    oflow = huge_dp
!    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
!    npts2 = 5
!    leniw = limit*2 + npts2
!    lenw = limit*4 + npts2
!    points(1) = p1
!    points(2) = p2
!    points(3) = 0.3E+01_DP
!    CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 .AND. &
!      last==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGP FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGP PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAGP
  !** CDQAGS
  SUBROUTINE CDQAGS(Lun,Kprint,Ipass)
    !> Quick check for DQAGS.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   911114  Modified test on IER=4 to allow IER=5.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAGS

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(5), Lun
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, result, work(800)
    INTEGER :: ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    REAL(DP), PARAMETER :: exact0 = 2._DP
    REAL(DP), PARAMETER :: exact1 = 0.115470066904E+01_DP
    REAL(DP), PARAMETER :: exact2 = 0.909864525656E-02_DP
    REAL(DP), PARAMETER :: exact3 = 0.31415926535897932E+01_DP
    REAL(DP), PARAMETER :: exact4 = 0.19984914554328673E+04_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAGS
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAGS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    limit = 200
    lenw = limit*4
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    a = 0._DP
    b = 1._DP
    CALL DQAGS(DF0S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
      iwork,work)
    error = ABS(result-exact0)
    ierv(1) = ier
    ip = 0
    IF( ier==0 .AND. error<=epsrel*ABS(exact0) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    CALL DQAGS(DF1S,a,b,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
!    uflow = tiny_dp
!    a = 0.1_DP
!    CALL DQAGS(DF2S,a,b,uflow,0._DP,result,abserr,neval,ier,limit,lenw,last,&
!      iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==4 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
!    a = 0._DP
!    b = 5._DP
!    CALL DQAGS(DF3S,a,b,uflow,0._DP,result,abserr,neval,ier,limit,lenw,last,&
!      iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 2
!    ip = 0
!    IF( ier==3 .OR. ier==4 .OR. ier==1 .OR. ier==2 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 4, OR 5 OR 3 OR 1 OR 0
    !
!    b = 1._DP
!    CALL DQAGS(DF4S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
!      iwork,work)
!    ierv(1) = ier
!    ierv(2) = 5
!    ierv(3) = 3
!    ierv(4) = 1
!    ierv(5) = 0
!    ip = 0
!    IF( ier==5 .OR. ier==4 .OR. ier==3 .OR. ier==1 .OR. ier==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,5)
    !
    ! TEST ON IER = 5
    !
!    oflow = huge_dp
!    CALL DQAGS(DF5S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
!      iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==5 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
!    CALL DQAGS(DF1S,a,b,epsabs,0._DP,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 .AND. &
!      last==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGS FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGS PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAGS
  !** CDQAWC
  SUBROUTINE CDQAWC(Lun,Kprint,Ipass)
    !> Quick check for DQAWC.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAWC

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(2), Lun
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, c, result, work(800)
    INTEGER :: ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval
    REAL(DP), PARAMETER :: exact0 = -0.6284617285065624E+03_DP
    REAL(DP), PARAMETER :: exact1 = 0.1855802E+01_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAWC
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAWC QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    c = 0.5_DP
    a = -1._DP
    b = 1._DP
    limit = 200
    lenw = limit*4
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
      last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF( ier==0 .AND. error<=abserr .AND. abserr<=epsrel*ABS(exact0) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
!      iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
!    uflow = tiny_dp
!    CALL DQAWC(DF0C,a,b,c,uflow,0._DP,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
!    CALL DQAWC(DF1C,0._DP,b,c,uflow,0._DP,result,abserr,neval,ier,limit,&
!      lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==3 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
!    epsabs = 0._DP
!    epsrel = 0._DP
!    CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
!      last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWC FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWC PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAWC
  !** CDQAWF
  SUBROUTINE CDQAWF(Lun,Kprint,Ipass)
    !> Quick check for DQAWF.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAWF

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(4), integr, iwork(450), leniw, Lun, maxp1
    REAL(DP) :: a, abserr, epsabs, epmach, error, omega, result, work(1425)
    INTEGER :: ier, ip, Ipass, Kprint, lenw, limit, limlst, lst, neval
    REAL(DP), PARAMETER :: exact0 = 0.1422552162575912E+01_DP
    REAL(DP), PARAMETER :: pi = 0.31415926535897932E+01_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAWF
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAWF QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    limlst = 50
    limit = 200
    leniw = limit*2 + limlst
    lenw = leniw*2 + maxp1*25
    epmach = eps_dp
    epsabs = MAX(SQRT(epmach),0.1E-2_DP)
    a = 0._DP
    omega = 0.8E+01_DP
    integr = 2
    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
      leniw,maxp1,lenw,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF( ier==0 .AND. error<=abserr .AND. abserr<=epsabs ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    limlst = 3
!    leniw = 403
!    lenw = leniw*2 + maxp1*25
!    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
!      leniw,maxp1,lenw,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
!    limlst = 50
!    leniw = limit*2 + limlst
!    lenw = leniw*2 + maxp1*25
!    uflow = tiny_dp
!    CALL DQAWF(DF1F,a,0._DP,1,uflow,result,abserr,neval,ier,limlst,lst,&
!      leniw,maxp1,lenw,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 2
!    ip = 0
!    IF( ier==3 .OR. ier==4 .OR. ier==1 .OR. ier==2 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 6
    !
!    limlst = 50
!    leniw = 20
!    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
!      leniw,maxp1,lenw,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 7
    !
!    limlst = 50
!    leniw = 52
!    lenw = leniw*2 + maxp1*25
!    CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
!      leniw,maxp1,lenw,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==7 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,7,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWF FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWF PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAWF
  !** CDQAWO
  SUBROUTINE CDQAWO(Lun,Kprint,Ipass)
    !> Quick check for DQAWO.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAWO

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: leniw
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, omega, result, work(1325)
    INTEGER :: ier, ierv(4), integr, ip, Ipass, iwork(400), Kprint, last, lenw, &
      Lun, maxp1, neval
    REAL(DP), PARAMETER :: exact0 = 0.1042872789432789E+05_DP
    REAL(DP), PARAMETER :: pi = 0.31415926535897932E+01_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAWO
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAWO QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    maxp1 = 21
    leniw = 400
    lenw = leniw*2 + maxp1*25
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    a = 0._DP
    b = pi
    omega = 1._DP
    integr = 2
    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
      leniw,maxp1,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF( ier==0 .AND. error<=abserr .AND. abserr<=epsrel*ABS(exact0) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    leniw = 2
!    lenw = leniw*2 + maxp1*25
!    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,maxp1,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 4 OR 1
    !
!    uflow = tiny_dp
!    leniw = 400
!    lenw = leniw*2 + maxp1*25
!    CALL DQAWO(DF0O,a,b,omega,integr,uflow,0._DP,result,abserr,neval,ier,&
!      leniw,maxp1,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==4 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,3)
    !
    ! TEST ON IER = 3 OR 4 OR 1 OR 2
    !
!    b = 5._DP
!    omega = 0._DP
!    integr = 1
!    CALL DQAWO(DF1O,a,b,omega,integr,uflow,0._DP,result,abserr,neval,ier,&
!      leniw,maxp1,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 4
!    ierv(3) = 1
!    ierv(4) = 2
!    ip = 0
!    IF( ier==3 .OR. ier==4 .OR. ier==1 .OR. ier==2 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
    !
    ! TEST ON IER = 5
    !
!    b = 1._DP
!    oflow = huge_dp
!    CALL DQAWO(DF2O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,maxp1,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==5 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
!    integr = 3
!    CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
!      leniw,maxp1,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 .AND. &
!      last==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWO FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWO PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAWO
  !** CDQAWS
  SUBROUTINE CDQAWS(Lun,Kprint,Ipass)
    !> Quick check for DQAWS.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : eps_dp
    USE diff_integ, ONLY : DQAWS

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: ierv(2), Lun
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, alfa, beta, &
      result, work(800)
    INTEGER :: ier, ip, Ipass, iwork(200), Kprint, last, lenw, limit, neval, integr
    REAL(DP), PARAMETER :: exact0 = 0.5350190569223644_DP
    REAL(DP), PARAMETER :: exact1 = 0.1998491554328673E+04_DP
    !* FIRST EXECUTABLE STATEMENT  CDQAWS
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQAWS QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    alfa = -0.5_DP
    beta = -0.5_DP
    integr = 1
    a = 0._DP
    b = 1._DP
    limit = 200
    lenw = limit*4
    epsabs = 0._DP
    epmach = eps_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
      ier,limit,lenw,last,iwork,work)
    ierv(1) = ier
    ip = 0
    error = ABS(exact0-result)
    IF( ier==0 .AND. error<=epsrel*ABS(exact0) ) ip = 1
    IF( ip==0 ) Ipass = 0
    CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
!      ier,2,8,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 2 OR 1
    !
!    uflow = tiny_dp
!    CALL DQAWS(DF0WS,a,b,alfa,beta,integr,uflow,0._DP,result,abserr,neval,&
!      ier,limit,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==2 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 3 OR 1
    !
!    CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
!      ier,limit,lenw,last,iwork,work)
!    ierv(1) = ier
!    ierv(2) = 1
!    ip = 0
!    IF( ier==3 .OR. ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
    !
    ! TEST ON IER = 6
    !
!    integr = 0
!    CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
!      ier,limit,lenw,last,iwork,work)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWS FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWS PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQAWS
  !** CDQNG
  SUBROUTINE CDQNG(Lun,Kprint,Ipass)
    !> Quick check for DQNG.
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
    !   901205  Added PASS/FAIL message and changed the name of the first argument.  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    USE service, ONLY : tiny_dp, eps_dp
    USE diff_integ, ONLY : DQNG

    ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
    INTEGER :: Lun
    REAL(DP) :: a, abserr, b, epmach, epsabs, epsrel, error, result, uflow
    INTEGER :: ier, ierv(1), ip, Ipass, Kprint, neval
    REAL(DP), PARAMETER :: exact1 = 0.7281029132255818_DP
    REAL(DP), PARAMETER :: exact2 = 0.1E+02_DP
    !* FIRST EXECUTABLE STATEMENT  CDQNG
    IF( Kprint>=2 ) WRITE (Lun,'(''1DQNG QUICK CHECK''/)')
    !
    ! TEST ON IER = 0
    !
    Ipass = 1
    epsabs = 0._DP
    epmach = eps_dp
    uflow = tiny_dp
    epsrel = MAX(SQRT(epmach),0.1E-7_DP)
    a = 0._DP
    b = 1._DP
    CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
    CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
    ierv(1) = ier
    ip = 0
    error = ABS(exact1-result)
    IF( ier==0 .AND. error<=abserr .AND. abserr<=epsrel*ABS(exact1) ) ip = 1
    IF( ip==0 ) Ipass = 0
    IF( Kprint/=0 ) CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 1
    !
!    CALL DQNG(DF2N,a,b,uflow,0._DP,result,abserr,neval,ier)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==1 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    IF( Kprint/=0 ) CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
    !
    ! TEST ON IER = 6
    !
!    epsabs = 0._DP
!    epsrel = 0._DP
!    CALL DQNG(DF1N,a,b,epsabs,0._DP,result,abserr,neval,ier)
!    ierv(1) = ier
!    ip = 0
!    IF( ier==6 .AND. result==0._DP .AND. abserr==0._DP .AND. neval==0 ) ip = 1
!    IF( ip==0 ) Ipass = 0
!    IF( Kprint/=0 ) CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
    !
    IF( Kprint>=1 ) THEN
      IF( Ipass==0 ) THEN
        WRITE (Lun,'(/'' SOME TEST(S) IN CDQNG FAILED''/)')
      ELSEIF( Kprint>=2 ) THEN
        WRITE (Lun,'(/'' ALL TEST(S) IN CDQNG PASSED''/)')
      END IF
    END IF
  END SUBROUTINE CDQNG
  !** DF0C
  REAL(DP) PURE FUNCTION DF0C(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0C
    DF0C = 1._DP/(X*X+1.E-4_DP)
  END FUNCTION DF0C
  !** DF0F
  REAL(DP) PURE FUNCTION DF0F(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0F
    DF0F = 0._DP
    IF( X/=0._DP ) DF0F = SIN(0.5E+02_DP*X)/(X*SQRT(X))
  END FUNCTION DF0F
  !** DF0O
  REAL(DP) PURE FUNCTION DF0O(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0O
    DF0O = (2._DP*SIN(X))**14
  END FUNCTION DF0O
  !** DF0S
  REAL(DP) PURE FUNCTION DF0S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0S
    DF0S = 0._DP
    IF( X/=0._DP ) DF0S = 1._DP/SQRT(X)
  END FUNCTION DF0S
  !** DF0WS
  REAL(DP) PURE FUNCTION DF0WS(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF0WS
    DF0WS = SIN(0.1E+02_DP*X)
  END FUNCTION DF0WS
  !** DF1C
  REAL(DP) PURE FUNCTION DF1C(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1C
    DF1C = 0._DP
    IF( X/=0.33_DP ) DF1C = (X-0.5_DP)*ABS(X-0.33_DP)**(-0.9_DP)
  END FUNCTION DF1C
  !** DF1F
  REAL(DP) PURE FUNCTION DF1F(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: x1, y
    !* FIRST EXECUTABLE STATEMENT  DF1F
    x1 = X + 1._DP
    DF1F = 5._DP/x1/x1
    y = 5._DP/x1
    IF( y>3.1415926535897932_DP ) DF1F = 0._DP
  END FUNCTION DF1F
  !** DF1G
  REAL(DP) PURE FUNCTION DF1G(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP), PARAMETER :: pi = 3.1415926535897932_DP
    !* FIRST EXECUTABLE STATEMENT  DF1G
    DF1G = 2._DP/(2._DP+SIN(10._DP*pi*X))
  END FUNCTION DF1G
  !** DF1N
  REAL(DP) PURE FUNCTION DF1N(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1N
    DF1N = 1._DP/(X**4+X**2+1._DP)
  END FUNCTION DF1N
  !** DF1O
  REAL(DP) PURE FUNCTION DF1O(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1O
    DF1O = 1._DP
    IF( X>0.31415926535897932E+01_DP ) DF1O = 0._DP
  END FUNCTION DF1O
  !** DF1P
  REAL(DP) PURE FUNCTION DF1P(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: alfa1, alfa2, d1, d2
    !* FIRST EXECUTABLE STATEMENT  DF1P
    !  P1 = 1/7, P2 = 2/3
    REAL(DP), PARAMETER :: p1 = 0.1428571428571428_DP
    REAL(DP), PARAMETER :: p2 = 0.6666666666666667_DP
    alfa1 = -0.25_DP
    alfa2 = -0.5_DP
    d1 = ABS(X-p1)
    d2 = ABS(X-p2)
    DF1P = 0._DP
    IF( d1/=0._DP .AND. d2/=0._DP ) DF1P = d1**alfa1 + d2**alfa2
  END FUNCTION DF1P
  !** DF1S
  REAL(DP) PURE FUNCTION DF1S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1S
    DF1S = 2._DP/(2._DP+SIN(0.314159E+02_DP*X))
  END FUNCTION DF1S
  !** DF1WS
  REAL(DP) PURE FUNCTION DF1WS(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF1WS
    DF1WS = 0.00_DP
    IF( X-0.33_DP/=0._DP ) DF1WS = ABS(X-0.33_DP)**(-0.999_DP)
  END FUNCTION DF1WS
  !** DF2G
  REAL(DP) PURE FUNCTION DF2G(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2G
    DF2G = X*SIN(0.3E+02_DP*X)*COS(0.5E+02_DP*X)
  END FUNCTION DF2G
  !** DF2N
  REAL(DP) PURE FUNCTION DF2N(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2N
    DF2N = X**(-0.9_DP)
  END FUNCTION DF2N
  !** DF2O
  REAL(DP) PURE FUNCTION DF2O(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2O
    DF2O = 0._DP
    IF( X/=0._DP ) DF2O = 1._DP/(X*X*SQRT(X))
  END FUNCTION DF2O
  !** DF2P
  REAL(DP) PURE FUNCTION DF2P(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2P
    DF2P = SIN(0.314159E+03_DP*X)/(0.314159E+01_DP*X)
  END FUNCTION DF2P
  !** DF2S
  REAL(DP) PURE FUNCTION DF2S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF2S
    DF2S = 0.1E+03_DP
    IF( X/=0._DP ) DF2S = SIN(0.314159E+03_DP*X)/(0.314159E+01_DP*X)
  END FUNCTION DF2S
  !** DF3G
  REAL(DP) PURE FUNCTION DF3G(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3G
    DF3G = ABS(X-0.33_DP)**(-.90_DP)
  END FUNCTION DF3G
  !** DF3P
  REAL(DP) PURE FUNCTION DF3P(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3P
    DF3P = 1._DP
    IF( X>0.31415926535897932E+01_DP ) DF3P = 0._DP
  END FUNCTION DF3P
  !** DF3S
  REAL(DP) PURE FUNCTION DF3S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF3S
    DF3S = 1._DP
    IF( X>0.31415926535897932E+01_DP ) DF3S = 0._DP
  END FUNCTION DF3S
  !** DF4P
  REAL(DP) PURE FUNCTION DF4P(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF4P
    DF4P = 0._DP
    IF( X>0._DP ) DF4P = 1._DP/(X*SQRT(X))
  END FUNCTION DF4P
  !** DF4S
  REAL(DP) PURE FUNCTION DF4S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF4S
    DF4S = 0.00_DP
    IF( X-0.33_DP/=0._DP ) DF4S = ABS(X-0.33_DP)**(-0.999_DP)
  END FUNCTION DF4S
  !** DF5S
  REAL(DP) PURE FUNCTION DF5S(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    !* FIRST EXECUTABLE STATEMENT  DF5S
    DF5S = 0._DP
    IF( X/=0._DP ) DF5S = 1._DP/(X*SQRT(X))
  END FUNCTION DF5S
  !** DT0
  REAL(DP) PURE FUNCTION DT0(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF0S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT0
    a = 0._DP
    b = 1._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT0 = (b-a)*DF0S(y)/x1/x1
  END FUNCTION DT0
  !** DT1
  REAL(DP) PURE FUNCTION DT1(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF1S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT1
    a = 0._DP
    b = 1._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT1 = (b-a)*DF1S(y)/x1/x1
  END FUNCTION DT1
  !** DT2
  REAL(DP) PURE FUNCTION DT2(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF2S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT2
    a = 0.1_DP
    b = 1._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT2 = (b-a)*DF2S(y)/x1/x1
  END FUNCTION DT2
  !** DT3
  REAL(DP) PURE FUNCTION DT3(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF3S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT3
    a = 0._DP
    b = 5._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT3 = (b-a)*DF3S(y)/x1/x1
  END FUNCTION DT3
  !** DT4
  REAL(DP) PURE FUNCTION DT4(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF4S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT4
    a = 0._DP
    b = 1._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT4 = (b-a)*DF4S(y)/x1/x1
  END FUNCTION DT4
  !** DT5
  REAL(DP) PURE FUNCTION DT5(X)
    !> Subsidiary to
    !***
    ! **Library:**   SLATEC
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  DF5S

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)

    REAL(DP), INTENT(IN) :: X
    REAL(DP) :: a, b, x1, y
    !* FIRST EXECUTABLE STATEMENT  DT5
    a = 0._DP
    b = 1._DP
    x1 = X + 1._DP
    y = (b-a)/x1 + a
    DT5 = (b-a)*DF5S(y)/x1/x1
  END FUNCTION DT5
  !** DPRIN
  SUBROUTINE DPRIN(Lun,Num1,Kprint,Ip,Exact,Result,Abserr,Neval,Ierv,Lierv)
    !> Subsidiary to CDQAG, CDQAG, CDQAGI, CDQAGP, CDQAGS, CDQAWC,
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
    REAL(DP) :: Abserr, Exact, Result
    INTEGER :: Ip, Kprint, Lierv, Lun, Neval, Num1
    !     .. Array Arguments ..
    INTEGER :: Ierv(Lierv)
    !     .. Local Scalars ..
    REAL(DP) :: error
    INTEGER :: ier, k
    !     .. Intrinsic Functions ..
    INTRINSIC ABS
    !* FIRST EXECUTABLE STATEMENT  DPRIN
    ier = Ierv(1)
    error = ABS(Exact-Result)
    !
    IF( Kprint>=2 ) THEN
      IF( Ip/=1 ) THEN
        !
        !         Write failure messages.
        !
        WRITE (UNIT=Lun,FMT=99002) Num1
        IF( Num1==0 ) WRITE (UNIT=Lun,FMT=99003)
        IF( Num1>0 ) WRITE (UNIT=Lun,FMT=99004) Num1
        IF( Lierv>1 ) WRITE (UNIT=Lun,FMT=99005) (Ierv(k),k=2,Lierv)
        IF( Num1==6 ) WRITE (UNIT=Lun,FMT=99006)
        WRITE (UNIT=Lun,FMT=99007)
        WRITE (UNIT=Lun,FMT=99008)
        IF( Num1/=5 ) THEN
          WRITE (UNIT=Lun,FMT=99009) Exact, Result, error, Abserr, ier, Neval
        ELSE
          WRITE (Lun,FMT=99010) Result, Abserr, ier, Neval
        END IF
      ELSEIF( Kprint>=3 ) THEN
        !
        !           Write PASS message.
        !
        WRITE (UNIT=Lun,FMT=99001) Num1
      END IF
    END IF
    !
    RETURN
    !
    99001 FORMAT (' TEST ON IER = ',I2,' PASSED')
    99002 FORMAT (' TEST ON IER = ',I1,' FAILED.')
    99003 FORMAT (' WE MUST HAVE IER = 0, ERROR<=ABSERR AND ABSERR.LE',&
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
  USE TEST40_MOD, ONLY : CDQAG, CDQAGI, CDQAGP, CDQAGS, CDQAWC, CDQAWF, CDQAWO, &
    CDQAWS, CDQNG
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
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
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST40
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test double precision QUADPACK routines
  !
  !     Test DQAG.
  !
  CALL CDQAG(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGS.
  !
  CALL CDQAGS(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGP.
  !
  CALL CDQAGP(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAGI.
  !
  CALL CDQAGI(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWO.
  !
  CALL CDQAWO(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWF.
  !
  CALL CDQAWF(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWS.
  !
  CALL CDQAWS(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQAWC.
  !
  CALL CDQAWC(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQNG.
  !
  CALL CDQNG(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST40 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST40 *************')
  END IF
  STOP
END PROGRAM TEST40
