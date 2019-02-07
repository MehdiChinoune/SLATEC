!*==FCNQX1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FCNQX1
SUBROUTINE FCNQX1(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--FCNQX15
  !*** Start of declarations inserted by SPAG
  INTEGER i, I1MACH, ic1, ic2, id, ierr, ierror, ip, ipn, iq, ir ,&
    irad, isig, isum, ix11, ix12, ix13, ix21, ix22, ix23
  INTEGER mu, mu1, mu2, n, nbits, ndec, nerr, nradpl, nu1, nudiff ,&
    NUMXER
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  FCNQX1
  !***SUBSIDIARY
  !***PURPOSE  THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES XLEGF
  !            AND XNRMP WHICH CALCULATE LEGENDRE FUNCTIONS
  !***LIBRARY   SLATEC
  !***CATEGORY  C3A2, C9
  !***TYPE      SINGLE PRECISION (FCNQX1-S, FCNQX2-D)
  !***KEYWORDS  LEGENDRE FUNCTIONS, QUICK CHECK
  !***AUTHOR  LOZIER, DANIEL W., (NIST)
  !           SMITH, JOHN M., (NIST AND GEORGE MASON UNIVERSITY)
  !***REFERENCES  OLVER AND SMITH,J.COMPUT.PHYSICS,51(1983),NO.3,502-518.
  !               SMITH, OLVER AND LOZIER,ACM TRANS MATH SOFTW,7(1981),
  !                 NO.1,93-105.
  !***ROUTINES CALLED  XCON, XCSRT, XERCLR, XLEGF, XNRMP, XSET, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   881020  DATE WRITTEN
  !   900306  Added SLATEC prologue to this routine. (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !   910104  Changed to print variable number of decimals. (DWL and JMS)
  !***END PROLOGUE  FCNQX1
  !
  CHARACTER(34) :: fmt, fmtf, fmti
  INTEGER Lun, Kprint, Ipass
  DIMENSION p(10), q(10), r(10), c1(10), c2(10), ip(10), iq(10) ,&
    ir(10)
  DIMENSION ic1(10), ic2(10), pn(10), ipn(10)
  REAL p, q, r, c1, c2, pn
  REAL deg, theta, dnu1, dzero
  REAL x11, x12, x13, x21, x22, x23
  REAL nu
  !
  !***FIRST EXECUTABLE STATEMENT  FCNQX1
  !
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (' ** TEST SINGLE PRECISION LEGENDRE FUNCTION ROUTINES',&
    ' IN FCNPAK ** ',/)
  Ipass = 1
  irad = 0
  nradpl = 0
  dzero = 0.0
  nbits = 0
  CALL XSET(irad,nradpl,dzero,nbits,ierror)
  IF ( ierror/=0 ) Ipass = 0
  ierr = 0
  dnu1 = 2000.4
  IF ( I1MACH(13)*LOG10(REAL(I1MACH(10)))<150. ) dnu1 = 100.4
  IF ( Kprint>2 ) THEN
    IF ( I1MACH(13)<500 ) WRITE (Lun,99002)
    99002   FORMAT (' ON COMPUTERS WITH MAXIMUM EXPONENT LESS THAN 500, SMALL'/&
      ' TEST VALUES FOR NU, MU ARE USED. IF LARGER THAN OR EQUAL 500,'&
      /&
      ' LARGER VALUES ARE USED. THIS COMPUTER USES THE SMALLER VALUES.'&
      )
    IF ( I1MACH(13)>=500 ) WRITE (Lun,99003)
    99003   FORMAT (' ON COMPUTERS WITH MAXIMUM EXPONENT LESS THAN 500, SMALL'/&
      ' TEST VALUES FOR NU, MU ARE USED. IF LARGER THAN OR EQUAL 500,'&
      /&
      ' LARGER VALUES ARE USED. THIS COMPUTER USES THE LARGER VALUES.'&
      )
  ENDIF
  nudiff = 5
  mu1 = dnu1
  mu2 = mu1
  deg = 0.1
  theta = deg*4.*ATAN(1.0)/180.0
  !
  ! In TEST 1 the Legendre functions P (of both positive and negative
  ! order) and Q are calculated.  Large values of mu and nu are used
  ! so that it is necessary to use extended range arithmetic.  The
  ! values of the Casoratians should be approximately equal to 1.0.
  ! The check which is applied is to verify that the difference between
  ! the Casoratians and 1.0 is less that 10.**(6-NDEC), where NDEC =
  ! INT((D-1)*LOG10(R)), D = I1MACH(11) = significand length, R =
  ! I1MACH(10) = radix. The value of IERROR should always be returned
  ! as zero. This test uses the programs
  ! XLEGF, XPQNU, XPSI, XQNU, XPMUP, XSET, XADD,
  ! XADJ, XCSRT, XRED, XC210, and XCON.
  !
  isum = 0
  ndec = (I1MACH(11)-1)*LOG10(REAL(I1MACH(10)))
  ! Formats that depend on NDEC ...
  fmt(1:20) = '(1X, 6X, 4H   (,E30.'
  WRITE (fmt(21:22),'(I2)') ndec
  fmt(23:34) = ',1H,,I8,1H))'
  fmtf(1:20) = '(1X,F6.1,4H   (,E30.'
  WRITE (fmtf(21:22),'(I2)') ndec
  fmtf(23:34) = ',1H,,I8,1H))'
  fmti(1:20) = '(1X, I6, 4H   (,E30.'
  WRITE (fmti(21:22),'(I2)') ndec
  fmti(23:34) = ',1H,,I8,1H))'
  IF ( Kprint>2 ) WRITE (Lun,99004) mu1, deg
  99004 FORMAT (/' TEST 1, FIXED MU = ',I4,' AND THETA = ',F3.1,&
    ' DEGREES, RECURRENCE IN NU,'/'         CASORATIS SHOULD = 1.0')
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,1,p,ip,ierror)
  isum = isum + ierror
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,2,q,iq,ierror)
  isum = isum + ierror
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,3,r,ir,ierror)
  isum = isum + ierror
  CALL XCSRT(dnu1,nudiff,mu1,mu2,theta,p,q,r,ip,iq,ir,c1,ic1,c2,ic2,ierror)
  isum = isum + ierror
  DO i = 1, 6
    CALL XCON(p(i),ip(i),ierror)
    isum = isum + ierror
    CALL XCON(q(i),iq(i),ierror)
    isum = isum + ierror
    CALL XCON(r(i),ir(i),ierror)
    isum = isum + ierror
  ENDDO
  x11 = p(1)
  ix11 = ip(1)
  x12 = r(1)
  ix12 = ir(1)
  x13 = q(1)
  ix13 = iq(1)
  IF ( Kprint>2 ) THEN
    WRITE (Lun,'(A)') '     NU   CASORATI 1'
    nu = dnu1
    DO i = 1, 5
      WRITE (Lun,fmtf) nu, c1(i), ic1(i)
      nu = nu + 1.
    ENDDO
    WRITE (Lun,'(A)') '     NU   CASORATI 2'
    nu = dnu1
    DO i = 1, 5
      WRITE (Lun,fmtf) nu, c2(i), ic2(i)
      nu = nu + 1.
    ENDDO
  ENDIF
  DO i = 1, 5
    IF ( ABS(1.0-c1(i))>=10.0E0**(6-ndec) ) GOTO 100
    IF ( ABS(1.0-c2(i))>=10.0E0**(6-ndec) ) GOTO 100
  ENDDO
  IF ( isum==0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99005)
    99005   FORMAT (' ***** TEST 1 (SINGLE PRECISION) PASSED *****')
    GOTO 200
  ENDIF
  100 CONTINUE
  IF ( Kprint>=1 ) WRITE (Lun,99006)
  99006 FORMAT (' ***** TEST 1 (SINGLE PRECISION) FAILED *****')
  ierr = ierr + 1
  Ipass = 0
  200  nudiff = 0
  mu1 = mu2 - 5
  !
  ! In TEST 2 P (of positive and negative order) and Q are again
  ! calculated but in this test the recurrence is in the mu-wise direction
  ! rather than in the nu-wise direction as was the case before.  The same
  ! programs are used except that XQNU is not used and XQMU and XPMU
  ! are used. Again the criterion for passing the test is that the
  ! Casoratians differ from 1.0 by less than 10.0**(6-NDEC). The value
  ! of IERROR should always be returned as zero.
  !
  isum = 0
  IF ( Kprint>2 ) WRITE (Lun,99007) dnu1, deg
  99007 FORMAT (/' TEST 2, FIXED NU = ',F6.1,' AND THETA = ',F3.1,&
    ' DEGREES, RECURRENCE IN MU,'/'         CASORATIS SHOULD = 1.0')
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,1,p,ip,ierror)
  isum = isum + ierror
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,2,q,iq,ierror)
  isum = isum + ierror
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,3,r,ir,ierror)
  isum = isum + ierror
  CALL XCSRT(dnu1,nudiff,mu1,mu2,theta,p,q,r,ip,iq,ir,c1,ic1,c2,ic2,ierror)
  isum = isum + ierror
  DO i = 1, 6
    CALL XCON(p(i),ip(i),ierror)
    isum = isum + ierror
    CALL XCON(q(i),iq(i),ierror)
    isum = isum + ierror
    CALL XCON(r(i),ir(i),ierror)
    isum = isum + ierror
  ENDDO
  x21 = p(6)
  ix21 = ip(6)
  x22 = r(6)
  ix22 = ir(6)
  x23 = q(6)
  ix23 = iq(6)
  IF ( Kprint>2 ) THEN
    WRITE (Lun,'(A)') '     MU   CASORATI 3'
    mu = mu1
    DO i = 1, 5
      WRITE (Lun,fmti) mu, c1(i), ic1(i)
      mu = mu + 1
    ENDDO
    WRITE (Lun,'(A)') '     MU   CASORATI 4'
    mu = mu1
    DO i = 1, 5
      WRITE (Lun,fmti) mu, c2(i), ic2(i)
      mu = mu + 1
    ENDDO
  ENDIF
  DO i = 1, 5
    IF ( ABS(1.0-c1(i))>=10.0E0**(6-ndec) ) GOTO 300
    IF ( ABS(1.0-c2(i))>=10.0E0**(6-ndec) ) GOTO 300
    IF ( isum/=0 ) GOTO 300
  ENDDO
  IF ( Kprint>=2 ) WRITE (Lun,99008)
  99008 FORMAT (' ***** TEST 2 (SINGLE PRECISION) PASSED *****')
  GOTO 400
  300 CONTINUE
  IF ( Kprint>=1 ) WRITE (Lun,99009)
  99009 FORMAT (' ***** TEST 2 (SINGLE PRECISION) FAILED *****')
  ierr = ierr + 1
  Ipass = 0
  !
  ! In TEST 3 values of P and Q which were calculated in two different
  ! manners, one by nu-wise recurrence in TEST 1 and one by mu-wise
  ! recurrence in TEST 2, are compared.  Again, the criterion for success
  ! is a relative error of less than 10.0**(6-NDEC).
  !
  400 CONTINUE
  IF ( Kprint>2 ) THEN
    WRITE (Lun,99010) deg, mu2, dnu1
    99010   FORMAT (/' TEST 3, COMPARISON OF VALUES FROM TEST 1 AND TEST 2',&
      ' WITH THETA = ',F3.1,' DEGREES,'/'         MU = ',I4,&
      ' AND NU = ',F6.1)
    WRITE (Lun,'(A)') '          P(-MU,NU)'
    WRITE (Lun,fmt) x11, ix11
    WRITE (Lun,fmt) x21, ix21
    WRITE (Lun,'(A)') '          P(MU,NU)'
    WRITE (Lun,fmt) x12, ix12
    WRITE (Lun,fmt) x22, ix22
    WRITE (Lun,'(A)') '          Q(MU,NU)'
    WRITE (Lun,fmt) x13, ix13
    WRITE (Lun,fmt) x23, ix23
  ENDIF
  IF ( ABS((x11-x21)/x11)<10.0E0**(6-ndec) ) THEN
    IF ( ABS((x12-x22)/x12)<10.0E0**(6-ndec) ) THEN
      IF ( ABS((x13-x13)/x13)<10.0E0**(6-ndec) ) THEN
        IF ( ix11==ix21 ) THEN
          IF ( ix12==ix22 ) THEN
            IF ( ix13==ix23 ) THEN
              IF ( Kprint>=2 ) WRITE (Lun,99011)
              99011             FORMAT (' ***** TEST 3 (SINGLE PRECISION) PASSED *****')
              GOTO 500
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  IF ( Kprint>=1 ) WRITE (Lun,99012)
  99012 FORMAT (' ***** TEST 3 (SINGLE PRECISION) FAILED *****')
  ierr = ierr + 1
  Ipass = 0
  !
  ! In TEST 4, the value of the normalized Legendre function as
  ! calculated by XLEGF and XPNRM is compared to the same value
  ! as calculated by the program XNRMP.  Again the criterion is a
  ! relative error of less than 10.0**(6-NDEC). The value of IERROR
  ! should always be returned as zero.
  !
  500  isum = 0
  dnu1 = 100.0
  nudiff = 0
  mu1 = 10
  mu2 = 10
  IF ( Kprint>2 ) WRITE (Lun,99013) deg, mu1, dnu1
  99013 FORMAT (/' TEST 4, COMPARISON OF VALUES FROM XLEGF AND XNRMP',&
    ' WITH THETA = ',F3.1,' DEGREES,'/'         MU = ',I4,&
    ' AND NU = ',F6.1)
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,4,pn,ipn,ierror)
  isum = isum + ierror
  x11 = pn(1)
  ix11 = ipn(1)
  nu1 = 100
  CALL XNRMP(nu1,mu1,mu2,theta,2,pn,ipn,isig,ierror)
  isum = isum + ierror
  x21 = pn(1)
  ix21 = ipn(1)
  IF ( Kprint>2 ) THEN
    WRITE (Lun,'(A)') '          NORMALIZED P'
    WRITE (Lun,fmt) x11, ix11
    WRITE (Lun,fmt) x21, ix21
  ENDIF
  IF ( ABS((x11-x21)/x11)<10.0E0**(6-ndec) ) THEN
    IF ( ix11==ix21 ) THEN
      IF ( isum==0 ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99014)
        99014       FORMAT (' ***** TEST 4 (SINGLE PRECISION) PASSED *****')
        GOTO 600
      ENDIF
    ENDIF
  ENDIF
  IF ( Kprint>=1 ) WRITE (Lun,99015)
  99015 FORMAT (' ***** TEST 4 (SINGLE PRECISION) FAILED *****')
  ierr = ierr + 1
  Ipass = 0
  !
  ! In TEST 5 errors are purposely made in input so as to test error
  ! handling capability. First, an incorrect value of ID is given. Then
  ! both NUDIFF and MU2-MU1 are non-zero. Finally, an incorrect value
  ! of THETA is given. In each case the value of the error indicator
  ! IERROR should equal the error number as returned by the error
  ! handling package XERROR (which includes XSETF, XERCLR, and NUMXER).
  !
  600  CALL XSETF(-1)
  IF ( Kprint<=2 ) CALL XSETF(0)
  IF ( Kprint>2 ) WRITE (Lun,99016)
  99016 FORMAT (/' TEST 5, TEST OF ERROR HANDLING. 3 ERROR MESSAGES',&
    ' SHOULD BE PRINTED.')
  nudiff = 0
  mu2 = mu1
  id = 5
  CALL XERCLR
  CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,id,p,ip,ierror)
  n = NUMXER(nerr)
  IF ( n==ierror ) THEN
    mu2 = mu1 + 5
    nudiff = 5
    CALL XERCLR
    CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,1,p,ip,ierror)
    n = NUMXER(nerr)
    IF ( n==ierror ) THEN
      nudiff = 0
      theta = 2.0
      CALL XERCLR
      CALL XLEGF(dnu1,nudiff,mu1,mu2,theta,1,p,ip,ierror)
      n = NUMXER(nerr)
      IF ( n==ierror ) THEN
        IF ( Kprint>=2 ) WRITE (Lun,99017)
        99017       FORMAT (' ***** TEST 5 (SINGLE PRECISION) PASSED *****')
        GOTO 700
      ENDIF
    ENDIF
  ENDIF
  IF ( Kprint>=1 ) WRITE (Lun,99018)
  99018 FORMAT (' ***** TEST 5 (SINGLE PRECISION) FAILED *****')
  ierr = ierr + 1
  Ipass = 0
  700 CONTINUE
  IF ( ierr/=0 ) THEN
    IF ( Kprint>=2 ) WRITE (Lun,99019) ierr
    99019   FORMAT (/'  TESTS COMPLETED, NUMBER OF TESTS FAILED = ',I2)
  ENDIF
END SUBROUTINE FCNQX1
