!DECK CQCBI
SUBROUTINE CQCBI(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CQCBI
  !***SUBSIDIARY
  !***PURPOSE  Quick check for SLATEC subroutine
  !            CBESI
  !***LIBRARY   SLATEC
  !***CATEGORY  C10B4
  !***TYPE      COMPLEX (CQCBI-C, ZQCBI-Z)
  !***KEYWORDS  QUICK CHECK, CBESI
  !***AUTHOR  Amos, Don, (SNL)
  !           Goudy, Sue, (SNL)
  !           Walton, Lee, (SNL)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL CQCBI (LUN, KPRINT, IPASS)
  !
  ! *Arguments:
  !
  !     LUN    :IN  is the unit number to which output is to be written.
  !
  !     KPRINT :IN  controls the amount of output, as specified in the
  !                 SLATEC Guidelines.
  !
  !     IPASS  :OUT indicates whether the test passed or failed.
  !                 A value of one is good, indicating no failures.
  !
  ! *Description:
  !
  !   CQCBI is a quick check routine for the complex I Bessel function
  !    generated by subroutine CBESI.
  !
  !   CQCBI generates sequences crossing formula boundaries which
  !    are started by one formula and terminated in a region where
  !    another formula applies. The terminated value is checked by
  !    the formula appropriate to that region.
  !
  !***REFERENCES  Abramowitz, M. and Stegun, I. A., Handbook
  !                 of Mathematical Functions, Dover Publications,
  !                 New York, 1964.
  !               Amos, D. E., A Subroutine Package for Bessel
  !                 Functions of a Complex Argument and Nonnegative
  !                 Order, SAND85-1018, May, 1985.
  !***ROUTINES CALLED  CBESI, CBESK, CWRSK, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890831  Revised to meet new SLATEC standards
  !
  !***END PROLOGUE  CQCBI
  !
  !*Internal Notes:
  !   Machine constants are defined by functions I1MACH and R1MACH.
  !
  !   The parameter MQC can have values 1 (the default) for a faster,
  !   less definitive test or 2 for a slower, more definitive test.
  !
  !**End
  !
  !  Set test complexity parameter.
  !
  INTEGER MQC
  PARAMETER (MQC=1)
  !
  !  Declare arguments.
  !
  INTEGER Lun, Kprint, Ipass
  !
  !  Declare external functions.
  !
  INTEGER I1MACH
  REAL R1MACH
  EXTERNAL I1MACH, R1MACH
  !
  !  Declare local variables.
  !
  COMPLEX ck, cone, csgn, cw, cy, w, y, z, zn, zsc, zt, zz
  REAL aa, ab, aer, alim, arg, atol, aw, carg, ct, dig, elim, &
    eps, er, ertol, film, fnu, fnul, gnu, hpi, pi, r, rl, &
    rlt, rm, r1, r1m4, r1m5, r2, sarg, slak, st, t, tol, ts, &
    xx, yy
  INTEGER i, icase, ierr, il, inu, iprnt, ir, it, itl, k, kdo, &
    keps, kk, kode, k1, k2, lflg, mflg, n, nl, nzi, nzk, &
    nz1, nz2, n1
  DIMENSION aer(20), ck(2), kdo(20), keps(20), t(20), w(20), y(20)
  !
  !***FIRST EXECUTABLE STATEMENT  CQCBI
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99001)
    99001   FORMAT (' QUICK CHECK ROUTINE FOR THE I BESSEL FUNCTION FROM ','CBESI'/)
  ENDIF
  !-----------------------------------------------------------------------
  !     Set parameters related to machine constants.
  !     TOL is the approximate unit roundoff limited to 1.0E-18.
  !     ELIM is the approximate exponential over- and underflow limit.
  !     exp(-ELIM).lt.exp(-ALIM)=exp(-ELIM)/TOL    and
  !     exp(ELIM).gt.exp(ALIM)=exp(ELIM)*TOL       are intervals near
  !     underflow and overflow limits where scaled arithmetic is done.
  !     RL is the lower boundary of the asymptotic expansion for large Z.
  !     DIG = number of base 10 digits in TOL = 10**(-DIG).
  !     FNUL is the lower boundary of the asymptotic series for large FNU.
  !-----------------------------------------------------------------------
  r1m4 = R1MACH(4)
  tol = MAX(r1m4,1.0E-18)
  atol = 100.0E0*tol
  aa = -LOG10(r1m4)
  k1 = I1MACH(12)
  k2 = I1MACH(13)
  r1m5 = R1MACH(5)
  k = MIN(ABS(k1),ABS(k2))
  elim = 2.303E0*(k*r1m5-3.0E0)
  ab = aa*2.303E0
  alim = elim + MAX(-ab,-41.45E0)
  dig = MIN(aa,18.0E0)
  fnul = 10.0E0 + 6.0E0*(dig-3.0E0)
  rl = 1.2E0*dig + 3.0E0
  slak = 3.0E0 + 4.0E0*(-LOG10(tol)-7.0E0)/11.0E0
  slak = MAX(slak,3.0E0)
  ertol = tol*10.0E0**slak
  rm = 0.5E0*(alim+elim)
  rm = MIN(rm,200.0E0)
  rm = MAX(rm,rl+10.0E0)
  r2 = MIN(rm,fnul)
  r1 = 2.0E0*SQRT(fnul+1.0E0)
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99002)
    99002   FORMAT (' PARAMETERS'/5X,'TOL ',8X,'ELIM',8X,'ALIM',8X,'RL  ',8X,'FNUL',&
      8X,'DIG')
    WRITE (Lun,99003) tol, elim, alim, rl, fnul, dig
    99003   FORMAT (1X,6E12.4/)
  ENDIF
  !-----------------------------------------------------------------------
  !     Set other constants needed in the tests.
  !-----------------------------------------------------------------------
  zz = CMPLX(1.0E0/tol,0.0E0)
  cone = CMPLX(1.0E0,0.0E0)
  hpi = 2.0E0*ATAN(1.0E0)
  pi = hpi + hpi
  !-----------------------------------------------------------------------
  !     Generate angles for construction of complex Z to be used in tests.
  !-----------------------------------------------------------------------
  !     KDO(K), K = 1,IL  determines which of the IL angles in -PI to PI
  !     are used to compute values of Z.
  !       KDO(K) = 0  means that the index K will be used for one or two
  !                   values of Z, depending on the choice of KEPS(K)
  !              = 1  means that the index K and the corresponding angle
  !                   will be skipped
  !     KEPS(K), K = 1,IL determines which of the angles get incremented
  !     up and down to put values of Z in regions where different
  !     formulae are used.
  !       KEPS(K)  = 0  means that the angle will be used without change
  !                = 1  means that the angle will be incremented up and
  !                   down by EPS
  !     The angles to be used are stored in the T(I) array, I = 1,ITL.
  !-----------------------------------------------------------------------
  IF ( MQC/=2 ) THEN
    nl = 2
    il = 5
    DO i = 1, il
      keps(i) = 0
      kdo(i) = 0
    ENDDO
  ELSE
    nl = 4
    il = 13
    DO i = 1, il
      kdo(i) = 0
      keps(i) = 0
    ENDDO
    kdo(2) = 1
    kdo(6) = 1
    kdo(8) = 1
    kdo(12) = 1
    keps(3) = 1
    keps(4) = 1
    keps(5) = 1
    keps(9) = 1
    keps(10) = 1
    keps(11) = 1
  ENDIF
  i = 2
  eps = 0.01E0
  film = il - 1
  t(1) = -pi + eps
  DO k = 2, il
    IF ( kdo(k)==0 ) THEN
      t(i) = pi*(-il+2*k-1)/film
      IF ( keps(k)/=0 ) THEN
        ts = t(i)
        t(i) = ts - eps
        i = i + 1
        t(i) = ts + eps
      ELSE
        i = i + 1
      ENDIF
    ENDIF
  ENDDO
  itl = i - 1
  !-----------------------------------------------------------------------
  !     Test values of Z in -PI.lt.arg(Z).le.PI near formula boundaries.
  !-----------------------------------------------------------------------
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (' CHECKS ACROSS FORMULA BOUNDARIES')
  ENDIF
  lflg = 0
  DO icase = 1, 6
    DO kode = 1, 2
      DO n = 1, nl
        n1 = n + 2
        !-----------------------------------------------------------------------
        !     Set values for R = magnitude of Z and for FNU to test computation
        !     methods for the various regions of the (Z,FNU) plane.
        !-----------------------------------------------------------------------
        DO ir = 1, 3
          !------------ switch (icase)
          SELECT CASE (icase)
            CASE (2)
              r = (rl*(3-ir)+r2*(ir-1))/2.0E0
              gnu = SQRT(r+r) - 0.2E0 - (n-1)
              fnu = MAX(0.0E0,gnu)
            CASE (3)
              IF ( r2>=rm ) GOTO 100
              r = (r2*(3-ir)+rm*(ir-1))/2.0E0
              gnu = SQRT(r+r) - 0.2E0 - (n-1)
              fnu = MAX(0.0E0,gnu)
            CASE (4)
              IF ( r1>=rl ) GOTO 100
              r = (r1*(3-ir)+rl*(ir-1))/2.0E0
              fnu = fnul - 0.2E0 - (n-1)
            CASE (5)
              r = (rl*(3-ir)+r2*(ir-1))/2.0E0
              fnu = fnul - 0.2E0 - (n-1)
            CASE (6)
              IF ( r2>=rm ) GOTO 100
              r = (r2*(3-ir)+rm*(ir-1))/2.0E0
              fnu = fnul - 0.2E0 - (n-1)
            CASE DEFAULT
              r = (2.0E0*(3-ir)+rl*(ir-1))/2.0E0
              gnu = r*r/4.0E0 - 0.2E0 - (n-1)
              fnu = MAX(0.0E0,gnu)
          END SELECT
          !------------ end switch
          DO it = 1, itl
            ct = COS(t(it))
            st = SIN(t(it))
            IF ( ABS(ct)<atol ) ct = 0.0E0
            IF ( ABS(st)<atol ) st = 0.0E0
            z = CMPLX(r*ct,r*st)
            xx = REAL(z)
            yy = AIMAG(z)
            CALL CBESI(z,fnu,kode,n1,y,nz1,ierr)
            IF ( nz1==0 ) THEN
              !-----------------------------------------------------------------------
              !     Compare values from CBESI with values from CWRSK, an alternative
              !     method for calculating the complex Bessel I function.
              !-----------------------------------------------------------------------
              zn = z
              IF ( xx>=0.0E0 ) THEN
                CALL CWRSK(zn,fnu,kode,n,w,nz2,ck,tol,elim,alim)
                IF ( nz2/=0 ) CYCLE
              ELSE
                zn = -z
                inu = INT(fnu)
                arg = (fnu-inu)*pi
                IF ( yy<0.0E0 ) arg = -arg
                carg = COS(arg)
                sarg = SIN(arg)
                csgn = CMPLX(carg,sarg)
                IF ( MOD(inu,2)==1 ) csgn = -csgn
                CALL CWRSK(zn,fnu,kode,n,w,nz2,ck,tol,elim,alim)
                IF ( nz2/=0 ) CYCLE
                DO i = 1, n
                  w(i) = w(i)*csgn
                  csgn = -csgn
                ENDDO
              ENDIF
              mflg = 0
              DO i = 1, n
                ab = fnu + i - 1
                aa = MAX(2.0E0,ab)
                zt = w(i)
                IF ( CABS(zt)>1.0E0 ) THEN
                  zsc = CMPLX(tol,0.0E0)
                ELSE
                  zsc = zz
                  !------------------ ZZ = CMPLX(1.0/TOL,0.0)
                ENDIF
                cw = w(i)*zsc
                cy = y(i)*zsc
                er = CABS(cy-cw)
                aw = CABS(cw)
                IF ( aw==0.0E0 ) THEN
                  er = CABS(y(i))
                ELSEIF ( xx/=0.0E0 ) THEN
                  er = er/aw
                ELSEIF ( ABS(yy)<aa ) THEN
                  er = er/aw
                ELSE
                  er = CABS(w(i)-y(i))
                ENDIF
                aer(i) = er
                IF ( er>=ertol ) mflg = 1
              ENDDO
              !-----------------------------------------------------------------------
              !     Write failure reports for KPRINT.ge.2 and KPRINT.ge.3
              !-----------------------------------------------------------------------
              IF ( mflg/=0 ) THEN
                IF ( lflg==0 ) THEN
                  IF ( Kprint>=2 ) THEN
                    WRITE (Lun,99005) ertol
                    99005                   FORMAT (/' CASES WHICH UNDERFLOW OR VIOLATE THE ',&
                      'RELATIVE ERROR TEST'/' WITH ERTOL = ',E12.4/)
                    WRITE (Lun,99006)
                    99006                   FORMAT (' INPUT TO CBESI   Z, FNU, KODE, N')
                  ENDIF
                  IF ( Kprint>=3 ) THEN
                    WRITE (Lun,99007)
                    99007                   FORMAT (' ERROR TEST ON RESULTS FROM CBESI AND ',&
                      'CWRSK   AER(K)')
                    WRITE (Lun,99008)
                    99008                   FORMAT (' RESULTS FROM CBESI   NZ1, Y(KK)'/,&
                      ' RESULTS FROM CWRSK   NZ2, W(KK)')
                    WRITE (Lun,99009)
                    99009                   FORMAT (' TEST CASE INDICES   IT, IR, ICASE'/)
                  ENDIF
                ENDIF
                lflg = lflg + 1
                IF ( Kprint>=2 ) THEN
                  WRITE (Lun,99010) z, fnu, kode, n
                  99010                 FORMAT ('   INPUT:    Z=',2E12.4,4X,'FNU=',E12.4,4X,&
                    'KODE=',I3,4X,'N=',I3)
                ENDIF
                IF ( Kprint>=3 ) THEN
                  WRITE (Lun,99011) (aer(k),k=1,n)
                  99011                 FORMAT ('   ERROR:  AER(K)=',4E12.4)
                  kk = MAX(nz1,nz2) + 1
                  kk = MIN(n,kk)
                  WRITE (Lun,99012) nz1, y(kk), nz2, w(kk)
                  99012                 FORMAT (' RESULTS:  NZ1=',I3,4X,'Y(KK)=',2E12.4,/11X,&
                    'NZ2=',I3,4X,'W(KK)=',2E12.4)
                  WRITE (Lun,99013) it, ir, icase
                  99013                 FORMAT ('    CASE:   IT=',I3,4X,'IR=',I3,4X,'ICASE=',I3/)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    100  ENDDO
    IF ( Kprint>=2 ) THEN
      IF ( lflg==0 ) THEN
        WRITE (Lun,99019)
      ELSE
        WRITE (Lun,99014) lflg
        99014     FORMAT (' ***',I5,' FAILURE(S) FOR CBESI CHECKS NEAR FORMULA ',&
          'BOUNDARIES')
      ENDIF
    ENDIF
    !
    !
    iprnt = 0
    IF ( MQC/=1 ) THEN
      !-----------------------------------------------------------------------
      !     Checks near underflow limits on series(I=1) and uniform
      !     asymptotic expansion(I=2)
      !     Compare 1/Z with I(Z,FNU)*K(Z,FNU+1) + I(Z,FNU+1)*K(Z,FNU) and
      !     report cases for which the relative error is greater than ERTOL.
      !-----------------------------------------------------------------------
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99015)
        99015     FORMAT (/' CHECKS NEAR UNDERFLOW AND OVERFLOW LIMITS'/)
      ENDIF
      z = (1.4E0,1.4E0)
      kode = 1
      n = 20
      DO i = 1, 2
        fnu = 10.2E0
        DO
          !-----------------------------------------------------------------------
          !       Adjust FNU by repeating until 0.lt.NZI.lt.10
          !-----------------------------------------------------------------------
          CALL CBESI(z,fnu,kode,n,y,nzi,ierr)
          IF ( nzi==0 ) THEN
            fnu = fnu + 5.0E0
            CYCLE
          ELSEIF ( nzi>=10 ) THEN
            fnu = fnu - 1.0E0
            CYCLE
          ENDIF
          !------ End repeat
          CALL CBESK(z,fnu,kode,2,w,nzk,ierr)
          zt = cone/z
          cy = w(1)*y(2)
          cw = w(2)*y(1)
          cw = cw + cy - zt
          er = ABS(cw)/ABS(zt)
          !-----------------------------------------------------------------------
          !     Write failure reports for KPRINT.ge.2 and KPRINT.ge.3
          !-----------------------------------------------------------------------
          IF ( er>=ertol ) THEN
            IF ( iprnt==0 ) THEN
              IF ( Kprint>=2 ) WRITE (Lun,99020)
              IF ( Kprint>=3 ) WRITE (Lun,99021)
            ENDIF
            iprnt = 1
            IF ( Kprint>=2 ) WRITE (Lun,99022) z, fnu, kode, n
            IF ( Kprint>=3 ) THEN
              WRITE (Lun,99023) zt, cw + cy
              WRITE (Lun,99024) er
            ENDIF
          ENDIF
          rlt = rl + rl
          z = CMPLX(rlt,0.0E0)
          EXIT
        ENDDO
      ENDDO
      !-----------------------------------------------------------------------
      !     Check near overflow limits
      !     Compare 1/Z with I(Z,FNU)*K(Z,FNU+1) + I(Z,FNU+1)*K(Z,FNU) and
      !     report cases for which the relative error is greater than ERTOL.
      !-----------------------------------------------------------------------
      z = CMPLX(elim,0.0E0)
      fnu = 0.0E0
      DO
        !-----------------------------------------------------------------------
        !     Adjust FNU by repeating until NZK.lt.10
        !     N = 20 set before DO 280 loop
        !-----------------------------------------------------------------------
        CALL CBESK(z,fnu,kode,n,y,nzk,ierr)
        IF ( nzk>=10 ) THEN
          IF ( nzk==n ) THEN
            fnu = fnu + 3.0E0
          ELSE
            fnu = fnu + 2.0E0
          ENDIF
          CYCLE
        ENDIF
        !---- End repeat
        gnu = fnu + (n-2)
        CALL CBESI(z,gnu,kode,2,w,nzi,ierr)
        zt = cone/z
        cy = y(n-1)*w(2)
        cw = y(n)*w(1)
        cw = cw + cy - zt
        er = ABS(cw)/ABS(zt)
        IF ( er>=ertol ) THEN
          IF ( iprnt==0 ) THEN
            IF ( Kprint>=2 ) WRITE (Lun,99020)
            IF ( Kprint>=3 ) WRITE (Lun,99021)
          ENDIF
          iprnt = 1
          IF ( Kprint>=2 ) WRITE (Lun,99022) z, fnu, kode, n
          IF ( Kprint>=3 ) THEN
            WRITE (Lun,99023) zt, cw + cy
            WRITE (Lun,99024) er
          ENDIF
        ENDIF
        IF ( Kprint>=2 ) THEN
          IF ( iprnt==0 ) THEN
            WRITE (Lun,99019)
            ! 99986   FORMAT (' QUICK CHECKS OK')
          ELSE
            WRITE (Lun,99016)
            99016         FORMAT (' ***',5X,'FAILURE(S) FOR CBESI NEAR UNDERFLOW AND ',&
              'OVERFLOW LIMITS')
          ENDIF
        ENDIF
        EXIT
      ENDDO
    ENDIF
    Ipass = 0
    IF ( iprnt==0.AND.lflg==0 ) Ipass = 1
    IF ( Ipass==1.AND.Kprint>=2 ) THEN
      WRITE (Lun,99017)
      99017   FORMAT (/' ****** CBESI  PASSED ALL TESTS  ******'/)
    ENDIF
    IF ( Ipass==0.AND.Kprint>=1 ) THEN
      WRITE (Lun,99018)
      99018   FORMAT (/' ****** CBESI  FAILED SOME TESTS ******'/)
    ENDIF
    99019 FORMAT (' QUICK CHECKS OK')
    99020 FORMAT (' INPUT TO CBESI   Z, FNU, KODE, N')
    99021 FORMAT (' COMPARE 1/Z WITH WRONSKIAN(CBESI(Z,FNU),','CBESK(Z,FNU))'/)
    99022 FORMAT (' INPUT: Z=',2E12.4,3X,'FNU=',E12.4,3X,'KODE=',I3,3X,'N=',I3)
    99023 FORMAT (' RESULTS:',15X,'1/Z=',2E12.4/10X,'WRON(CBESI,CBESK)=',2E12.4)
    99024 FORMAT (' RELATIVE ERROR:',9X,'ER=',E12.4/)
END SUBROUTINE CQCBI
