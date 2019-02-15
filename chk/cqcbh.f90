!DECK CQCBH
SUBROUTINE CQCBH(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CQCBH
  !***SUBSIDIARY
  !***PURPOSE  Quick check for SLATEC subroutine
  !            CBESH
  !***LIBRARY   SLATEC
  !***CATEGORY  C10A4
  !***TYPE      COMPLEX (CQCBH-C, ZQCBH-Z)
  !***KEYWORDS  QUICK CHECK, CBESH
  !***AUTHOR  Amos, Don, (SNL)
  !           Goudy, Sue, (SNL)
  !           Walton, Lee, (SNL)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL CQCBH (LUN, KPRINT, IPASS)
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
  !   CQCBH is a quick check routine for the complex H Bessel functions
  !    generated by subroutine CBESH.
  !
  !   CQCBH generates sequences of H Bessel functions for kinds 1 and 2
  !    from CBESH and checks them against the Wronskian evaluation
  !    in the (Z,FNU) space.
  !
  !***REFERENCES  Abramowitz, M. and Stegun, I. A., Handbook
  !                 of Mathematical Functions, Dover Publications,
  !                 New York, 1964.
  !               Amos, D. E., A Subroutine Package for Bessel
  !                 Functions of a Complex Argument and Nonnegative
  !                 Order, SAND85-1018, May, 1985.
  !***ROUTINES CALLED  CBESH, CUOIK, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   890831  Revised to meet new SLATEC standards
  !
  !***END PROLOGUE  CQCBH
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
  COMPLEX cv, cw, cy, w, y, z, zn
  REAL aa, ab, acw, acy, aer, alim, atol, av, aw, ay, az, ct, &
    dig, elim, eps, er, ertol, film, fnu, fnul, fpi, hpi, pi, &
    r, rfpi, rl, rm, r1m4, r1m5, r2, slak, st, t, tol, ts, &
    xnu
  INTEGER i, icase, ierr, il, ir, irb, it, itl, k, kdo, keps, &
    kk, kode, k1, k2, lflg, mflg, n, nl, nu, nul, nz1, &
    nz2, n1
  DIMENSION aer(20), kdo(20), keps(20), t(20), w(20), xnu(20), y(20)
  !
  !***FIRST EXECUTABLE STATEMENT  CQCBH
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99001)
    99001   FORMAT (' QUICK CHECK ROUTINE FOR THE H BESSEL FUNCTIONS FROM ',&
      'CBESH'/)
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
  r2 = MIN(fnul,rm)
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99002)
    99002   FORMAT (' PARAMETERS'/5X,'TOL ',8X,'ELIM',8X,'ALIM',8X,'RL  ',8X,'FNUL',&
      8X,'DIG')
    WRITE (Lun,99003) tol, elim, alim, rl, fnul, dig
    99003   FORMAT (6E12.4/)
  ENDIF
  !-----------------------------------------------------------------------
  !     Set other constants needed in the tests.
  !-----------------------------------------------------------------------
  fpi = ATAN(1.0E0)
  hpi = fpi + fpi
  pi = hpi + hpi
  rfpi = 1.0E0/fpi
  zn = CMPLX(0.0E0,-rfpi)
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
    nul = 5
    xnu(1) = 0.0E0
    xnu(2) = 1.0E0
    xnu(3) = 2.0E0
    xnu(4) = 0.5E0*fnul
    xnu(5) = fnul + 1.1E0
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
    nul = 6
    xnu(1) = 0.0E0
    xnu(2) = 0.6E0
    xnu(3) = 1.3E0
    xnu(4) = 2.0E0
    xnu(5) = 0.5E0*fnul
    xnu(6) = fnul + 1.1E0
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
  !     Test values of Z in -PI.lt.arg(Z).le.PI.
  !-----------------------------------------------------------------------
  IF ( Kprint>=2 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (' CHECKS IN THE (Z,FNU) SPACE'/)
  ENDIF
  lflg = 0
  DO kode = 1, 2
    DO n = 1, nl
      n1 = n + 1
      DO nu = 1, nul
        fnu = xnu(nu)
        DO icase = 1, 3
          irb = MIN(icase,2)
          DO ir = irb, 3
            !-------------- switch (icase)
            SELECT CASE (icase)
              CASE (2)
                r = (2.0E0*(3-ir)+r2*(ir-1))/2.0E0
              CASE (3)
                IF ( r2>=rm ) EXIT
                r = (r2*(3-ir)+rm*(ir-1))/2.0E0
              CASE DEFAULT
                r = (eps*(3-ir)+2.0E0*(ir-1))/2.0E0
            END SELECT
            !-------------- end switch
            DO it = 1, itl
              ct = COS(t(it))
              st = SIN(t(it))
              IF ( ABS(ct)<atol ) ct = 0.0E0
              IF ( ABS(st)<atol ) st = 0.0E0
              z = CMPLX(r*ct,r*st)
              IF ( fnu>=2.0E0 ) THEN
                !------------------ Check for possible overflow condition
                cv = z*(0.0E0,1.0E0)
                CALL CUOIK(cv,fnu,kode,2,n1,w,nz2,tol,elim,alim)
                !------------------ Overflow detected? - skip test for this case
                IF ( nz2==(-1) ) CYCLE
                cv = -cv
                CALL CUOIK(cv,fnu,kode,2,n1,w,nz2,tol,elim,alim)
                !------------------ Overflow detected? - skip test for this case
                IF ( nz2==(-1) ) CYCLE
              ENDIF
              !---------------- No overflow - calculate H1(Z,FNU) and H2(Z,FNU)
              CALL CBESH(z,fnu,kode,1,n1,y,nz1,ierr)
              !---------------- Underflow? - skip test for this case
              IF ( nz1==0 ) THEN
                CALL CBESH(z,fnu,kode,2,n1,w,nz2,ierr)
                !---------------- Underflow? - skip test for this case
                IF ( nz2==0 ) THEN
                  !-----------------------------------------------------------------------
                  !     Compare ZN/Z with the Wronskian of H1(Z,FNU) and H2(Z,FNU).
                  !     ZN = -4i/PI
                  !-----------------------------------------------------------------------
                  cv = zn/z
                  mflg = 0
                  DO i = 1, n
                    !-----------------------------------------------------------------------
                    !     Error relative to maximum term
                    !-----------------------------------------------------------------------
                    aw = ABS(w(i+1))
                    ay = ABS(y(i))
                    az = LOG(aw) + LOG(ay)
                    az = ABS(az)
                    IF ( az<=alim ) THEN
                      !-------------------- No scaling problem - do error analysis
                      av = ABS(cv)
                      cw = w(i)*y(i+1)
                      cy = w(i+1)*y(i)
                      cy = cw - cy - cv
                      acy = aw*ay
                      acw = ABS(w(i))*ABS(y(i+1))
                      av = MAX(acw,acy,av)
                      er = ABS(cy)/av
                      aer(i) = er
                      IF ( er>ertol ) mflg = 1
                    ENDIF
                  ENDDO
                  IF ( mflg/=0 ) THEN
                    IF ( lflg==0 ) THEN
                      IF ( Kprint>=2 ) THEN
                        WRITE (Lun,99005) ertol
                        99005                       FORMAT (' CASES WHICH VIOLATE THE RELATIVE ',&
                          'ERROR TEST WITH ERTOL = ',E12.4/)
                        WRITE (Lun,99006)
                        99006                       FORMAT (' INPUT TO CBESH   Z, FNU, KODE, N')
                      ENDIF
                      IF ( Kprint>=3 ) THEN
                        WRITE (Lun,99007)
                        99007                       FORMAT (' COMPARE -4i/(PI*Z) WITH WRONSKIAN OF',&
                          ' H1(Z,FNU) AND H2(Z,FNU)')
                        WRITE (Lun,99008)
                        99008                       FORMAT (' RESULTS FROM CBESH FOR FUNCTION H1'/&
                          '         AND FUNCTION H2')
                        WRITE (Lun,99009)
                        99009                       FORMAT (' TEST CASE INDICES'/)
                      ENDIF
                    ENDIF
                    lflg = lflg + 1
                    IF ( Kprint>=2 ) THEN
                      WRITE (Lun,99010) z, fnu, kode, n
                      99010                     FORMAT ('   INPUT:    Z=',2E12.4,4X,'FNU=',E12.4,4X,&
                        'KODE=',I3,4X,'N=',I3)
                    ENDIF
                    IF ( Kprint>=3 ) THEN
                      WRITE (Lun,99011) (aer(k),k=1,n)
                      99011                     FORMAT ('   ERROR:   AER(K)=',4E12.4)
                      kk = MAX(nz1,nz2) + 1
                      kk = MIN(n,kk)
                      WRITE (Lun,99012) nz1, y(kk), nz2, w(kk)
                      99012                     FORMAT (' RESULTS:  NZ1=',I3,4X,'Y(KK)=',2E12.4/11X,&
                        'NZ2=',I3,4X,'W(KK)=',2E12.4)
                      WRITE (Lun,99013) it, ir, icase
                      99013                     FORMAT ('    CASE:   IT=',I3,4X,'IR=',I3,4X,'ICASE=',&
                        I3/)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  IF ( Kprint>=2 ) THEN
    IF ( lflg==0 ) THEN
      WRITE (Lun,99014)
      99014     FORMAT (' QUICK CHECKS OK')
    ELSE
      WRITE (Lun,99015) lflg
      99015     FORMAT (' ***',I5,' FAILURE(S) FOR CBESH IN THE (Z,FNU)',' PLANE')
    ENDIF
  ENDIF
  Ipass = 0
  IF ( lflg==0 ) Ipass = 1
  IF ( Ipass==1.AND.Kprint>=2 ) THEN
    WRITE (Lun,99016)
    99016   FORMAT (/' ****** CBESH  PASSED ALL TESTS  ******'/)
  ENDIF
  IF ( Ipass==0.AND.Kprint>=1 ) THEN
    WRITE (Lun,99017)
    99017   FORMAT (/' ****** CBESH  FAILED SOME TESTS ******'/)
  ENDIF
END SUBROUTINE CQCBH
