MODULE TEST51_MOD
  IMPLICIT NONE

CONTAINS
  !** FFTQX
  SUBROUTINE FFTQX(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for the NCAR FFT routines.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Swarztrauber, P. N., (NCAR)
    !***
    ! **Description:**
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !                       VERSION 4  APRIL 1985
    !
    !                         A TEST DRIVER FOR
    !          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER
    !           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES
    !
    !                              BY
    !
    !                       PAUL N SWARZTRAUBER
    !
    !    NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER, COLORADO 80307
    !
    !        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER
    !     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
    !     CERTAIN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
    !
    !     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
    !     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
    !     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
    !
    !     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB
    !     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
    !     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
    !
    !     7.   SINTI     INITIALIZE SINT
    !     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
    !
    !     9.   COSTI     INITIALIZE COST
    !     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
    !
    !     11.  SINQI     INITIALIZE SINQF AND SINQB
    !     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
    !     13.  SINQB     UNNORMALIZED INVERSE OF SINQF
    !
    !     14.  COSQI     INITIALIZE COSQF AND COSQB
    !     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
    !     16.  COSQB     UNNORMALIZED INVERSE OF COSQF
    !
    !     17.  CFFTI     INITIALIZE CFFTF AND CFFTB
    !     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
    !     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
    !
    !***
    ! **Routines called:**  CFFTB, CFFTF, CFFTI, COSQB, COSQF, COSQI, COST,
    !                    COSTI, EZFFTB, EZFFTF, EZFFTI, PIMACH, R1MACH,
    !                    RFFTB, RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT,
    !                    SINTI

    !* REVISION HISTORY  (YYMMDD)
    !   790601  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920211  Code cleaned up, an error in printing an error message fixed
    !           and comments on PASS/FAIL of individual tests added.  (WRB)
    !   920618  Code upgraded to "Version 4".  (BKS, WRB)
    !   930315  Modified RFFT* tests to compute "slow-transform" in DOUBLE
    !           PRECISION.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(8) :: arg, arg1, arg2, dt, pi, sum, sum1, sum2
    REAL azero, azeroh, cf, cosqbt, cosqfb, cosqft, costfb, costt, &
      dcfb, dcfftb, dcfftf, dezb1, dezf1, dezfb, errmax, rftb, &
      rftf, rftfb, sign, sinqbt, sinqfb, sinqft, sintfb, sintt, sqrt2, tpi
    INTEGER i, j, k, modn, n, nm1, nns, np1, ns2, ns2m, nz
    !     .. Local Arrays ..
    COMPLEX cx(200), cy(200)
    REAL a(100), ah(100), b(100), bh(100), w(2000), x(200), xh(200), y(200)
    !     .. External Functions ..
    REAL, EXTERNAL :: R1MACH
    !     .. External Subroutines ..
    EXTERNAL :: CFFTB, CFFTF, CFFTI, COSQB, COSQF, COSQI, COST, COSTI, EZFFTB, &
      EZFFTF, EZFFTI, RFFTB, RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT, SINTI
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, CABS, CMPLX, COS, MAX, MOD, SIN, SQRT
    !     .. Data statements ..
    INTEGER, PARAMETER :: nd(7) = [ 120, 54, 49, 32, 4, 3, 2 ]
    !* FIRST EXECUTABLE STATEMENT  FFTQX
    sqrt2 = SQRT(2.0)
    errmax = 2.0*SQRT(R1MACH(4))
    nns = 7
    pi = 4.0D0*ATAN(1.0D0)
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1'/' FFT QUICK CHECK')
    Ipass = 1
    DO nz = 1, nns
      n = nd(nz)
      IF ( Kprint>=2 ) WRITE (Lun,99002) n
      99002 FORMAT (/' Test FFT routines with a sequence of length ',I3)
      modn = MOD(n,2)
      np1 = n + 1
      nm1 = n - 1
      DO j = 1, np1
        x(j) = SIN(j*sqrt2)
        y(j) = x(j)
        xh(j) = x(j)
      ENDDO
      !
      !       Test Subroutines RFFTI, RFFTF and RFFTB
      !
      CALL RFFTI(n,w)
      dt = (pi+pi)/n
      ns2 = (n+1)/2
      IF ( ns2>=2 ) THEN
        DO k = 2, ns2
          sum1 = 0.0D0
          sum2 = 0.0D0
          arg = (k-1)*dt
          DO i = 1, n
            arg1 = (i-1)*arg
            sum1 = sum1 + x(i)*COS(arg1)
            sum2 = sum2 + x(i)*SIN(arg1)
          ENDDO
          y(2*k-2) = REAL( sum1, 4 )
          y(2*k-1) = REAL( -sum2, 4 )
        ENDDO
      ENDIF
      sum1 = 0.0D0
      sum2 = 0.0D0
      DO i = 1, nm1, 2
        sum1 = sum1 + x(i)
        sum2 = sum2 + x(i+1)
      ENDDO
      IF ( modn==1 ) sum1 = sum1 + x(n)
      y(1) = REAL( sum1 + sum2, 4 )
      IF ( modn==0 ) y(n) = REAL( sum1 - sum2, 4 )
      CALL RFFTF(n,x,w)
      rftf = 0.0
      DO i = 1, n
        rftf = MAX(rftf,ABS(x(i)-y(i)))
        x(i) = xh(i)
      ENDDO
      rftf = rftf/n
      IF ( rftf<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99003)
        99003 FORMAT (' Test of RFFTF PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99004)
        99004 FORMAT (' Test of RFFTF FAILED')
      ENDIF
      sign = 1.0
      DO i = 1, n
        sum = 0.5D0*x(1)
        arg = (i-1)*dt
        IF ( ns2>=2 ) THEN
          DO k = 2, ns2
            arg1 = (k-1)*arg
            sum = sum + x(2*k-2)*COS(arg1) - x(2*k-1)*SIN(arg1)
          ENDDO
        ENDIF
        IF ( modn==0 ) sum = sum + 0.5D0*sign*x(n)
        y(i) = REAL( sum + sum, 4 )
        sign = -sign
      ENDDO
      CALL RFFTB(n,x,w)
      rftb = 0.0
      DO i = 1, n
        rftb = MAX(rftb,ABS(x(i)-y(i)))
        x(i) = xh(i)
        y(i) = xh(i)
      ENDDO
      rftb = rftb/n
      IF ( rftb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99005)
        99005 FORMAT (' Test of RFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99006)
        99006 FORMAT (' Test of RFFTB FAILED')
      ENDIF
      !
      CALL RFFTB(n,y,w)
      CALL RFFTF(n,y,w)
      cf = 1.0/n
      rftfb = 0.0
      DO i = 1, n
        rftfb = MAX(rftfb,ABS(cf*y(i)-x(i)))
      ENDDO
      IF ( rftfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99007)
        99007 FORMAT (' Test of RFFTF and RFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99008)
        99008 FORMAT (' Test of RFFTF and RFFTB FAILED')
      ENDIF
      !
      !       Test Subroutines SINTI and SINT
      !
      dt = pi/n
      DO i = 1, nm1
        x(i) = xh(i)
      ENDDO
      DO i = 1, nm1
        y(i) = 0.0
        arg1 = i*dt
        DO k = 1, nm1
          y(i) = y(i) + x(k)*REAL( SIN((k)*arg1), 4 )
        ENDDO
        y(i) = y(i) + y(i)
      ENDDO
      CALL SINTI(nm1,w)
      CALL SINT(nm1,x,w)
      cf = 0.5/n
      sintt = 0.0
      DO i = 1, nm1
        sintt = MAX(sintt,ABS(x(i)-y(i)))
        x(i) = xh(i)
        y(i) = x(i)
      ENDDO
      sintt = cf*sintt
      IF ( sintt<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99009)
        99009 FORMAT (' First test of SINT PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99010)
        99010 FORMAT (' First test of SINT FAILED')
      ENDIF
      CALL SINT(nm1,x,w)
      CALL SINT(nm1,x,w)
      sintfb = 0.0
      DO i = 1, nm1
        sintfb = MAX(sintfb,ABS(cf*x(i)-y(i)))
      ENDDO
      IF ( sintfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99011)
        99011 FORMAT (' Second test of SINT PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99012)
        99012 FORMAT (' Second test of SINT FAILED')
      ENDIF
      !
      !       Test Subroutines COSTI and COST
      !
      DO i = 1, np1
        x(i) = xh(i)
      ENDDO
      sign = 1.0
      DO i = 1, np1
        y(i) = 0.5*(x(1)+sign*x(n+1))
        arg = (i-1)*dt
        DO k = 2, n
          y(i) = y(i) + x(k)*REAL( COS((k-1)*arg), 4 )
        ENDDO
        y(i) = y(i) + y(i)
        sign = -sign
      ENDDO
      CALL COSTI(np1,w)
      CALL COST(np1,x,w)
      costt = 0.0
      DO i = 1, np1
        costt = MAX(costt,ABS(x(i)-y(i)))
        x(i) = xh(i)
        y(i) = xh(i)
      ENDDO
      costt = cf*costt
      IF ( costt<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99013)
        99013 FORMAT (' First test of COST PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99014)
        99014 FORMAT (' First test of COST FAILED')
      ENDIF
      !
      CALL COST(np1,x,w)
      CALL COST(np1,x,w)
      costfb = 0.0
      DO i = 1, np1
        costfb = MAX(costfb,ABS(cf*x(i)-y(i)))
      ENDDO
      IF ( costfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99015)
        99015 FORMAT (' Second test of COST PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99016)
        99016 FORMAT (' Second test of COST FAILED')
      ENDIF
      !
      !       Test Subroutines SINQI, SINQF and SINQB
      !
      cf = 0.25/n
      DO i = 1, n
        y(i) = xh(i)
      ENDDO
      dt = pi/(n+n)
      DO i = 1, n
        x(i) = 0.0
        arg = i*dt
        DO k = 1, n
          x(i) = x(i) + y(k)*REAL( SIN((k+k-1)*arg), 4 )
        ENDDO
        x(i) = 4.0*x(i)
      ENDDO
      CALL SINQI(n,w)
      CALL SINQB(n,y,w)
      sinqbt = 0.0
      DO i = 1, n
        sinqbt = MAX(sinqbt,ABS(y(i)-x(i)))
        x(i) = xh(i)
      ENDDO
      sinqbt = cf*sinqbt
      IF ( sinqbt<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99017)
        99017 FORMAT (' Test of SINQB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99018)
        99018 FORMAT (' Test of SINQB FAILED')
      ENDIF
      !
      sign = 1.0
      DO i = 1, n
        arg = (i+i-1)*dt
        y(i) = 0.5*sign*x(n)
        DO k = 1, nm1
          y(i) = y(i) + x(k)*REAL( SIN((k)*arg), 4 )
        ENDDO
        y(i) = y(i) + y(i)
        sign = -sign
      ENDDO
      CALL SINQF(n,x,w)
      sinqft = 0.0
      DO i = 1, n
        sinqft = MAX(sinqft,ABS(x(i)-y(i)))
        y(i) = xh(i)
        x(i) = xh(i)
      ENDDO
      IF ( sinqft<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99019)
        99019 FORMAT (' Test of SINQF PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99020)
        99020 FORMAT (' Test of SINQF FAILED')
      ENDIF
      !
      CALL SINQF(n,y,w)
      CALL SINQB(n,y,w)
      sinqfb = 0.0
      DO i = 1, n
        sinqfb = MAX(sinqfb,ABS(cf*y(i)-x(i)))
      ENDDO
      IF ( sinqfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99021)
        99021 FORMAT (' Test of SINQF and SINQB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99022)
        99022 FORMAT (' Test of SINQF and SINQB FAILED')
      ENDIF
      !
      !       Test Subroutines COSQI, COSQF and COSQB
      !
      DO i = 1, n
        y(i) = xh(i)
      ENDDO
      DO i = 1, n
        x(i) = 0.0
        arg = (i-1)*dt
        DO k = 1, n
          x(i) = x(i) + y(k)*REAL( COS((k+k-1)*arg), 4 )
        ENDDO
        x(i) = 4.0*x(i)
      ENDDO
      CALL COSQI(n,w)
      CALL COSQB(n,y,w)
      cosqbt = 0.0
      DO i = 1, n
        cosqbt = MAX(cosqbt,ABS(x(i)-y(i)))
        x(i) = xh(i)
      ENDDO
      cosqbt = cf*cosqbt
      IF ( cosqbt<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99023)
        99023 FORMAT (' Test of COSQB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99024)
        99024 FORMAT (' Test of COSQB FAILED')
      ENDIF
      !
      DO i = 1, n
        y(i) = 0.5*x(1)
        arg = (i+i-1)*dt
        DO k = 2, n
          y(i) = y(i) + x(k)*REAL( COS((k-1)*arg), 4 )
        ENDDO
        y(i) = y(i) + y(i)
      ENDDO
      CALL COSQF(n,x,w)
      cosqft = 0.0
      DO i = 1, n
        cosqft = MAX(cosqft,ABS(y(i)-x(i)))
        x(i) = xh(i)
        y(i) = xh(i)
      ENDDO
      cosqft = cf*cosqft
      IF ( cosqft<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99025)
        99025 FORMAT (' Test of COSQF PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99026)
        99026 FORMAT (' Test of COSQF FAILED')
      ENDIF
      !
      CALL COSQB(n,x,w)
      CALL COSQF(n,x,w)
      cosqfb = 0.0
      DO i = 1, n
        cosqfb = MAX(cosqfb,ABS(cf*x(i)-y(i)))
      ENDDO
      IF ( cosqfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99027)
        99027 FORMAT (' Test of COSQF and COSQB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99028)
        99028 FORMAT (' Test of COSQF and COSQB FAILED')
      ENDIF
      !
      !       Test Subroutines EZFFTI, EZFFTF and EZFFTB
      !
      CALL EZFFTI(n,w)
      DO i = 1, n
        x(i) = xh(i)
      ENDDO
      tpi = REAL( 2.0*pi, 4 )
      dt = tpi/n
      ns2 = (n+1)/2
      cf = 2.0/n
      ns2m = ns2 - 1
      IF ( ns2m>0 ) THEN
        DO k = 1, ns2m
          sum1 = 0.0D0
          sum2 = 0.0D0
          arg = k*dt
          DO i = 1, n
            arg1 = (i-1)*arg
            sum1 = sum1 + x(i)*COS(arg1)
            sum2 = sum2 + x(i)*SIN(arg1)
          ENDDO
          a(k) = REAL( cf*sum1, 4 )
          b(k) = REAL( cf*sum2, 4 )
        ENDDO
      ENDIF
      nm1 = n - 1
      sum1 = 0.0D0
      sum2 = 0.0D0
      DO i = 1, nm1, 2
        sum1 = sum1 + x(i)
        sum2 = sum2 + x(i+1)
      ENDDO
      IF ( modn==1 ) sum1 = sum1 + x(n)
      azero = REAL( 0.5*cf*(sum1+sum2), 4 )
      IF ( modn==0 ) a(ns2) = REAL( 0.5*cf*(sum1-sum2), 4 )
      CALL EZFFTF(n,x,azeroh,ah,bh,w)
      dezf1 = ABS(azeroh-azero)
      IF ( modn==0 ) dezf1 = MAX(dezf1,ABS(a(ns2)-ah(ns2)))
      IF ( ns2m>0 ) THEN
        DO i = 1, ns2m
          dezf1 = MAX(dezf1,ABS(ah(i)-a(i)),ABS(bh(i)-b(i)))
        ENDDO
        IF ( dezf1<=errmax ) THEN
          IF ( Kprint>=3 ) WRITE (Lun,99029)
          99029 FORMAT (' Test of EZFFTF PASSED')
        ELSE
          Ipass = 0
          IF ( Kprint>=2 ) WRITE (Lun,99030)
          99030 FORMAT (' Test of EZFFTF FAILED')
        ENDIF
      ENDIF
      !
      ns2 = n/2
      IF ( modn==0 ) b(ns2) = 0.0
      DO i = 1, n
        sum = azero
        arg1 = (i-1)*dt
        DO k = 1, ns2
          arg2 = k*arg1
          sum = sum + a(k)*COS(arg2) + b(k)*SIN(arg2)
        ENDDO
        x(i) = REAL( sum, 4 )
      ENDDO
      CALL EZFFTB(n,y,azero,a,b,w)
      dezb1 = 0.0
      DO i = 1, n
        dezb1 = MAX(dezb1,ABS(x(i)-y(i)))
        x(i) = xh(i)
      ENDDO
      IF ( dezb1<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99031)
        99031 FORMAT (' Test of EZFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99032)
        99032 FORMAT (' Test of EZFFTB FAILED')
      ENDIF
      !
      CALL EZFFTF(n,x,azero,a,b,w)
      CALL EZFFTB(n,y,azero,a,b,w)
      dezfb = 0.0
      DO i = 1, n
        dezfb = MAX(dezfb,ABS(x(i)-y(i)))
      ENDDO
      IF ( dezfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99033)
        99033 FORMAT (' Test of EZFFTF and EZFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99034)
        99034 FORMAT (' Test of EZFFTF and EZFFTB FAILED')
      ENDIF
      !
      !       Test Subroutines CFFTI, CFFTF and CFFTB
      !
      DO i = 1, n
        cx(i) = CMPLX(COS(sqrt2*i),SIN(sqrt2*(i*i)))
      ENDDO
      dt = (pi+pi)/n
      DO i = 1, n
        arg1 = -(i-1)*dt
        cy(i) = (0.0,0.0)
        DO k = 1, n
          arg2 = (k-1)*arg1
          cy(i) = cy(i) + CMPLX(COS(arg2),SIN(arg2),4)*cx(k)
        ENDDO
      ENDDO
      CALL CFFTI(n,w)
      CALL CFFTF(n,cx,w)
      dcfftf = 0.0
      DO i = 1, n
        dcfftf = MAX(dcfftf,CABS(cx(i)-cy(i)))
        cx(i) = cx(i)/n
      ENDDO
      dcfftf = dcfftf/n
      IF ( dcfftf<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99035)
        99035 FORMAT (' Test of CFFTF PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99036)
        99036 FORMAT (' Test of CFFTF FAILED')
      ENDIF
      !
      DO i = 1, n
        arg1 = (i-1)*dt
        cy(i) = (0.0,0.0)
        DO k = 1, n
          arg2 = (k-1)*arg1
          cy(i) = cy(i) + CMPLX(COS(arg2),SIN(arg2),4)*cx(k)
        ENDDO
      ENDDO
      CALL CFFTB(n,cx,w)
      dcfftb = 0.0
      DO i = 1, n
        dcfftb = MAX(dcfftb,CABS(cx(i)-cy(i)))
        cx(i) = cy(i)
      ENDDO
      IF ( dcfftb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99037)
        99037 FORMAT (' Test of CFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99038)
        99038 FORMAT (' Test of CFFTB FAILED')
      ENDIF
      !
      cf = 1.0/n
      CALL CFFTF(n,cx,w)
      CALL CFFTB(n,cx,w)
      dcfb = 0.0
      DO i = 1, n
        dcfb = MAX(dcfb,CABS(cf*cx(i)-cy(i)))
      ENDDO
      IF ( dcfb<=errmax ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99039)
        99039 FORMAT (' Test of CFFTF and CFFTB PASSED')
      ELSE
        Ipass = 0
        IF ( Kprint>=2 ) WRITE (Lun,99040)
        99040 FORMAT (' Test of CFFTF and CFFTB FAILED')
      ENDIF
      IF ( Kprint>=3 ) THEN
        WRITE (Lun,99041) n, rftf, rftb, rftfb, sintt, sintfb, costt, &
          costfb, sinqft, sinqbt, sinqfb, cosqft, &
          cosqbt, cosqfb, dezf1, dezb1, dezfb, dcfftf, dcfftb, dcfb
        99041 FORMAT ('0N',I5,'  RFFTF  ',E10.3,'  RFFTB  ',E10.3,'  RFFTFB ',E10.3/7X,&
          '  SINT   ',E10.3,'  SINTFB ',E10.3/7X,'  COST   ',E10.3,&
          '  COSTFB ',E10.3/7X,'  SINQF  ',E10.3,'  SINQB  ',E10.3,&
          '  SINQFB ',E10.3/7X,'  COSQF  ',E10.3,'  COSQB  ',E10.3,&
          '  COSQFB ',E10.3/7X,'  DEZF1  ',E10.3,'  DEZB1  ',E10.3,&
          '  DEZFB  ',E10.3/7X,'  CFFTF  ',E10.3,'  CFFTB  ',E10.3,&
          '  CFFTFB ',E10.3)
      ENDIF
    ENDDO
    IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99042)
    99042 FORMAT (/' ***********FFT ROUTINES PASSED ALL TESTS************')
    IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99043)
    99043 FORMAT (/' ***********FFT ROUTINES FAILED SOME TESTS***********')
    RETURN
  END SUBROUTINE FFTQX
END MODULE TEST51_MOD
!** TEST51
PROGRAM TEST51
  USE TEST51_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  J1
  !***
  ! **Type:**      SINGLE PRECISION (TEST51-S)
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
  !        COSQB    COSQF    COSQI    COST     COSTI    EZFFTB
  !        EZFFTF   RFFTB    RFFTF    RFFTI    SINQB    SINQF
  !        SINQI    SINT     SINTI
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  FFTQX, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST51
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
  !     Test FFT package
  !
  CALL FFTQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST51 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST51 *************')
  ENDIF
  STOP
END PROGRAM TEST51
