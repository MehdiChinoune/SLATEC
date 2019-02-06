*DECK FCNQX1
      SUBROUTINE FCNQX1 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  FCNQX1
C***SUBSIDIARY
C***PURPOSE  THIS IS A QUICK CHECK PROGRAM FOR THE SUBROUTINES XLEGF
C            AND XNRMP WHICH CALCULATE LEGENDRE FUNCTIONS
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      SINGLE PRECISION (FCNQX1-S, FCNQX2-D)
C***KEYWORDS  LEGENDRE FUNCTIONS, QUICK CHECK
C***AUTHOR  LOZIER, DANIEL W., (NIST)
C           SMITH, JOHN M., (NIST AND GEORGE MASON UNIVERSITY)
C***REFERENCES  OLVER AND SMITH,J.COMPUT.PHYSICS,51(1983),NO.3,502-518.
C               SMITH, OLVER AND LOZIER,ACM TRANS MATH SOFTW,7(1981),
C                 NO.1,93-105.
C***ROUTINES CALLED  XCON, XCSRT, XERCLR, XLEGF, XNRMP, XSET, XSETF
C***REVISION HISTORY  (YYMMDD)
C   881020  DATE WRITTEN
C   900306  Added SLATEC prologue to this routine. (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C   910104  Changed to print variable number of decimals. (DWL and JMS)
C***END PROLOGUE  FCNQX1
C
      CHARACTER*34 FMT, FMTF, FMTI
      INTEGER LUN,KPRINT,IPASS
      DIMENSION P(10),Q(10),R(10),C1(10),C2(10),IP(10),IQ(10),IR(10)
      DIMENSION IC1(10),IC2(10),PN(10),IPN(10)
      REAL P,Q,R,C1,C2,PN
      REAL DEG,THETA,DNU1,DZERO
      REAL X11,X12,X13,X21,X22,X23
      REAL NU
C
C***FIRST EXECUTABLE STATEMENT  FCNQX1
C
      IF(KPRINT.GE.2) WRITE(LUN,1)
    1 FORMAT(' ** TEST SINGLE PRECISION LEGENDRE FUNCTION ROUTINES',
     2' IN FCNPAK ** ',/)
      IPASS=1
      IRAD=0
      NRADPL=0
      DZERO=0.0
      NBITS=0
      CALL XSET(IRAD,NRADPL,DZERO,NBITS,IERROR)
      IF(IERROR.NE.0) IPASS=0
      IERR=0
      DNU1=2000.4
      IF(I1MACH(13)*LOG10(REAL(I1MACH(10))).LT.150.) DNU1=100.4
      IF (KPRINT.LE.2) GO TO 150
      IF (I1MACH(13).LT.500) WRITE(LUN,24)
   24 FORMAT(' ON COMPUTERS WITH MAXIMUM EXPONENT LESS THAN 500, SMALL'/
     1' TEST VALUES FOR NU, MU ARE USED. IF LARGER THAN OR EQUAL 500,'/
     2' LARGER VALUES ARE USED. THIS COMPUTER USES THE SMALLER VALUES.')
      IF (I1MACH(13).GE.500) WRITE(LUN,26)
   26 FORMAT(' ON COMPUTERS WITH MAXIMUM EXPONENT LESS THAN 500, SMALL'/
     1' TEST VALUES FOR NU, MU ARE USED. IF LARGER THAN OR EQUAL 500,'/
     2' LARGER VALUES ARE USED. THIS COMPUTER USES THE LARGER VALUES.')
  150 CONTINUE
      NUDIFF=5
      MU1=DNU1
      MU2=MU1
      DEG=0.1
      THETA=DEG*4.*ATAN(1.0)/180.0
C
C In TEST 1 the Legendre functions P (of both positive and negative
C order) and Q are calculated.  Large values of mu and nu are used
C so that it is necessary to use extended range arithmetic.  The
C values of the Casoratians should be approximately equal to 1.0.
C The check which is applied is to verify that the difference between
C the Casoratians and 1.0 is less that 10.**(6-NDEC), where NDEC =
C INT((D-1)*LOG10(R)), D = I1MACH(11) = significand length, R =
C I1MACH(10) = radix. The value of IERROR should always be returned
C as zero. This test uses the programs
C XLEGF, XPQNU, XPSI, XQNU, XPMUP, XSET, XADD,
C XADJ, XCSRT, XRED, XC210, and XCON.
C
      ISUM=0
      NDEC = (I1MACH(11)-1) * LOG10(REAL(I1MACH(10)))
C Formats that depend on NDEC ...
      FMT(1:20)='(1X, 6X, 4H   (,E30.'
      WRITE(FMT(21:22),'(I2)') NDEC
      FMT(23:34)=',1H,,I8,1H))'
      FMTF(1:20)='(1X,F6.1,4H   (,E30.'
      WRITE(FMTF(21:22),'(I2)') NDEC
      FMTF(23:34)=',1H,,I8,1H))'
      FMTI(1:20)='(1X, I6, 4H   (,E30.'
      WRITE(FMTI(21:22),'(I2)') NDEC
      FMTI(23:34)=',1H,,I8,1H))'
      IF (KPRINT.GT.2) WRITE(LUN,2) MU1, DEG
    2 FORMAT(/
     1' TEST 1, FIXED MU = ',I4,' AND THETA = ',F3.1,
     1' DEGREES, RECURRENCE IN NU,'/
     2'         CASORATIS SHOULD = 1.0')
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,1,P,IP,IERROR)
      ISUM=ISUM+IERROR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,2,Q,IQ,IERROR)
      ISUM=ISUM+IERROR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,3,R,IR,IERROR)
      ISUM=ISUM+IERROR
      CALL XCSRT(DNU1,NUDIFF,MU1,MU2,THETA,P,Q,R,IP,IQ,IR,
     1 C1,IC1,C2,IC2,IERROR)
      ISUM=ISUM+IERROR
      DO 20 I=1,6
      CALL XCON(P(I),IP(I),IERROR)
      ISUM=ISUM+IERROR
      CALL XCON(Q(I),IQ(I),IERROR)
      ISUM=ISUM+IERROR
      CALL XCON(R(I),IR(I),IERROR)
      ISUM=ISUM+IERROR
   20 CONTINUE
      X11=P(1)
      IX11=IP(1)
      X12=R(1)
      IX12=IR(1)
      X13=Q(1)
      IX13=IQ(1)
      IF(KPRINT.GT.2) THEN
        WRITE(LUN,'(A)') '     NU   CASORATI 1'
        NU=DNU1
        DO 25 I=1,5
        WRITE(LUN,FMTF) NU,C1(I),IC1(I)
        NU=NU+1.
   25   CONTINUE
        WRITE(LUN,'(A)') '     NU   CASORATI 2'
        NU=DNU1
        DO 30 I=1,5
        WRITE(LUN,FMTF) NU,C2(I),IC2(I)
        NU=NU+1.
   30   CONTINUE
      ENDIF
      DO 35 I=1,5
      IF(ABS(1.0-C1(I)).GE.10.0E0**(6-NDEC)) GO TO 40
      IF(ABS(1.0-C2(I)).GE.10.0E0**(6-NDEC)) GO TO 40
   35 CONTINUE
      IF(ISUM.NE.0) GO TO 40
      IF (KPRINT.GE.2) WRITE(LUN,8)
    8 FORMAT(' ***** TEST 1 (SINGLE PRECISION) PASSED *****')
      GO TO 50
   40 IF(KPRINT.GE.1) WRITE(LUN,7)
    7 FORMAT(' ***** TEST 1 (SINGLE PRECISION) FAILED *****')
      IERR=IERR+1
      IPASS=0
   50 NUDIFF=0
      MU1=MU2-5
C
C In TEST 2 P (of positive and negative order) and Q are again
C calculated but in this test the recurrence is in the mu-wise direction
C rather than in the nu-wise direction as was the case before.  The same
C programs are used except that XQNU is not used and XQMU and XPMU
C are used. Again the criterion for passing the test is that the
C Casoratians differ from 1.0 by less than 10.0**(6-NDEC). The value
C of IERROR should always be returned as zero.
C
      ISUM=0
      IF (KPRINT.GT.2) WRITE(LUN,9) DNU1, DEG
    9 FORMAT(/
     1' TEST 2, FIXED NU = ',F6.1,' AND THETA = ',F3.1,
     1' DEGREES, RECURRENCE IN MU,'/
     2'         CASORATIS SHOULD = 1.0')
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,1,P,IP,IERROR)
      ISUM=ISUM+IERROR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,2,Q,IQ,IERROR)
      ISUM=ISUM+IERROR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,3,R,IR,IERROR)
      ISUM=ISUM+IERROR
      CALL XCSRT(DNU1,NUDIFF,MU1,MU2,THETA,P,Q,R,IP,IQ,IR,
     1 C1,IC1,C2,IC2,IERROR)
      ISUM=ISUM+IERROR
      DO 60 I=1,6
      CALL XCON(P(I),IP(I),IERROR)
      ISUM=ISUM+IERROR
      CALL XCON(Q(I),IQ(I),IERROR)
      ISUM=ISUM+IERROR
      CALL XCON(R(I),IR(I),IERROR)
      ISUM=ISUM+IERROR
   60 CONTINUE
      X21=P(6)
      IX21=IP(6)
      X22=R(6)
      IX22=IR(6)
      X23=Q(6)
      IX23=IQ(6)
      IF(KPRINT.GT.2) THEN
        WRITE(LUN,'(A)') '     MU   CASORATI 3'
        MU=MU1
        DO 65 I=1,5
        WRITE(LUN,FMTI) MU,C1(I),IC1(I)
        MU=MU+1
   65   CONTINUE
        WRITE(LUN,'(A)') '     MU   CASORATI 4'
        MU=MU1
        DO 70 I=1,5
        WRITE(LUN,FMTI) MU,C2(I),IC2(I)
        MU=MU+1
   70   CONTINUE
      ENDIF
      DO 75 I=1,5
      IF(ABS(1.0-C1(I)).GE.10.0E0**(6-NDEC)) GO TO 80
      IF(ABS(1.0-C2(I)).GE.10.0E0**(6-NDEC)) GO TO 80
      IF(ISUM.NE.0) GO TO 80
   75 CONTINUE
      IF(KPRINT.GE.2) WRITE(LUN,12)
   12 FORMAT(' ***** TEST 2 (SINGLE PRECISION) PASSED *****')
      GO TO 85
   80 IF(KPRINT.GE.1) WRITE(LUN,11)
   11 FORMAT(' ***** TEST 2 (SINGLE PRECISION) FAILED *****')
      IERR=IERR+1
      IPASS=0
C
C In TEST 3 values of P and Q which were calculated in two different
C manners, one by nu-wise recurrence in TEST 1 and one by mu-wise
C recurrence in TEST 2, are compared.  Again, the criterion for success
C is a relative error of less than 10.0**(6-NDEC).
C
   85 IF(KPRINT.GT.2) THEN
      WRITE(LUN,13) DEG, MU2, DNU1
   13 FORMAT(/
     1' TEST 3, COMPARISON OF VALUES FROM TEST 1 AND TEST 2',
     1' WITH THETA = ',F3.1,' DEGREES,'/
     2'         MU = ',I4,' AND NU = ',F6.1)
        WRITE(LUN,'(A)') '          P(-MU,NU)'
        WRITE(LUN,FMT)   X11,IX11
        WRITE(LUN,FMT)   X21,IX21
        WRITE(LUN,'(A)') '          P(MU,NU)'
        WRITE(LUN,FMT)   X12,IX12
        WRITE(LUN,FMT)   X22,IX22
        WRITE(LUN,'(A)') '          Q(MU,NU)'
        WRITE(LUN,FMT)   X13,IX13
        WRITE(LUN,FMT)   X23,IX23
      ENDIF
      IF(ABS((X11-X21)/X11).GE.10.0E0**(6-NDEC)) GO TO 90
      IF(ABS((X12-X22)/X12).GE.10.0E0**(6-NDEC)) GO TO 90
      IF(ABS((X13-X13)/X13).GE.10.0E0**(6-NDEC)) GO TO 90
      IF(IX11.NE.IX21) GO TO 90
      IF(IX12.NE.IX22) GO TO 90
      IF(IX13.NE.IX23) GO TO 90
      IF(KPRINT.GE.2) WRITE(LUN,15)
   15 FORMAT(' ***** TEST 3 (SINGLE PRECISION) PASSED *****')
      GO TO 100
   90 IF(KPRINT.GE.1) WRITE(LUN,16)
   16 FORMAT(' ***** TEST 3 (SINGLE PRECISION) FAILED *****')
      IERR=IERR+1
      IPASS=0
  100 CONTINUE
C
C In TEST 4, the value of the normalized Legendre function as
C calculated by XLEGF and XPNRM is compared to the same value
C as calculated by the program XNRMP.  Again the criterion is a
C relative error of less than 10.0**(6-NDEC). The value of IERROR
C should always be returned as zero.
C
      ISUM=0
      DNU1=100.0
      NUDIFF=0
      MU1=10
      MU2=10
      IF(KPRINT.GT.2) WRITE(LUN,17) DEG, MU1, DNU1
   17 FORMAT(/
     1' TEST 4, COMPARISON OF VALUES FROM XLEGF AND XNRMP',
     1' WITH THETA = ',F3.1,' DEGREES,'/
     2'         MU = ',I4,' AND NU = ',F6.1)
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,4,PN,IPN,IERROR)
      ISUM=ISUM+IERROR
      X11=PN(1)
      IX11=IPN(1)
      NU1=100
      CALL XNRMP(NU1,MU1,MU2,THETA,2,PN,IPN,ISIG,IERROR)
      ISUM=ISUM+IERROR
      X21=PN(1)
      IX21=IPN(1)
      IF(KPRINT.GT.2) THEN
        WRITE(LUN,'(A)') '          NORMALIZED P'
        WRITE(LUN,FMT) X11,IX11
        WRITE(LUN,FMT) X21,IX21
      ENDIF
      IF(ABS((X11-X21)/X11).GE.10.0E0**(6-NDEC)) GO TO 110
      IF(IX11.NE.IX21) GO TO 110
      IF(ISUM.NE.0) GO TO 110
      IF(KPRINT.GE.2) WRITE(LUN,18)
   18 FORMAT(' ***** TEST 4 (SINGLE PRECISION) PASSED *****')
      GO TO 120
  110 IF(KPRINT.GE.1) WRITE(LUN,19)
   19 FORMAT(' ***** TEST 4 (SINGLE PRECISION) FAILED *****')
      IERR=IERR+1
      IPASS=0
  120 CONTINUE
C
C In TEST 5 errors are purposely made in input so as to test error
C handling capability. First, an incorrect value of ID is given. Then
C both NUDIFF and MU2-MU1 are non-zero. Finally, an incorrect value
C of THETA is given. In each case the value of the error indicator
C IERROR should equal the error number as returned by the error
C handling package XERROR (which includes XSETF, XERCLR, and NUMXER).
C
      CALL XSETF(-1)
      IF (KPRINT.LE.2) CALL XSETF(0)
      IF (KPRINT.GT.2) WRITE(LUN,23)
   23 FORMAT(/' TEST 5, TEST OF ERROR HANDLING. 3 ERROR MESSAGES',
     1' SHOULD BE PRINTED.')
      NUDIFF=0
      MU2=MU1
      ID=5
      CALL XERCLR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,ID,P,IP,IERROR)
      N=NUMXER(NERR)
      IF (N.NE.IERROR) GO TO 125
      MU2=MU1+5
      NUDIFF=5
      CALL XERCLR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,1,P,IP,IERROR)
      N=NUMXER(NERR)
      IF(N.NE.IERROR) GO TO 125
      NUDIFF=0
      THETA=2.0
      CALL XERCLR
      CALL XLEGF(DNU1,NUDIFF,MU1,MU2,THETA,1,P,IP,IERROR)
      N=NUMXER(NERR)
      IF(N.NE.IERROR) GO TO 125
      IF(KPRINT.GE.2) WRITE(LUN,28)
   28 FORMAT(' ***** TEST 5 (SINGLE PRECISION) PASSED *****')
      GO TO 135
  125 IF(KPRINT.GE.1) WRITE(LUN,29)
   29 FORMAT(' ***** TEST 5 (SINGLE PRECISION) FAILED *****')
      IERR=IERR+1
      IPASS=0
  135 CONTINUE
      IF(IERR.EQ.0) GO TO 140
      IF(KPRINT.GE.2) WRITE(LUN,21) IERR
   21 FORMAT(/'  TESTS COMPLETED, NUMBER OF TESTS FAILED = ',I2)
  140 CONTINUE
      RETURN
      END
