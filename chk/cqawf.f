*DECK CQAWF
      SUBROUTINE CQAWF (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  CQAWF
C***PURPOSE  Quick check for QAWF.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CQAWF-S, CDQAWF-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  CPRIN, F0F, F1F, QAWF, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Added PASS/FAIL message and changed the name of the first
C           argument.  (RWC)
C   910501  Added PURPOSE and TYPE records.  (WRB)
C***END PROLOGUE  CQAWF
C
C FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
C
      REAL A,ABSERR,R1MACH,EPSABS,EPMACH,
     *  ERROR,EXACT0,F0F,F1F,OMEGA,PI,RESULT,UFLOW,WORK
      INTEGER IER,IP,IPASS,KPRINT,LENW,LIMIT,LIMLST,LST,NEVAL
      DIMENSION IERV(3),IWORK(450),WORK(1425)
      EXTERNAL F0F,F1F
      DATA EXACT0/0.1422552162575912E+01/
      DATA PI/0.31415926535897932E+01/
C***FIRST EXECUTABLE STATEMENT  CQAWF
      IF (KPRINT.GE.2) WRITE (LUN, '(''1QAWF QUICK CHECK''/)')
C
C TEST ON IER = 0
C
      IPASS = 1
      MAXP1 = 21
      LIMLST = 50
      LIMIT = 200
      LENIW = LIMIT*2+LIMLST
      LENW = LENIW*2+MAXP1*25
      EPMACH = R1MACH(4)
      EPSABS = MAX(SQRT(EPMACH),0.1E-02)
      A = 0.0E+00
      OMEGA = 0.8E+01
      INTEGR = 2
      CALL QAWF(F0F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
     *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
      IERV(1) = IER
      IP = 0
      ERROR = ABS(EXACT0-RESULT)
      IF(IER.EQ.0.AND.ERROR.LE.ABSERR.AND.ABSERR.LE.EPSABS)
     *  IP = 1
      IF(IP.EQ.0) IPASS = 0
      CALL CPRIN(LUN,0,KPRINT,IP,EXACT0,RESULT,ABSERR,NEVAL,IERV,1)
C
C TEST ON IER = 1
C
      LIMLST = 3
      LENIW = 403
      LENW = LENIW*2+MAXP1*25
      CALL QAWF(F0F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
     *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
      IERV(1) = IER
      IP = 0
      IF(IER.EQ.1) IP = 1
      IF(IP.EQ.0) IPASS = 0
      CALL CPRIN(LUN,1,KPRINT,IP,EXACT0,RESULT,ABSERR,NEVAL,IERV,1)
C
C TEST ON IER = 3 OR 4 OR 1
C
      LIMLST = 50
      LENIW = LIMIT*2+LIMLST
      LENW = LENIW*2+MAXP1*25
      UFLOW = R1MACH(1)
      CALL QAWF(F1F,A,0.0E+00,1,UFLOW,RESULT,ABSERR,NEVAL,
     *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
      IERV(1) = IER
      IERV(2) = 4
      IERV(3) = 1
      IP = 0
      IF(IER.EQ.3.OR.IER.EQ.4.OR.IER.EQ.1) IP = 1
      IF(IP.EQ.0) IPASS = 0
      CALL CPRIN(LUN,3,KPRINT,IP,PI,RESULT,ABSERR,NEVAL,IERV,3)
C
C TEST ON IER = 6
C
      LIMLST = 50
      LENIW = 20
      CALL QAWF(F0F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
     *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
      IERV(1) = IER
      IP = 0
      IF(IER.EQ.6) IP = 1
      IF(IP.EQ.0) IPASS = 0
      CALL CPRIN(LUN,6,KPRINT,IP,EXACT0,RESULT,ABSERR,NEVAL,IERV,1)
C
C TEST ON IER = 7
C
      LIMLST = 50
      LENIW = 52
      LENW = LENIW*2+MAXP1*25
      CALL QAWF(F0F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
     *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
      IERV(1) = IER
      IP = 0
      IF(IER.EQ.7) IP = 1
      IF(IP.EQ.0) IPASS = 0
      CALL CPRIN(LUN,7,KPRINT,IP,EXACT0,RESULT,ABSERR,NEVAL,IERV,1)
C
      IF (KPRINT.GE.1) THEN
         IF (IPASS.EQ.0) THEN
            WRITE(LUN, '(/'' SOME TEST(S) IN CQAWF FAILED''/)')
         ELSEIF (KPRINT.GE.2) THEN
            WRITE(LUN, '(/'' ALL TEST(S) IN CQAWF PASSED''/)')
         ENDIF
      ENDIF
      RETURN
      END
