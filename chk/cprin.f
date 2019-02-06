*DECK CPRIN
      SUBROUTINE CPRIN (LUN, NUM1, KPRINT, IP, EXACT, RESULT, ABSERR,
     +   NEVAL, IERV, LIERV)
C***BEGIN PROLOGUE  CPRIN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CQAG, CQAG, CQAGI, CQAGP, CQAGS, CQAWC,
C            CQAWF, CQAWO, CQAWS, and CQNG.
C***LIBRARY   SLATEC
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C   This program is called by the (single precision) Quadpack quick
C   check routines for printing out their messages.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   910627  Code completely rewritten.  (WRB)
C***END PROLOGUE  CPRIN
C     .. Scalar Arguments ..
      REAL ABSERR, EXACT, RESULT
      INTEGER IP, KPRINT, LIERV, LUN, NEVAL, NUM1
C     .. Array Arguments ..
      INTEGER IERV(*)
C     .. Local Scalars ..
      REAL ERROR
      INTEGER IER, K
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C***FIRST EXECUTABLE STATEMENT  CPRIN
      IER = IERV(1)
      ERROR = ABS(EXACT-RESULT)
C
      IF (KPRINT .GE. 2) THEN
        IF (IP.EQ.1) THEN
          IF (KPRINT .GE. 3) THEN
C
C           Write PASS message.
C
            WRITE (UNIT=LUN, FMT=9000) NUM1
          ENDIF
        ELSE
C
C         Write failure messages.
C
          WRITE (UNIT=LUN, FMT=9010) NUM1
          IF (NUM1 .EQ. 0) WRITE (UNIT=LUN, FMT=9020)
          IF (NUM1 .GT. 0) WRITE (UNIT=LUN, FMT=9030) NUM1
          IF (LIERV .GT. 1) WRITE (UNIT=LUN, FMT=9040) (IERV(K),
     +                      K=2,LIERV)
          IF (NUM1 .EQ. 6) WRITE (UNIT=LUN, FMT=9050)
          WRITE (UNIT=LUN, FMT=9060)
          WRITE (UNIT=LUN, FMT=9070)
          IF (NUM1 .NE. 5) THEN
            WRITE (UNIT=LUN, FMT=9080) EXACT,RESULT,ERROR,ABSERR,IER,
     +                                 NEVAL
          ELSE
            WRITE (LUN,FMT=9090) RESULT,ABSERR,IER,NEVAL
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
C
 9000 FORMAT (' TEST ON IER = ', I2, ' PASSED')
 9010 FORMAT (' TEST ON IER = ', I1, ' FAILED.')
 9020 FORMAT (' WE MUST HAVE IER = 0, ERROR.LE.ABSERR AND ABSERR.LE',
     +        '.MAX(EPSABS,EPSREL*ABS(EXACT))')
 9030 FORMAT (' WE MUST HAVE IER = ', I1)
 9040 FORMAT (' OR IER =     ', 8(I1,2X))
 9050 FORMAT (' RESULT, ABSERR, NEVAL AND EVENTUALLY LAST SHOULD BE',
     +        ' ZERO')
 9060 FORMAT (' WE HAVE   ')
 9070 FORMAT (7X, 'EXACT', 11X, 'RESULT', 6X, 'ERROR', 4X, 'ABSERR',
     +        4X, 'IER     NEVAL', /, ' ', 42X,
     +        '(EST.ERR.)(FLAG)(NO F-EVAL)')
 9080 FORMAT (' ', 2(E15.7,1X), 2(E9.2,1X), I4, 4X, I6)
 9090 FORMAT (5X, 'INFINITY', 4X, E15.7, 11X, E9.2, I5, 4X, I6)
      END
