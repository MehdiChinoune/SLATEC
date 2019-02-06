*DECK DGEQC
      SUBROUTINE DGEQC (LUN, KPRINT, NERR)
C***BEGIN PROLOGUE  DGEQC
C***PURPOSE  Quick check for DGEFS.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SGEQC-S, DGEQC-D, CGEQC-C)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Jacobsen, Nancy, (LANL)
C***DESCRIPTION
C
C   Let A*X=B be a DOUBLE PRECISION linear system where the
C   matrix is of the proper type for the Linpack subroutines
C   being called.  The values of A and B and the pre-computed
C   values of BXEX (the solution vector) are given in DATA
C   statements.  The computed test results for X are compared to
C   the stored pre-computed values.  Failure of the test occurs
C   when there is less than 80% agreement between the absolute
C   values.  There are 2 tests - one for the normal case and one
C   for the singular case.  A message is printed indicating
C   whether each subroutine has passed or failed for each case.
C
C   On return, NERR (INTEGER type) contains the total count of
C   all failures detected.
C
C***ROUTINES CALLED  D1MACH, DGEFS
C***REVISION HISTORY  (YYMMDD)
C   801022  DATE WRITTEN
C   891009  Removed unreferenced statement label.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920601  Code reworked and TYPE section added.  (RWC, WRB)
C***END PROLOGUE  DGEQC
C     .. Scalar Arguments ..
      INTEGER KPRINT, LUN, NERR
C     .. Local Scalars ..
      DOUBLE PRECISION ERRCMP, ERRMAX
      INTEGER I, IND, ITASK, J, KPROG, LDA, N
C     .. Local Arrays ..
      DOUBLE PRECISION A(5,4), ATEMP(5,4), B(4), BTEMP(4), BXEX(4),
     +                 WORK(20)
      INTEGER IWORK(4)
      CHARACTER LIST(2)*4
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     .. External Subroutines ..
      EXTERNAL DGEFS
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX
C     .. Data statements ..
      DATA A /5.0D0,  1.0D0,  0.3D0, 2.1D0, 0.0D0,
     +       -1.0D0, -0.5D0,  1.0D0, 1.0D0, 0.0D0,
     +        4.5D0, -1.0D0, -1.7D0, 2.0D0, 0.0D0,
     +        0.5D0,  2.0D0,  0.6D0, 1.3D0, 0.0D0/
      DATA B /0.0D0, 3.5D0, 3.6D0, 2.4D0/
      DATA BXEX /0.10D+01, 0.10D+01, -0.10D+01, 0.10D+01/
      DATA LIST /'GEFS', 'GEIR'/
C***FIRST EXECUTABLE STATEMENT  DGEQC
      N = 4
      LDA = 5
      NERR = 0
      ERRCMP = D1MACH(4)**0.8D0
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
      KPROG = 1
C
C     First test case - normal
C
      ITASK = 1
      DO 100 I=1,N
        BTEMP(I) = B(I)
  100 CONTINUE
      DO 120 J=1,N
        DO 110 I=1,N
          ATEMP(I,J) = A(I,J)
  110   CONTINUE
  120 CONTINUE
      CALL DGEFS (ATEMP, LDA, N, BTEMP, ITASK, IND, WORK, IWORK)
      IF (IND .LT. 0) THEN
        IF (KPRINT .GE. 2) WRITE (LUN, FMT=9020) LIST(KPROG), IND
        NERR = NERR + 1
      ENDIF
C
C     Calculate error for first test
C
      ERRMAX = 0.0D0
C
      DO 130 I=1,N
        ERRMAX = MAX(ERRMAX,ABS(BTEMP(I)-BXEX(I)))
  130 CONTINUE
      IF (ERRCMP .GT. ERRMAX) THEN
        IF (KPRINT .GE. 3) WRITE (LUN, FMT=9010) LIST(KPROG)
      ELSE
        IF (KPRINT .GE. 2) WRITE (LUN, FMT=9020) LIST(KPROG), ERRMAX
        NERR = NERR + 1
      ENDIF
C
C     Second test case - singular matrix
C
      ITASK = 1
      DO 140 I=1,N
        BTEMP(I) = B(I)
  140 CONTINUE
      DO 160 J=1,N
        DO 150 I=1,N
          ATEMP(I,J) = A(I,J)
  150   CONTINUE
  160 CONTINUE
      DO 170 J=1,N
        ATEMP(1,J) = 0.0D0
  170 CONTINUE
      CALL DGEFS (ATEMP, LDA, N, BTEMP, ITASK, IND,  WORK, IWORK)
      IF (IND .EQ. -4) THEN
        IF (KPRINT .GE. 3) WRITE (LUN, FMT=9030) LIST(KPROG)
      ELSE
        IF (KPRINT .GE. 2) WRITE (LUN, FMT=9040) LIST(KPROG), IND
        NERR = NERR + 1
      ENDIF
C
      IF (KPRINT.GE.3 .AND. NERR.EQ.0) WRITE (LUN,9050)
      IF (KPRINT.GE.2 .AND. NERR.NE.0) WRITE (LUN,9060)
      RETURN
C
 9000 FORMAT (//, 2X, 'DGEFS Quick Check' /)
 9010 FORMAT (/, 5X, 'D', A, ' Normal test PASSED')
 9020 FORMAT (/, 5X, 'D', A, ' Test FAILED, MAX ABS(ERROR) is', E13.5)
 9030 FORMAT (/, 5X, 'D', A, ' Singular test PASSED')
 9040 FORMAT (/, 5X, 'D', A, ' Singular test FAILED, IND=', I3)
 9050 FORMAT (/, 2X, 'DGEFS Quick Check PASSED' /)
 9060 FORMAT (/, 2X, 'SGEFS and SGEIR Quick Check FAILED' /)
      END
