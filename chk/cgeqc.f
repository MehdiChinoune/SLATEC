*DECK CGEQC
      SUBROUTINE CGEQC (LUN, KPRINT, NERR)
C***BEGIN PROLOGUE  CGEQC
C***PURPOSE  Quick check for CGEFS and CGEIR.
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (SGEQC-S, DGEQC-D, CGEQC-C)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Jacobsen, Nancy, (LANL)
C***DESCRIPTION
C
C   Let A*X=B be a COMPLEX linear system where the
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
C***ROUTINES CALLED  CGEFS, CGEIR
C***REVISION HISTORY  (YYMMDD)
C   801029  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920601  Code reworked and TYPE section added.  (RWC, WRB)
C***END PROLOGUE  CGEQC
C     .. Scalar Arguments ..
      INTEGER KPRINT, LUN, NERR
C     .. Local Scalars ..
      COMPLEX XA, XB
      INTEGER I, IND, INDX, ITASK, J, KPROG, LDA, N
C     .. Local Arrays ..
      COMPLEX A(3,3), ATEMP(5,3), B(3), BTEMP(3), BXEX(3), WORK(12)
      INTEGER IWORK(3)
      CHARACTER LIST(2)*4
C     .. External Subroutines ..
      EXTERNAL CGEFS, CGEIR
C     .. Intrinsic Functions ..
      INTRINSIC ABS, AIMAG, REAL
C     .. Statement Functions ..
      REAL DELX
C     .. Data statements ..
      DATA A /(2., 3.), (1., 1.),  (1., 2.),
     +        (2., 0.), (1., -1.), (0., 0.),
     +        (0., 0.), (2., 5.),  (3., 2.)/
      DATA B /(-1., 1.), (-5., 4.), (-4., 7.)/
      DATA BXEX /(.21459E-01, .209012E+01), (.261373E+01, -.162231E+01),
     +           (.785407E+00, .109871E+01)/
      DATA LIST /'GEFS', 'GEIR'/
C     .. Statement Function definitions ..
      DELX(XA,XB) = ABS(REAL(XA-XB)) + ABS(AIMAG(XA-XB))
C***FIRST EXECUTABLE STATEMENT  CGEQC
      N = 3
      LDA = 5
      NERR = 0
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
      DO 180 KPROG=1,2
C
C     First test case - normal
C
        ITASK = 1
        DO 100 I=1,N
          BTEMP(I) = B(I)
  100   CONTINUE
        DO 120 J=1,N
          DO 110 I=1,N
            ATEMP(I,J) = A(I,J)
  110     CONTINUE
  120   CONTINUE
        IF (KPROG .EQ. 1) THEN
          CALL CGEFS (ATEMP, LDA, N, BTEMP, ITASK, IND, WORK, IWORK)
        ELSE
          CALL CGEIR (ATEMP, LDA, N, BTEMP, ITASK, IND, WORK, IWORK)
        ENDIF
        IF (IND .LT. 0) THEN
          IF (KPRINT .GE. 2) WRITE (LUN, FMT=9020) LIST(KPROG), IND
          NERR = NERR + 1
        ENDIF
C
C       Calculate error for first test
C
        INDX = 0
        DO 130 I=1,N
          IF (DELX(BXEX(I),BTEMP(I)) .GT. .0001) INDX = INDX + 1
  130   CONTINUE
        IF (INDX .EQ. 0) THEN
          IF(KPRINT .GE. 3) WRITE (LUN, FMT=9010) LIST(KPROG)
        ELSE
          IF(KPRINT .GE. 2) WRITE (LUN, FMT=9020) LIST(KPROG)
          NERR = NERR + 1
        ENDIF
C
C       Second test case - singular matrix
C
        ITASK = 1
        DO 140 I=1,N
          BTEMP(I) = B(I)
  140   CONTINUE
        DO 160 J=1,N
          DO 150 I=1,N
            ATEMP(I,J) = A(I,J)
  150     CONTINUE
  160   CONTINUE
        DO 170 J=1,N
          ATEMP(1,J) = (0.E0,0.E0)
  170   CONTINUE
        IF (KPROG .EQ. 1) THEN
          CALL CGEFS (ATEMP, LDA, N, BTEMP, ITASK, IND,  WORK, IWORK)
        ELSE
          CALL CGEIR (ATEMP, LDA, N, BTEMP, ITASK, IND, WORK, IWORK)
        ENDIF
        IF (IND .EQ. -4) THEN
          IF (KPRINT .GE. 3) WRITE (LUN, FMT=9030) LIST(KPROG)
        ELSE
          IF (KPRINT .GE. 2) WRITE (LUN, FMT=9040) LIST(KPROG), IND
          NERR = NERR + 1
        ENDIF
  180 CONTINUE
C
      IF (KPRINT.GE.3 .AND. NERR.EQ.0) WRITE (LUN,9050)
      IF (KPRINT.GE.2 .AND. NERR.NE.0) WRITE (LUN,9060)
      RETURN
C
 9000 FORMAT (//, 2X, 'CGEFS and CGEIR Quick Check' /)
 9010 FORMAT (/, 5X, 'C', A, ' Normal test PASSED')
 9020 FORMAT (/, 5X, 'C', A, ' Test FAILED')
 9030 FORMAT (/, 5X, 'C', A, ' Singular test PASSED')
 9040 FORMAT (/, 5X, 'C', A, ' Singular test FAILED, IND=', I3)
 9050 FORMAT (/, 2X, 'CGEFS and CGEIR Quick Check PASSED' /)
 9060 FORMAT (/, 2X, 'CGEFS and CGEIR Quick Check FAILED' /)
      END
