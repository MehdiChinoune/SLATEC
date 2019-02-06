*DECK EISQX1
      SUBROUTINE EISQX1 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  EISQX1
C***PURPOSE  Quick check for SGEEV and CGEEV.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
C     SGEEV AND CGEEV.  THE EIGENVALUES OF INPUT MATRIX A(.,.)
C     ARE STORED IN EK(.).  RELERR IS THE RELATIVE ACCURACY
C     REQUIRED FOR THEM TO PASS.
C
C***ROUTINES CALLED  CGEEV, R1MACH, SGEEV
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900405  CALL to XERROR replaced by message to LUN.  (WRB)
C***END PROLOGUE  EISQX1
      INTEGER KPRINT,IPASS,LUN
      INTEGER LDA,N,LDV,JOB,I,J,ID
      REAL A(3,3),EK(3),W(9)
      REAL ERR,ERRI,RELERR,RECJ
      COMPLEX AC(3,3),EC(3),VC(3,3)
      DATA LDA,N,LDV / 3*3 /
      DATA A / 1., -2., 6., -1., 0., -3., 2., 5., 6. /
      DATA EK / -1., 3., 5. /
C***FIRST EXECUTABLE STATEMENT  EISQX1
      IPASS = 1
      RELERR = SQRT(R1MACH(4))
      DO 20 J=1,N
         DO 10 I=1,N
            AC(I,J) = CMPLX(A(I,J),0.)
   10       CONTINUE
   20    CONTINUE
      JOB = 1
      CALL CGEEV(AC,LDA,N,EC,VC,LDV,W,JOB,INFO)
      IF (INFO .NE. 0) THEN
         IF (KPRINT .GE. 2) WRITE (LUN, 688) 'CGEEV', INFO
         IPASS = 0
      ENDIF
      DO 40 J=1,N
         ERR = ABS(AIMAG(EC(J)))
         IF (ERR .GE. RELERR) IPASS = 0
         RECJ = REAL(EC(J))
         ERR = ABS(RECJ - EK(1))
         ID = 1
         DO 30 I=2,N
            ERRI = ABS(RECJ - EK(I))
            IF (ERRI .LT. ERR) ID = I
            ERR = MIN(ERRI,ERR)
   30       CONTINUE
         IF (ABS(RECJ-EK(ID))/ABS(EK(ID)) .GE. RELERR) IPASS = 0
   40    CONTINUE
      JOB = 0
      CALL SGEEV(A,LDA,N,EC,VC,LDV,W,JOB,INFO)
      IF (INFO .NE. 0) THEN
         IF (KPRINT .GE. 2) WRITE (LUN, 688) 'SGEEV', INFO
         IPASS = 0
      ENDIF
      DO 60 J=1,N
         ERR = ABS(AIMAG(EC(J)))
         IF (ERR .GE. RELERR) IPASS = 0
         RECJ = REAL(EC(J))
         ERR = ABS(RECJ - EK(1))
         ID = 1
         DO 50 I=2,N
            ERRI = ABS(RECJ - EK(I))
            IF (ERRI .LT. ERR) ID = I
            ERR = MIN(ERRI,ERR)
   50       CONTINUE
         IF (ABS(RECJ-EK(ID))/ABS(EK(ID)) .GE. RELERR) IPASS = 0
   60    CONTINUE
      IF (KPRINT.GE.2 .AND. IPASS.NE.0) WRITE (LUN,670)
  670 FORMAT(25H EISQX1 PASSES ALL TESTS.)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (LUN,680)
  680 FORMAT(25H EISQX1 FAILS SOME TESTS.)
  688 FORMAT (1X, 'Eigenvalue iteration failed to converge in ', A5,
     +        ', INFO = ', I4)
      RETURN
      END
