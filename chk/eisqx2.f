*DECK EISQX2
      SUBROUTINE EISQX2 (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  EISQX2
C***PURPOSE  Quick check for SSIEV, CHIEV and SSPEV.
C***LIBRARY   SLATEC
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Kahaner, D. K., (NBS)
C***DESCRIPTION
C
C     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
C     SSIEV, CHIEV AND SSPEV.  THE EIGENVALUES OF INPUT MATRIX
C     A(.,.) ARE STORED IN EK(.).  RELERR IS THE RELATIVE
C     ACCURACY REQUIRED FOR THEM TO PASS.
C
C***ROUTINES CALLED  CHIEV, R1MACH, SSIEV, SSPEV
C***REVISION HISTORY  (YYMMDD)
C   800808  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900405  CALL to XERROR replaced by message to LUN.  (WRB)
C***END PROLOGUE  EISQX2
      INTEGER KPRINT,IPASS,LUN
      INTEGER LDA,N,LDV,JOB,I,J,ID
      REAL A1(4,4),A2(10),AP(10),E(4),V(4,4),EK(4),W(16)
      REAL ERR,ERRI,RELERR
      COMPLEX AC(4,4),VC(4,4)
      EQUIVALENCE (V,VC)
      DATA LDA,N,LDV / 3*4 /
      DATA AP / 5., 4., 5., 1., 1., 4., 1., 1., 2., 4. /
      DATA EK / 1., 2., 5., 10. /
C***FIRST EXECUTABLE STATEMENT  EISQX2
      IPASS = 1
      RELERR = SQRT(R1MACH(4))
      ID = 0
      DO 20 J=1,N
         DO 10 I=1,J
            ID = ID + 1
            A1(I,J) = AP(ID)
            A2(ID) = AP(ID)
            AC(I,J) = CMPLX(AP(ID),0.)
   10       CONTINUE
   20    CONTINUE
      JOB = 1
      CALL CHIEV(AC,LDA,N,E,VC,LDV,W,JOB,INFO)
      IF (INFO .NE. 0) THEN
         IF (KPRINT .GE. 2) WRITE (LUN, 688) 'CHIEV', INFO
         IPASS = 0
      ENDIF
      DO 40 J=1,N
         ERR = ABS(E(J) - EK(1))
         ID = 1
         DO 30 I=2,N
            ERRI = ABS(E(J) - EK(I))
            IF (ERRI .LT. ERR) ID = I
            ERR = MIN(ERRI,ERR)
   30       CONTINUE
         IF (ABS(E(J)-EK(ID))/ABS(EK(ID)) .GE. RELERR) IPASS = 0
   40    CONTINUE
      CALL SSIEV(A1,LDA,N,E,W,JOB,INFO)
      IF (INFO .NE. 0) THEN
         IF (KPRINT .GE. 2) WRITE (LUN, 688) 'SSIEV', INFO
         IPASS = 0
      ENDIF
      DO 60 J=1,N
         ERR = ABS(E(J) - EK(1))
         ID = 1
         DO 50 I=2,N
            ERRI = ABS(E(J) - EK(I))
            IF (ERRI .LT. ERR) ID = I
            ERR = MIN(ERRI,ERR)
   50       CONTINUE
         IF (ABS(E(J)-EK(ID))/ABS(EK(ID)) .GE. RELERR) IPASS = 0
   60    CONTINUE
      JOB = 0
      CALL SSPEV(A2,N,E,V,LDV,W,JOB,INFO)
      IF (INFO .NE. 0) THEN
         IF (KPRINT .GE. 2) WRITE (LUN, 688) 'SSPEV', INFO
         IPASS = 0
      ENDIF
      DO 80 J=1,N
         ERR = ABS(E(J) - EK(1))
         ID = 1
         DO 70 I=2,N
            ERRI = ABS(E(J) - EK(I))
            IF (ERRI .LT. ERR) ID = I
            ERR = MIN(ERRI,ERR)
   70       CONTINUE
         IF (ABS(E(J)-EK(ID))/ABS(EK(ID)) .GE. RELERR) IPASS = 0
   80    CONTINUE
      IF (KPRINT.GE.2 .AND. IPASS.NE.0) WRITE (LUN,684)
  684 FORMAT(25H EISQX2 PASSES ALL TESTS.)
      IF (KPRINT.GE.1 .AND. IPASS.EQ.0) WRITE (LUN,686)
  686 FORMAT(25H EISQX2 FAILS SOME TESTS.)
  688 FORMAT (1X, 'Eigenvalue iteration failed to converge in ', A5,
     +        ', INFO = ', I4)
      RETURN
      END
