*DECK CHKXER
      SUBROUTINE CHKXER (SRNAMT, INFOT, NOUT, FATAL, KPRINT)
C***BEGIN PROLOGUE  CHKXER
C***SUBSIDIARY
C***PURPOSE  Test whether an error has been detected.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Tests whether an error has been detected.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  NUMXER
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CHKXER
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      INTEGER            INFOT, KPRINT, NOUT
      CHARACTER*6        SRNAMT
C     .. Local Scalars ..
      INTEGER            NERR
C     .. External Functions ..
      INTEGER            NUMXER
      EXTERNAL           NUMXER
C***FIRST EXECUTABLE STATEMENT  CHKXER
      IF( NUMXER(NERR) .NE. INFOT) THEN
         FATAL = .TRUE.
         IF (KPRINT .GE. 3) THEN
           WRITE( NOUT, FMT = 9999 )INFOT, SRNAMT
         ENDIF
      END IF
      RETURN
C
 9999 FORMAT( ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ', I2, ' NOT D',
     $      'ETECTED BY ', A6, ' *****' )
C
C     End of CHKXER.
C
      END
