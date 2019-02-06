*DECK SCHKE3
      SUBROUTINE SCHKE3 (ISNUM, SRNAMT, NOUT, KPRINT, FATAL)
C***BEGIN PROLOGUE  SCHKE3
C***SUBSIDIARY
C***PURPOSE  Test the error exits from the Level 3 Blas.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Tests the error exits from the Level 3 Blas.
C  ALPHA, BETA, A, X and Y should not need to be defined.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHKXER, SGEMM, SSYMM, SSYR2K, SSYRK, STRMM, STRSM,
C                    XERCLR, XERDMP, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C   930701  Name changed from SCHKE5 to SCHKE3.  (BKS)
C***END PROLOGUE  SCHKE3
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      INTEGER            ISNUM, KPRINT, NOUT
      CHARACTER*6        SRNAMT
C     .. Scalars in Common ..
      INTEGER            INFOT
C     .. Local Scalars ..
      REAL               ALPHA, BETA
      INTEGER            KONTRL
C     .. Local Arrays ..
      REAL               A( 1, 1), B( 1, 1), C( 1, 1)
C     .. External Subroutines ..
      EXTERNAL           CHKXER, SGEMM, SSYMM, SSYR2K, SSYRK, STRMM,
     $                   STRSM
C***FIRST EXECUTABLE STATEMENT  SCHKE3
       CALL XGETF (KONTRL)
       IF (KPRINT .LE. 2) THEN
         CALL XSETF(0)
       ELSE
         CALL XSETF(1)
       ENDIF
      GO TO ( 10, 20, 30, 40, 50, 60 )ISNUM
   10 INFOT = 1
      CALL XERCLR
      CALL SGEMM( '/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 1
      CALL XERCLR
      CALL SGEMM( '/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL SGEMM( 'N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL SGEMM( 'T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SGEMM( 'N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SGEMM( 'N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SGEMM( 'T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SGEMM( 'T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SGEMM( 'N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SGEMM( 'N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SGEMM( 'T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SGEMM( 'T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL SGEMM( 'N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL SGEMM( 'N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL SGEMM( 'T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL SGEMM( 'T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL SGEMM( 'N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL SGEMM( 'N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL SGEMM( 'T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL SGEMM( 'T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SGEMM( 'N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SGEMM( 'T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SGEMM( 'N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SGEMM( 'T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL SGEMM( 'N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL SGEMM( 'N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL SGEMM( 'T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL SGEMM( 'T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   20 INFOT = 1
      CALL XERCLR
      CALL SSYMM( '/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL SSYMM( 'L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYMM( 'L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYMM( 'R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYMM( 'L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYMM( 'R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYMM( 'L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYMM( 'R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYMM( 'L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYMM( 'R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYMM( 'L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYMM( 'R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYMM( 'L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYMM( 'R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYMM( 'L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYMM( 'R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYMM( 'L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYMM( 'R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYMM( 'L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYMM( 'R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYMM( 'L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYMM( 'R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   30 INFOT = 1
      CALL XERCLR
      CALL STRMM( '/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL STRMM( 'L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL STRMM( 'L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRMM( 'R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   40 INFOT = 1
      CALL XERCLR
      CALL STRSM( '/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL STRSM( 'L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL STRSM( 'L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL STRSM( 'R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   50 INFOT = 1
      CALL XERCLR
      CALL SSYRK( '/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL SSYRK( 'U', '/', 0, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYRK( 'U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYRK( 'U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYRK( 'L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYRK( 'L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYRK( 'U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYRK( 'U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYRK( 'L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYRK( 'L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYRK( 'U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYRK( 'U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYRK( 'L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYRK( 'L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SSYRK( 'U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SSYRK( 'U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SSYRK( 'L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL SSYRK( 'L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   60 INFOT = 1
      CALL XERCLR
      CALL SSYR2K( '/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL SSYR2K( 'U', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYR2K( 'U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYR2K( 'U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYR2K( 'L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL SSYR2K( 'L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYR2K( 'U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYR2K( 'U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYR2K( 'L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL SSYR2K( 'L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYR2K( 'U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYR2K( 'U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYR2K( 'L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL SSYR2K( 'L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYR2K( 'U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYR2K( 'U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYR2K( 'L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL SSYR2K( 'L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYR2K( 'U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYR2K( 'U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYR2K( 'L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL SSYR2K( 'L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
C
   70 IF( KPRINT .GE. 2) THEN
        CALL XERDMP
        IF(.NOT.FATAL) THEN
         WRITE( NOUT, FMT = 9999 )SRNAMT
        ELSE
         WRITE( NOUT, FMT = 9998 )SRNAMT
        ENDIF
      END IF
      CALL XSETF (KONTRL)
      RETURN
C
 9999 FORMAT( ' ', A6, ' PASSED THE TESTS OF ERROR-EXITS' )
 9998 FORMAT( ' ******* ', A6, ' FAILED THE TESTS OF ERROR-EXITS *****',
     $      '**' )
C
C     End of SCHKE3.
C
      END
