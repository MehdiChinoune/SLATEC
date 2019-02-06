*DECK DCHKE3
      SUBROUTINE DCHKE3 (ISNUM, SRNAMT, NOUT, KPRINT, FATAL)
C***BEGIN PROLOGUE  DCHKE3
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
C***ROUTINES CALLED  CHKXER, DGEMM, DSYMM, DSYR2K, DSYRK, DTRMM, DTRSM,
C                    XERCLR, XERDMP, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  DCHKE3
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      INTEGER            ISNUM, NOUT
      CHARACTER*6        SRNAMT
      INTEGER            INFOT, KPRINT
C     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            KONTRL
C     .. Local Arrays ..
      DOUBLE PRECISION   A( 1, 1), B( 1, 1), C( 1, 1)
C     .. External Subroutines ..
      EXTERNAL           CHKXER, DGEMM, DSYMM, DTRMM, DTRSM, DSYRK,
     $                   DSYR2K
C***FIRST EXECUTABLE STATEMENT  DCHKE3
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
        CALL XSETF(0)
      ELSE
        CALL XSETF(1)
      ENDIF
      GO TO ( 10, 20, 30, 40, 50, 60 )ISNUM
   10 INFOT = 1
      CALL XERCLR
      CALL DGEMM( '/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 1
      CALL XERCLR
      CALL DGEMM( '/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DGEMM( 'N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DGEMM( 'T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DGEMM( 'N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DGEMM( 'N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DGEMM( 'T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DGEMM( 'T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DGEMM( 'N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DGEMM( 'N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DGEMM( 'T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DGEMM( 'T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DGEMM( 'N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DGEMM( 'N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DGEMM( 'T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DGEMM( 'T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL DGEMM( 'N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL DGEMM( 'N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL DGEMM( 'T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 8
      CALL XERCLR
      CALL DGEMM( 'T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DGEMM( 'N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DGEMM( 'T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DGEMM( 'N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DGEMM( 'T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL DGEMM( 'N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL DGEMM( 'N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL DGEMM( 'T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 13
      CALL XERCLR
      CALL DGEMM( 'T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   20 INFOT = 1
      CALL XERCLR
      CALL DSYMM( '/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DSYMM( 'L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYMM( 'L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYMM( 'R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYMM( 'L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYMM( 'R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYMM( 'L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYMM( 'R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYMM( 'L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYMM( 'R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYMM( 'L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYMM( 'R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYMM( 'L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYMM( 'R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYMM( 'L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYMM( 'R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYMM( 'L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYMM( 'R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYMM( 'L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYMM( 'R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYMM( 'L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYMM( 'R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   30 INFOT = 1
      CALL XERCLR
      CALL DTRMM( '/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DTRMM( 'L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DTRMM( 'L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRMM( 'R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   40 INFOT = 1
      CALL XERCLR
      CALL DTRSM( '/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DTRSM( 'L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DTRSM( 'L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 5
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 6
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 11
      CALL XERCLR
      CALL DTRSM( 'R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   50 INFOT = 1
      CALL XERCLR
      CALL DSYRK( '/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DSYRK( 'U', '/', 0, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYRK( 'U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYRK( 'U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYRK( 'L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYRK( 'L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYRK( 'U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYRK( 'U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYRK( 'L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYRK( 'L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYRK( 'U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYRK( 'U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYRK( 'L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYRK( 'L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DSYRK( 'U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DSYRK( 'U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DSYRK( 'L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 10
      CALL XERCLR
      CALL DSYRK( 'L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      GO TO 70
   60 INFOT = 1
      CALL XERCLR
      CALL DSYR2K( '/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 2
      CALL XERCLR
      CALL DSYR2K( 'U', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYR2K( 'U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYR2K( 'U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYR2K( 'L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 3
      CALL XERCLR
      CALL DSYR2K( 'L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYR2K( 'U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYR2K( 'U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYR2K( 'L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 4
      CALL XERCLR
      CALL DSYR2K( 'L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYR2K( 'U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYR2K( 'U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYR2K( 'L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 7
      CALL XERCLR
      CALL DSYR2K( 'L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYR2K( 'U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYR2K( 'U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYR2K( 'L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 9
      CALL XERCLR
      CALL DSYR2K( 'L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYR2K( 'U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYR2K( 'U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYR2K( 'L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
      INFOT = 12
      CALL XERCLR
      CALL DSYR2K( 'L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1 )
      CALL CHKXER( SRNAMT, INFOT, NOUT, FATAL, KPRINT )
C
   70 IF( KPRINT .GE. 2) THEN
        CALL XERDMP
        IF( .NOT.FATAL )THEN
         WRITE( NOUT, FMT = 9999 )SRNAMT
        ELSE
         WRITE( NOUT, FMT = 9998 )SRNAMT
        END IF
      ENDIF
      CALL XSETF (KONTRL)
      RETURN
C
 9999 FORMAT( ' ', A6, ' PASSED THE TESTS OF ERROR-EXITS' )
 9998 FORMAT( ' ******* ', A6, ' FAILED THE TESTS OF ERROR-EXITS *****',
     $      '**' )
C
C     End of DCHKE3.
C
      END
