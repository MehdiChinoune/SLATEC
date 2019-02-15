!DECK DCHKE3
SUBROUTINE DCHKE3(Isnum,Srnamt,Nout,Kprint,Fatal)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCHKE3
  !***SUBSIDIARY
  !***PURPOSE  Test the error exits from the Level 3 Blas.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Tests the error exits from the Level 3 Blas.
  !  ALPHA, BETA, A, X and Y should not need to be defined.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CHKXER, DGEMM, DSYMM, DSYR2K, DSYRK, DTRMM, DTRSM,
  !                    XERCLR, XERDMP, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  DCHKE3
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Isnum, Nout
  CHARACTER(6) :: Srnamt
  INTEGER infot, Kprint
  !     .. Local Scalars ..
  REAL(8) :: alpha, beta
  INTEGER kontrl
  !     .. Local Arrays ..
  REAL(8) :: a(1,1), b(1,1), c(1,1)
  !     .. External Subroutines ..
  EXTERNAL CHKXER, DGEMM, DSYMM, DTRMM, DTRSM, DSYRK, DSYR2K
  !***FIRST EXECUTABLE STATEMENT  DCHKE3
  CALL XGETF(kontrl)
  IF ( Kprint<=2 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  SELECT CASE (Isnum)
    CASE (2)
      infot = 1
      CALL XERCLR
      CALL DSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (3)
      infot = 1
      CALL XERCLR
      CALL DTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (4)
      infot = 1
      CALL XERCLR
      CALL DTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL DTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL DTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (5)
      infot = 1
      CALL XERCLR
      CALL DSYRK('/','N',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DSYRK('U','/',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYRK('U','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYRK('U','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYRK('L','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYRK('L','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DSYRK('U','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DSYRK('U','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DSYRK('L','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DSYRK('L','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (6)
      infot = 1
      CALL XERCLR
      CALL DSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DSYR2K('U','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL DSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL DSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL DSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE DEFAULT
      infot = 1
      CALL XERCLR
      CALL DGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 1
      CALL XERCLR
      CALL DGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL DGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL DGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL DGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL DGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL DGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL DGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL DGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL DGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL DGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL DGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL DGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL DGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL DGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
  END SELECT
  !
  IF ( Kprint>=2 ) THEN
    CALL XERDMP
    IF ( .NOT.Fatal ) THEN
      WRITE (Nout,FMT=99001) Srnamt
    ELSE
      WRITE (Nout,FMT=99002) Srnamt
    ENDIF
  ENDIF
  CALL XSETF(kontrl)
  RETURN
  !
  99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
  99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****','**')
  !
  !     End of DCHKE3.
  !
END SUBROUTINE DCHKE3
