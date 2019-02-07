!*==SCHKE3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHKE3
SUBROUTINE SCHKE3(Isnum,Srnamt,Nout,Kprint,Fatal)
  IMPLICIT NONE
  !*--SCHKE35
  !***BEGIN PROLOGUE  SCHKE3
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
  !***ROUTINES CALLED  CHKXER, SGEMM, SSYMM, SSYR2K, SSYRK, STRMM, STRSM,
  !                    XERCLR, XERDMP, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !   930701  Name changed from SCHKE5 to SCHKE3.  (BKS)
  !***END PROLOGUE  SCHKE3
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Isnum, Kprint, Nout
  CHARACTER(6) :: Srnamt
  !     .. Scalars in Common ..
  INTEGER infot
  !     .. Local Scalars ..
  REAL alpha, beta
  INTEGER kontrl
  !     .. Local Arrays ..
  REAL a(1,1), b(1,1), c(1,1)
  !     .. External Subroutines ..
  EXTERNAL CHKXER, SGEMM, SSYMM, SSYR2K, SSYRK, STRMM, STRSM
  !***FIRST EXECUTABLE STATEMENT  SCHKE3
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
      CALL SSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (3)
      infot = 1
      CALL XERCLR
      CALL STRMM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STRMM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STRMM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STRMM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRMM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRMM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (4)
      infot = 1
      CALL XERCLR
      CALL STRSM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STRSM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STRSM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STRSM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STRSM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL STRSM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (5)
      infot = 1
      CALL XERCLR
      CALL SSYRK('/','N',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYRK('U','/',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYRK('U','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYRK('U','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYRK('L','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYRK('L','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SSYRK('U','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SSYRK('U','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SSYRK('L','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SSYRK('L','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (6)
      infot = 1
      CALL XERCLR
      CALL SSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYR2K('U','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL SSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE DEFAULT
      infot = 1
      CALL XERCLR
      CALL SGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 1
      CALL XERCLR
      CALL SGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL SGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL SGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL SGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL SGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
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
  !     End of SCHKE3.
  !
END SUBROUTINE SCHKE3
