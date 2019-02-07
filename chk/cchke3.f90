!*==CCHKE3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CCHKE3
SUBROUTINE CCHKE3(Isnum,Srnamt,Nout,Kprint,Fatal)
  IMPLICIT NONE
  !*--CCHKE35
  !***BEGIN PROLOGUE  CCHKE3
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
  !***ROUTINES CALLED  CGEMM, CHEMM, CHER2K, CHERK, CHKXER, CSYMM, CSYR2K,
  !                    CSYRK, CTRMM, CTRSM, XERCLR, XERDMP, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CCHKE3
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Isnum , Kprint , Nout
  CHARACTER(6) :: Srnamt
  !     .. Scalars in Common ..
  INTEGER infot
  !     .. Local Scalars ..
  COMPLEX alpha , beta
  REAL ralpha , rbeta
  INTEGER kontrl
  !     .. Local Arrays ..
  COMPLEX a(2,1) , b(2,1) , c(2,1)
  !     .. External Subroutines ..
  EXTERNAL CGEMM , CHEMM , CHER2K , CHERK , CHKXER , CSYMM , CSYR2K ,&
    CSYRK , CTRMM , CTRSM
  !***FIRST EXECUTABLE STATEMENT  CCHKE3
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
      CALL CHEMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHEMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHEMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHEMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHEMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHEMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHEMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHEMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHEMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHEMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHEMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHEMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHEMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHEMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHEMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHEMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHEMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHEMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHEMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHEMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHEMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHEMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (3)
      infot = 1
      CALL XERCLR
      CALL CSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (4)
      infot = 1
      CALL XERCLR
      CALL CTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','U','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','U','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','L','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','L','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','U','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','U','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','L','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','L','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','U','C','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','U','C','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','L','C','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','L','C','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','U','C','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','U','C','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','L','C','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','L','C','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (5)
      infot = 1
      CALL XERCLR
      CALL CTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','U','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','U','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','L','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','L','C','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','U','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','U','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','L','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','L','C','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','U','C','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','U','C','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','L','C','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','L','C','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','U','C','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','U','C','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','L','C','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','L','C','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (6)
      infot = 1
      CALL XERCLR
      CALL CHERK('/','N',0,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHERK('U','T',0,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHERK('U','N',-1,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHERK('U','C',-1,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHERK('L','N',-1,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHERK('L','C',-1,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHERK('U','N',0,-1,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHERK('U','C',0,-1,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHERK('L','N',0,-1,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHERK('L','C',0,-1,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHERK('U','N',2,0,ralpha,a,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHERK('U','C',0,2,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHERK('L','N',2,0,ralpha,a,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHERK('L','C',0,2,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CHERK('U','N',2,0,ralpha,a,2,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CHERK('U','C',2,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CHERK('L','N',2,0,ralpha,a,2,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CHERK('L','C',2,0,ralpha,a,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (7)
      infot = 1
      CALL XERCLR
      CALL CSYRK('/','N',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CSYRK('U','C',0,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYRK('U','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYRK('U','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYRK('L','N',2,0,alpha,a,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYRK('L','T',0,2,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CSYRK('U','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CSYRK('U','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CSYRK('L','N',2,0,alpha,a,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CSYRK('L','T',2,0,alpha,a,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (8)
      infot = 1
      CALL XERCLR
      CALL CHER2K('/','N',0,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHER2K('U','T',0,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHER2K('U','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHER2K('U','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHER2K('L','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHER2K('L','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHER2K('U','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHER2K('U','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHER2K('L','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CHER2K('L','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER2K('U','N',2,0,alpha,a,1,b,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER2K('U','C',0,2,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER2K('L','N',2,0,alpha,a,1,b,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER2K('L','C',0,2,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHER2K('U','N',2,0,alpha,a,2,b,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHER2K('U','C',0,2,alpha,a,2,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHER2K('L','N',2,0,alpha,a,2,b,1,rbeta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHER2K('L','C',0,2,alpha,a,2,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHER2K('U','N',2,0,alpha,a,2,b,2,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHER2K('U','C',2,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHER2K('L','N',2,0,alpha,a,2,b,2,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CHER2K('L','C',2,0,alpha,a,1,b,1,rbeta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (9)
      infot = 1
      CALL XERCLR
      CALL CSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CSYR2K('U','C',0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 12
      CALL XERCLR
      CALL CSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE DEFAULT
      infot = 1
      CALL XERCLR
      CALL CGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 1
      CALL XERCLR
      CALL CGEMM('/','C',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 1
      CALL XERCLR
      CALL CGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGEMM('C','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('N','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('C','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('C','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('C','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('T','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('N','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('C','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('C','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('C','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('T','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('N','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('C','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('C','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('C','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('T','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('N','C',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('C','N',0,0,2,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('C','C',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('C','T',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('T','C',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('C','N',0,0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('N','C',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('C','C',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('T','C',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('C','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('N','C',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('C','N',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('C','C',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('C','T',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('T','C',2,0,0,alpha,a,1,b,1,beta,c,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
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
  !     End of CCHKE3.
  !
END SUBROUTINE CCHKE3
