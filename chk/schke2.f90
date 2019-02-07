!*==SCHKE2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHKE2
SUBROUTINE SCHKE2(Isnum,Srnamt,Nout,Kprint,Fatal)
  IMPLICIT NONE
  !*--SCHKE25
  !***BEGIN PROLOGUE  SCHKE2
  !***SUBSIDIARY
  !***PURPOSE  Test the error exits from the Level 2 Blas.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Tests the error exits from the Level 2 Blas.
  !  ALPHA, BETA, A, X and Y should not need to be defined.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CHKXER, SGBMV, SGEMV, SGER, SSBMV, SSPMV, SSPR,
  !                    SSPR2, SSYMV, SSYR, SSYR2, STBMV, STBSV, STPMV,
  !                    STPSV, STRMV, STRSV, XERCLR, XERDMP, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  SCHKE2
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Isnum , Kprint , Nout
  CHARACTER(6) :: Srnamt
  !     .. Scalars in Common ..
  INTEGER infot
  !     .. Local Scalars ..
  REAL alpha , beta
  INTEGER kontrl
  !     .. Local Arrays ..
  REAL a(1,1) , x(1) , y(1)
  !     .. External Subroutines ..
  EXTERNAL CHKXER , SGBMV , SGEMV , SGER , SSBMV , SSPMV , SSPR , SSPR2 ,&
    SSYMV , SSYR , SSYR2 , STBMV , STBSV , STPMV , STPSV , STRMV ,&
    STRSV
  !***FIRST EXECUTABLE STATEMENT  SCHKE2
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
      CALL SGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL SGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL SGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (3)
      infot = 1
      CALL XERCLR
      CALL SSYMV('/',0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYMV('U',-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SSYMV('U',2,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYMV('U',0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL SSYMV('U',0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (4)
      infot = 1
      CALL XERCLR
      CALL SSBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SSBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL SSBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SSBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL SSBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (5)
      infot = 1
      CALL XERCLR
      CALL SSPMV('/',0,alpha,a,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSPMV('U',-1,alpha,a,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL SSPMV('U',0,alpha,a,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSPMV('U',0,alpha,a,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (6)
      infot = 1
      CALL XERCLR
      CALL STRMV('/','N','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STRMV('U','/','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STRMV('U','N','/',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STRMV('U','N','N',-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRMV('U','N','N',2,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL STRMV('U','N','N',0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (7)
      infot = 1
      CALL XERCLR
      CALL STBMV('/','N','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STBMV('U','/','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STBMV('U','N','/',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STBMV('U','N','N',-1,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STBMV('U','N','N',0,-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL STBMV('U','N','N',0,1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STBMV('U','N','N',0,0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (8)
      infot = 1
      CALL XERCLR
      CALL STPMV('/','N','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STPMV('U','/','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STPMV('U','N','/',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STPMV('U','N','N',-1,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL STPMV('U','N','N',0,a,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (9)
      infot = 1
      CALL XERCLR
      CALL STRSV('/','N','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STRSV('U','/','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STRSV('U','N','/',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STRSV('U','N','N',-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL STRSV('U','N','N',2,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL STRSV('U','N','N',0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (10)
      infot = 1
      CALL XERCLR
      CALL STBSV('/','N','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STBSV('U','/','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STBSV('U','N','/',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STBSV('U','N','N',-1,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL STBSV('U','N','N',0,-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL STBSV('U','N','N',0,1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL STBSV('U','N','N',0,0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (11)
      infot = 1
      CALL XERCLR
      CALL STPSV('/','N','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL STPSV('U','/','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL STPSV('U','N','/',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL STPSV('U','N','N',-1,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL STPSV('U','N','N',0,a,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (12)
      infot = 1
      CALL XERCLR
      CALL SGER(-1,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SGER(0,-1,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SGER(0,0,alpha,x,0,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SGER(0,0,alpha,x,1,y,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SGER(2,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (13)
      infot = 1
      CALL XERCLR
      CALL SSYR('/',0,alpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYR('U',-1,alpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SSYR('U',0,alpha,x,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR('U',2,alpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (14)
      infot = 1
      CALL XERCLR
      CALL SSPR('/',0,alpha,x,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSPR('U',-1,alpha,x,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SSPR('U',0,alpha,x,0,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (15)
      infot = 1
      CALL XERCLR
      CALL SSYR2('/',0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSYR2('U',-1,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SSYR2('U',0,alpha,x,0,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSYR2('U',0,alpha,x,1,y,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL SSYR2('U',2,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (16)
      infot = 1
      CALL XERCLR
      CALL SSPR2('/',0,alpha,x,1,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SSPR2('U',-1,alpha,x,1,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL SSPR2('U',0,alpha,x,0,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL SSPR2('U',0,alpha,x,1,y,0,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE DEFAULT
      infot = 1
      CALL XERCLR
      CALL SGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL SGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL SGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL SGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL SGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL SGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
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
  !     End of SCHKE2.
  !
END SUBROUTINE SCHKE2
