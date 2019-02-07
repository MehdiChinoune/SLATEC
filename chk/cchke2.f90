!*==CCHKE2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CCHKE2
SUBROUTINE CCHKE2(Isnum,Srnamt,Nout,Kprint,Fatal)
  IMPLICIT NONE
  !*--CCHKE25
  !***BEGIN PROLOGUE  CCHKE2
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
  !***ROUTINES CALLED  CGBMV, CGEMV, CGERC, CGERU, CHBMV, CHEMV, CHER,
  !                    CHER2, CHKXER, CHPMV, CHPR, CHPR2, CTBMV, CTBSV,
  !                    CTPMV, CTPSV, CTRMV, CTRSV, XERCLR, XERDMP, XGETF,
  !                    XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CCHKE2
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Isnum , Kprint , Nout
  CHARACTER(6) :: Srnamt
  !     .. Scalars in Common ..
  INTEGER infot
  !     .. Local Scalars ..
  COMPLEX alpha , beta
  REAL ralpha
  INTEGER kontrl
  !     .. Local Arrays ..
  COMPLEX a(1,1) , x(1) , y(1)
  !     .. External Subroutines ..
  EXTERNAL CGBMV , CGEMV , CGERC , CGERU , CHBMV , CHEMV , CHER , CHER2 ,&
    CHKXER , CHPMV , CHPR , CHPR2 , CTBMV , CTBSV , CTPMV , CTPSV ,&
    CTRMV , CTRSV
  !***FIRST EXECUTABLE STATEMENT  CCHKE2
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
      CALL CGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 13
      CALL XERCLR
      CALL CGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (3)
      infot = 1
      CALL XERCLR
      CALL CHEMV('/',0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHEMV('U',-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CHEMV('U',2,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHEMV('U',0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 10
      CALL XERCLR
      CALL CHEMV('U',0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (4)
      infot = 1
      CALL XERCLR
      CALL CHBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CHBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CHBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CHBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CHBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (5)
      infot = 1
      CALL XERCLR
      CALL CHPMV('/',0,alpha,a,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHPMV('U',-1,alpha,a,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CHPMV('U',0,alpha,a,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHPMV('U',0,alpha,a,x,1,beta,y,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (6)
      infot = 1
      CALL XERCLR
      CALL CTRMV('/','N','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTRMV('U','/','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTRMV('U','N','/',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTRMV('U','N','N',-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRMV('U','N','N',2,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CTRMV('U','N','N',0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (7)
      infot = 1
      CALL XERCLR
      CALL CTBMV('/','N','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTBMV('U','/','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTBMV('U','N','/',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTBMV('U','N','N',-1,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTBMV('U','N','N',0,-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CTBMV('U','N','N',0,1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTBMV('U','N','N',0,0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (8)
      infot = 1
      CALL XERCLR
      CALL CTPMV('/','N','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTPMV('U','/','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTPMV('U','N','/',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTPMV('U','N','N',-1,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CTPMV('U','N','N',0,a,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (9)
      infot = 1
      CALL XERCLR
      CALL CTRSV('/','N','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTRSV('U','/','N',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTRSV('U','N','/',0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTRSV('U','N','N',-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CTRSV('U','N','N',2,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CTRSV('U','N','N',0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (10)
      infot = 1
      CALL XERCLR
      CALL CTBSV('/','N','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTBSV('U','/','N',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTBSV('U','N','/',0,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTBSV('U','N','N',-1,0,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CTBSV('U','N','N',0,-1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CTBSV('U','N','N',0,1,a,1,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CTBSV('U','N','N',0,0,a,1,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (11)
      infot = 1
      CALL XERCLR
      CALL CTPSV('/','N','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CTPSV('U','/','N',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CTPSV('U','N','/',0,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 4
      CALL XERCLR
      CALL CTPSV('U','N','N',-1,a,x,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CTPSV('U','N','N',0,a,x,0)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (12)
      infot = 1
      CALL XERCLR
      CALL CGERC(-1,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGERC(0,-1,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGERC(0,0,alpha,x,0,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CGERC(0,0,alpha,x,1,y,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CGERC(2,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (13)
      infot = 1
      CALL XERCLR
      CALL CGERU(-1,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGERU(0,-1,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CGERU(0,0,alpha,x,0,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CGERU(0,0,alpha,x,1,y,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CGERU(2,0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (14)
      infot = 1
      CALL XERCLR
      CALL CHER('/',0,ralpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHER('U',-1,ralpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CHER('U',0,ralpha,x,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER('U',2,ralpha,x,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (15)
      infot = 1
      CALL XERCLR
      CALL CHPR('/',0,ralpha,x,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHPR('U',-1,ralpha,x,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CHPR('U',0,ralpha,x,0,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (16)
      infot = 1
      CALL XERCLR
      CALL CHER2('/',0,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHER2('U',-1,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CHER2('U',0,alpha,x,0,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHER2('U',0,alpha,x,1,y,0,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 9
      CALL XERCLR
      CALL CHER2('U',2,alpha,x,1,y,1,a,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE (17)
      infot = 1
      CALL XERCLR
      CALL CHPR2('/',0,alpha,x,1,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CHPR2('U',-1,alpha,x,1,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 5
      CALL XERCLR
      CALL CHPR2('U',0,alpha,x,0,y,1,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 7
      CALL XERCLR
      CALL CHPR2('U',0,alpha,x,1,y,0,a)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
    CASE DEFAULT
      infot = 1
      CALL XERCLR
      CALL CGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 2
      CALL XERCLR
      CALL CGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 3
      CALL XERCLR
      CALL CGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 6
      CALL XERCLR
      CALL CGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 8
      CALL XERCLR
      CALL CGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
      CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      infot = 11
      CALL XERCLR
      CALL CGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
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
  !     End of CCHKE2.
  !
END SUBROUTINE CCHKE2
