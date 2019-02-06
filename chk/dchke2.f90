!*==DCHKE2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DCHKE2
      SUBROUTINE DCHKE2(Isnum,Srnamt,Nout,Kprint,Fatal)
      IMPLICIT NONE
!*--DCHKE25
!***BEGIN PROLOGUE  DCHKE2
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
!***ROUTINES CALLED  CHKXER, DGBMV, DGEMV, DGER, DSBMV, DSPMV, DSPR,
!                    DSPR2, DSYMV, DSYR, DSYR2, DTBMV, DTBSV, DTPMV,
!                    DTPSV, DTRMV, DTRSV, XERCLR, XERDMP, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   870810  DATE WRITTEN
!   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  DCHKE2
!     .. Scalar Arguments ..
      INTEGER Isnum , Nout
      LOGICAL Fatal
      CHARACTER*6 Srnamt
      INTEGER infot , Kprint
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , beta
      INTEGER kontrl
!     .. Local Arrays ..
      DOUBLE PRECISION a(1,1) , x(1) , y(1)
!     .. External Subroutines ..
      EXTERNAL CHKXER , DGBMV , DGEMV , DGER , DSBMV , DSPMV , DSPR , DSPR2 , 
     &         DSYMV , DSYR , DSYR2 , DTBMV , DTBSV , DTPMV , DTPSV , DTRMV , 
     &         DTRSV
!***FIRST EXECUTABLE STATEMENT  DCHKE2
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
        CALL DGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 13
        CALL XERCLR
        CALL DGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (3)
        infot = 1
        CALL XERCLR
        CALL DSYMV('/',0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYMV('U',-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYMV('U',2,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYMV('U',0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 10
        CALL XERCLR
        CALL DSYMV('U',0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (4)
        infot = 1
        CALL XERCLR
        CALL DSBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DSBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DSBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DSBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DSBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (5)
        infot = 1
        CALL XERCLR
        CALL DSPMV('/',0,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPMV('U',-1,alpha,a,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DSPMV('U',0,alpha,a,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSPMV('U',0,alpha,a,x,1,beta,y,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (6)
        infot = 1
        CALL XERCLR
        CALL DTRMV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRMV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRMV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRMV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRMV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DTRMV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (7)
        infot = 1
        CALL XERCLR
        CALL DTBMV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTBMV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTBMV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTBMV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTBMV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTBMV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTBMV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (8)
        infot = 1
        CALL XERCLR
        CALL DTPMV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTPMV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTPMV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTPMV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTPMV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (9)
        infot = 1
        CALL XERCLR
        CALL DTRSV('/','N','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTRSV('U','/','N',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTRSV('U','N','/',0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTRSV('U','N','N',-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DTRSV('U','N','N',2,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DTRSV('U','N','N',0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (10)
        infot = 1
        CALL XERCLR
        CALL DTBSV('/','N','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTBSV('U','/','N',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTBSV('U','N','/',0,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTBSV('U','N','N',-1,0,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DTBSV('U','N','N',0,-1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTBSV('U','N','N',0,1,a,1,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DTBSV('U','N','N',0,0,a,1,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (11)
        infot = 1
        CALL XERCLR
        CALL DTPSV('/','N','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DTPSV('U','/','N',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DTPSV('U','N','/',0,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 4
        CALL XERCLR
        CALL DTPSV('U','N','N',-1,a,x,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DTPSV('U','N','N',0,a,x,0)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (12)
        infot = 1
        CALL XERCLR
        CALL DGER(-1,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGER(0,-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DGER(0,0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DGER(0,0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DGER(2,0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (13)
        infot = 1
        CALL XERCLR
        CALL DSYR('/',0,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYR('U',-1,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYR('U',0,alpha,x,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR('U',2,alpha,x,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (14)
        infot = 1
        CALL XERCLR
        CALL DSPR('/',0,alpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPR('U',-1,alpha,x,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSPR('U',0,alpha,x,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (15)
        infot = 1
        CALL XERCLR
        CALL DSYR2('/',0,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSYR2('U',-1,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSYR2('U',0,alpha,x,0,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSYR2('U',0,alpha,x,1,y,0,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 9
        CALL XERCLR
        CALL DSYR2('U',2,alpha,x,1,y,1,a,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE (16)
        infot = 1
        CALL XERCLR
        CALL DSPR2('/',0,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DSPR2('U',-1,alpha,x,1,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 5
        CALL XERCLR
        CALL DSPR2('U',0,alpha,x,0,y,1,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 7
        CALL XERCLR
        CALL DSPR2('U',0,alpha,x,1,y,0,a)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
      CASE DEFAULT
        infot = 1
        CALL XERCLR
        CALL DGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 2
        CALL XERCLR
        CALL DGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 3
        CALL XERCLR
        CALL DGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 6
        CALL XERCLR
        CALL DGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 8
        CALL XERCLR
        CALL DGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
        CALL CHKXER(Srnamt,infot,Nout,Fatal,Kprint)
        infot = 11
        CALL XERCLR
        CALL DGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
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
!     End of DCHKE2.
!
      END SUBROUTINE DCHKE2
