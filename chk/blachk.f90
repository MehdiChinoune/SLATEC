!*==BLACHK.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK BLACHK
SUBROUTINE BLACHK(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--BLACHK5
  !*** Start of declarations inserted by SPAG
  INTEGER ICAse , INCx , INCy , Kprint , Lun , MODe , N , NPRint
  REAL sdfac , sfac
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BLACHK
  !***PURPOSE  Quick check for Basic Linear Algebra Subprograms.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !     ********************************* TBLA ***************************
  !     TEST DRIVER FOR BASIC LINEAR ALGEBRA SUBPROGRAMS.
  !     C. L. LAWSON, JPL, 1974 DEC 10, 1975 MAY 28
  !
  !     UPDATED BY K. HASKELL - JUNE 23,1980
  !
  !***ROUTINES CALLED  CHECK0, CHECK1, CHECK2, HEADER
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   751210  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  BLACHK
  INTEGER Ipass , jtest(38)
  REAL(8) :: dfac , dqfac
  LOGICAL PASs
  COMMON /COMBLA/ NPRint , ICAse , N , INCx , INCy , MODe , PASs
  DATA sfac , sdfac , dfac , dqfac/.625E-1 , .50 , .625D-1 , 0.625D-1/
  DATA jtest/1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , &
    1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , &
    1 , 1 , 1 , 1 , 1/
  !***FIRST EXECUTABLE STATEMENT  BLACHK
  NPRint = Lun
  Ipass = 1
  !
  IF ( Kprint>=2 ) WRITE (NPRint,99001)
  99001 FORMAT ('1','QUICK CHECK OF 38 BASIC LINEAR ALGEBRA SUBROUTINES'/)
  DO ICAse = 1 , 38
    IF ( jtest(ICAse)/=0 ) THEN
      CALL HEADER(Kprint)
      !
      !         INITIALIZE  PASS, INCX, INCY, AND MODE FOR A NEW CASE.
      !         THE VALUE 9999 FOR INCX, INCY OR MODE WILL APPEAR IN THE
      !         DETAILED  OUTPUT, IF ANY, FOR CASES THAT DO NOT INVOLVE
      !         THESE PARAMETERS.
      !
      PASs = .TRUE.
      INCx = 9999
      INCy = 9999
      MODe = 9999
      SELECT CASE (ICAse)
        CASE (1,2,3,4,5,6,7,8,9,10,11,14,15,18,19,20,21,22,23,24,25)
          !                                       ICASE =  1-11, 14-15, OR 18-25
          CALL CHECK2(sfac,sdfac,dfac,dqfac,Kprint)
        CASE (26,27,28,29,30,31,32,33,34,35,36,37,38)
          !                                       ICASE = 26-38
          CALL CHECK1(sfac,dfac,Kprint)
        CASE DEFAULT
          !                                       ICASE = 12-13 OR 16-17
          CALL CHECK0(sfac,dfac,Kprint)
      END SELECT
      !                                                  PRINT
      IF ( Kprint>=2.AND.PASs ) WRITE (NPRint,99002)
      99002     FORMAT ('+',39X,'PASS')
      IF ( .NOT.PASs ) Ipass = 0
    ENDIF
  ENDDO
  IF ( Kprint>=2.AND.Ipass==1 ) WRITE (NPRint,99003)
  99003 FORMAT (/' ****************BLAS PASSED ALL TESTS****************')
  IF ( Kprint>=1.AND.Ipass==0 ) WRITE (NPRint,99004)
  99004 FORMAT (/' ****************BLAS FAILED SOME TESTS***************')
  RETURN
END SUBROUTINE BLACHK
