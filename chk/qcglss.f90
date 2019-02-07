!*==QCGLSS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QCGLSS
SUBROUTINE QCGLSS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QCGLSS5
  !*** Start of declarations inserted by SPAG
  INTEGER i , Ipass , j , kk , Kprint
  REAL R1MACH , rnorm
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  QCGLSS
  !***PURPOSE  Quick check for SGLSS.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QCGLSS-S, DQCGLS-D)
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !      QUICK CHECK SUBROUTINE  QCGLSS  TESTS THE EXECUTION
  !      OF THE GENERAL LINEAR SYSTEM SOLVER, SGLSS .  THE
  !      SGLSS  SUBROUTINE PACKAGE WAS WRITTEN BY T. MANTEUFFEL
  !      (LANL).
  !
  !      A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED
  !      BY QCGLSS.  THE SUMMARY LINE GIVES A COUNT OF THE
  !      NUMBER OF PROBLEMS DETECTED DURING THE TEST.
  !
  !      THE REAL QUANTITIES FOR THE COMPUTED SOLUTION VECTOR
  !      X  AND THE CORRESPONDING  RNORM  ARE COMPARED AGAINST
  !      STORED VALUES.  DISAGREEMENT OCCURS IF A DIFFERENCE
  !      IS SQRT(R1MACH(4) OR MORE.  THE RETURNED VALUE (INTEGER)
  !      OF  INFO  IS ALSO CHECKED.  FOUR CASES ARE RUN, TWO
  !      INVOLVING  LLSIA  AND TWO INVOLVING  ULSIA .
  !
  !      QCGLSS REQUIRES NO INPUT ARGUMENTS.  ON RETURN, NERR
  !      (INTEGER TYPE) CONTAINS THE COUNT OF THE NUMBER OF
  !      PROBLEMS DETECTED BY  QCGLSS .
  !
  !***ROUTINES CALLED  R1MACH, SGLSS
  !***REVISION HISTORY  (YYMMDD)
  !   811026  DATE WRITTEN
  !   820801  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
  !           including removing an illegal character from column 1, and
  !           editorial changes.  (RWC)
  !***END PROLOGUE  QCGLSS
  REAL aa(4,4,2) , a(4,4) , bb(4,2) , b(4) , xx(4,4) , delmax , delx , r
  REAL work(20)
  CHARACTER :: list(2)
  INTEGER inf(4) , nerr , kprog , kcase
  INTEGER iwork(7) , info , Lun
  DATA aa/1. , .5 , 1. , .25 , 0. , 2. , 0. , 1. , 2. , -1. , 1. , 0. , 0. ,&
    0. , 0. , 0. , 1. , 2. , -1. , 0. , 0. , 1. , 2. , 0. , -1. , 0. ,&
    1. , 0. , 1. , 0. , 1. , 0./
  DATA bb/3. , 1.5 , 2. , 1.25 , 1. , 3. , 3. , 0./
  DATA xx/.9999999999999787 , 1.000000000000007 , 1.000000000000007 , 0. ,&
    .8095238095238102 , 1.047619047619044 , 1.095238095238081 , 0. ,&
    .7777777777777857 , 1.444444444444429 , .3333333333333393 ,&
    .5555555555555500 , .3333333333333321 , 0.0 , -.3333333333333286 ,&
    .3333333333333286/
  DATA inf/0 , 1 , 0 , 2/
  DATA list/'L' , 'U'/
  !***FIRST EXECUTABLE STATEMENT  QCGLSS
  info = 0
  nerr = 0
  r = SQRT(R1MACH(4))
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (/' *    QCGLSS - QUICK CHECK FOR SGLSS (LLSIA AND ULSIA)'/)
  DO kprog = 1 , 2
    DO kcase = 1 , 2
      !
      !           FORM BASIC MATRIX  A  AND VECTOR  B .  (CASE 1)
      !
      DO i = 1 , 4
        DO j = 1 , 4
          a(i,j) = aa(i,j,kprog)
        ENDDO
        b(i) = bb(i,kprog)
      ENDDO
      !
      !           MAKE 3 ROWS IDENTICAL FOR CASE 2.
      !
      IF ( kcase/=1 ) THEN
        DO i = 2 , 3
          DO j = 1 , 4
            a(i,j) = a(1,j)
          ENDDO
          b(i) = b(1)
        ENDDO
      ENDIF
      !
      !           SOLVE FOR VECTOR  X .
      !
      info = 0
      IF ( kprog==1 ) CALL SGLSS(a,4,4,3,b,4,1,rnorm,work,20,iwork,7,info)
      IF ( kprog==2 ) CALL SGLSS(a,4,3,4,b,4,1,rnorm,work,20,iwork,7,info)
      !
      !           TEST COMPUTED  X , RNORM , AND  INFO .
      !
      kk = 2*(kprog-1) + kcase
      delmax = 0.0E0
      DO i = 1 , 4
        delx = ABS(b(i)-xx(i,kk))
        delmax = MAX(delmax,delx)
      ENDDO
      !
      IF ( Kprint>=3 ) WRITE (Lun,99002) list(kprog) , kcase , delmax
      !
      99002     FORMAT (3X,A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',E11.4/)
      IF ( delmax>=r ) THEN
        nerr = nerr + 1
        IF ( Kprint>=2 ) WRITE (Lun,99003) list(kprog) , kcase , delmax
        99003       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',&
          E11.4/)
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99004) list(kprog) , kcase , rnorm
      99004     FORMAT (3X,A,'LSIA, CASE ',I1,'.  RNORM IS ',E11.4/)
      IF ( rnorm>r ) THEN
        nerr = nerr + 1
        IF ( Kprint>=2 ) WRITE (Lun,99005) list(kprog) , kcase , rnorm
        99005       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,&
          '.  RNORM (TOO LARGE) IS',E11.4/)
      ENDIF
      !
      IF ( Kprint>=3 ) WRITE (Lun,99006) list(kprog) , kcase , info ,&
        inf(kk)
      99006     FORMAT (3X,A,'LSIA, CASE ',I1,'.  INFO=',I1,' (SHOULD = ',I1,')'/)
      IF ( info/=inf(kk) ) THEN
        nerr = nerr + 1
        IF ( Kprint>=2 ) WRITE (Lun,99007) list(kprog) , kcase , info ,&
          inf(kk)
        99007       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  INFO=',I1,&
          ' (SHOULD = ',I1,')'/)
      ENDIF
    ENDDO
  ENDDO
  !
  !     SUMMARY PRINT
  !
  Ipass = 0
  IF ( nerr==0 ) Ipass = 1
  IF ( nerr/=0.AND.Kprint/=0 ) WRITE (Lun,99008) nerr
  99008 FORMAT (/' **** QCGLSS DETECTED A TOTAL OF ',I2,&
    ' PROBLEMS WITH SGLSS. ****'/)
  IF ( nerr==0.AND.Kprint>1 ) WRITE (Lun,99009)
  99009 FORMAT ('     QCGLSS DETECTED NO PROBLEMS WITH SGLSS.'/)
  RETURN
END SUBROUTINE QCGLSS
