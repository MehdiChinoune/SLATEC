!*==DQCGLS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DQCGLS
SUBROUTINE DQCGLS(Lun,Kprint,Ipass)
  !***BEGIN PROLOGUE  DQCGLS
  !***PURPOSE  Quick check for DGLSS.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QCGLSS-S, DQCGLS-D)
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !      QUICK CHECK SUBROUTINE  DQCGLS  TESTS THE EXECUTION
  !      OF THE GENERAL LINEAR SYSTEM SOLVER, DGLSS .  THE
  !      DGLSS  SUBROUTINE PACKAGE WAS WRITTEN BY T. MANTEUFFEL
  !      (LANL).
  !
  !      A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED
  !      BY DQCGLS.  THE SUMMARY LINE GIVES A COUNT OF THE
  !      NUMBER OF PROBLEMS DETECTED DURING THE TEST.
  !
  !      THE REAL QUANTITIES FOR THE COMPUTED SOLUTION VECTOR
  !      X  AND THE CORRESPONDING  RNORM  ARE COMPARED AGAINST
  !      STORED VALUES.  DISAGREEMENT OCCURS IF A DIFFERENCE
  !      IS SQRT(D1MACH(4) OR MORE.  THE RETURNED VALUE (INTEGER)
  !      OF  INFO  IS ALSO CHECKED.  FOUR CASES ARE RUN, TWO
  !      INVOLVING  LLSIA  AND TWO INVOLVING  ULSIA .
  !
  !      DQCGLS REQUIRES NO INPUT ARGUMENTS.  ON RETURN, NERR
  !      (INTEGER TYPE) CONTAINS THE COUNT OF THE NUMBER OF
  !      PROBLEMS DETECTED BY  QCGLSS .
  !
  !***ROUTINES CALLED  D1MACH, DGLSS
  !***REVISION HISTORY  (YYMMDD)
  !   811026  DATE WRITTEN
  !   850601  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
  !           including removing an illegal character from column 1, and
  !           editorial changes.  (RWC)
  !***END PROLOGUE  DQCGLS
  !
  IMPLICIT NONE
  !*--DQCGLS42
  !*** Start of declarations inserted by SPAG
  REAL(8) :: a , aa , b , bb , D1MACH , delmax , delx , r , rnorm ,&
    work , xx
  INTEGER i , Ipass , j , kk , Kprint
  !*** End of declarations inserted by SPAG
  DIMENSION aa(4,4,2) , a(4,4) , bb(4,2) , b(4) , xx(4,4)
  DIMENSION work(50)
  CHARACTER :: list(2)
  INTEGER inf(4) , nerr , kprog , kcase
  INTEGER iwork(20) , info , Lun
  DATA aa/1.D0 , .5D0 , 1.D0 , .25D0 , 0.D0 , 2.D0 , 0.D0 , 1.D0 , 2.D0 ,&
    -1.D0 , 1.D0 , 0.D0 , 0.D0 , 0.D0 , 0.D0 , 0.D0 , 1.D0 , 2.D0 ,&
    -1.D0 , 0.D0 , 0.D0 , 1.D0 , 2.D0 , 0.D0 , -1.D0 , 0.D0 , 1.D0 ,&
    0.D0 , 1.D0 , 0.D0 , 1.D0 , 0.D0/
  DATA bb/3.D0 , 1.5D0 , 2.D0 , 1.25D0 , 1.D0 , 3.D0 , 3.D0 , 0.D0/
  DATA xx/.9999999999999787D0 , 1.000000000000007D0 , 1.000000000000007D0 ,&
    0.D0 , .8095238095238102D0 , 1.047619047619044D0 ,&
    1.095238095238081D0 , 0.D0 , .7777777777777857D0 ,&
    1.444444444444429D0 , .3333333333333393D0 , .5555555555555500D0 ,&
    .3333333333333321D0 , 0.0D0 , -.3333333333333286D0 ,&
    .3333333333333286D0/
  DATA inf/0 , 1 , 0 , 2/
  DATA list/'L' , 'U'/
  !***FIRST EXECUTABLE STATEMENT  DQCGLS
  info = 0
  nerr = 0
  r = MAX(SQRT(D1MACH(4)),1.D-12)
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT (/' *  DQCGLS - QUICK CHECK FOR DGLSS (DLLSIA AND DULSIA)'/)
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
      IF ( kprog==1 ) CALL DGLSS(a,4,4,3,b,4,1,rnorm,work,50,iwork,20,info)
      IF ( kprog==2 ) CALL DGLSS(a,4,3,4,b,4,1,rnorm,work,50,iwork,20,info)
      !
      !           TEST COMPUTED  X , RNORM , AND  INFO .
      !
      kk = 2*(kprog-1) + kcase
      delmax = 0.0D0
      DO i = 1 , 4
        delx = ABS(b(i)-xx(i,kk))
        delmax = MAX(delmax,delx)
      ENDDO
      !
      IF ( Kprint>=3 ) WRITE (Lun,99002) list(kprog) , kcase , delmax
      99002     FORMAT (3X,A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',D11.4/)
      IF ( delmax>=r ) THEN
        nerr = nerr + 1
        IF ( Kprint>=2 ) WRITE (Lun,99003) list(kprog) , kcase , delmax
        99003       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',&
          D11.4/)
      ENDIF
      !
      IF ( Kprint>=3 ) WRITE (Lun,99004) list(kprog) , kcase , rnorm
      99004     FORMAT (3X,A,'LSIA, CASE ',I1,'.  RNORM IS ',D11.4/)
      IF ( rnorm>=r ) THEN
        nerr = nerr + 1
        IF ( Kprint>=2 ) WRITE (Lun,99005) list(kprog) , kcase , rnorm
        99005       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,&
          '.  RNORM (TOO LARGE) IS',D11.4/)
      ENDIF
      IF ( Kprint>=3 ) WRITE (Lun,99006) list(kprog) , kcase , info ,&
        inf(kk)
      !
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
  99008 FORMAT (/' **** DQCGLS DETECTED A TOTAL OF ',I2,&
    ' PROBLEMS WITH DGLSS. ****'/)
  IF ( nerr==0.AND.Kprint>1 ) WRITE (Lun,99009)
  99009 FORMAT ('     DQCGLS DETECTED NO PROBLEMS WITH DGLSS.'/)
  RETURN
END SUBROUTINE DQCGLS
