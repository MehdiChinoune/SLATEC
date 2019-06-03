!** SMOUT
SUBROUTINE SMOUT(M,N,Lda,A,Ifmt,Idigit)
  !>
  !  Subsidiary to FC and SBOCLS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SMOUT-S, DMOUT-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     SINGLE PRECISION MATRIX OUTPUT ROUTINE.
  !
  !  INPUT..
  !
  !  M,N,LDA,A(*,*) PRINT THE SINGLE PRECISION ARRAY A(I,J),I  = 1,...,M,
  !                 J=1,...,N, ON OUTPUT UNIT LOUT=6. LDA IS THE DECLARED
  !                 FIRST DIMENSION OF A(*,*) AS SPECIFIED IN THE CALLING
  !                 PROGRAM. THE HEADING IN THE FORTRAN FORMAT STATEMENT
  !                 IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST STEP.
  !                 THE COMPONENTS A(I,J) ARE INDEXED, ON OUTPUT, IN A
  !                 PLEASANT FORMAT.
  !  IFMT(*)        A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON
  !                 OUTPUT UNIT LOUT=6 WITH THE VARIABLE FORMAT FORTRAN
  !                 STATEMENT
  !                       WRITE(LOUT,IFMT).
  !  IDIGIT         PRINT AT LEAST ABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
  !                 THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10, OR 14
  !                 WHICH WILL PRINT AT LEAST ABS(IDIGIT) NUMBER OF
  !                 PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE
  !                 UTILIZED TO WRITE EACH LINE OF OUTPUT OF THE ARRAY
  !                 A(*,*). (THIS CAN BE USED ON MOST TIME-SHARING
  !                 TERMINALS).  IF IDIGIT.GE.0, 133 PRINTING COLUMNS ARE
  !                 UTILIZED. (THIS CAN BE USED ON MOST LINE PRINTERS).
  !
  !  EXAMPLE..
  !
  !  PRINT AN ARRAY CALLED (SIMPLEX TABLEAU   ) OF SIZE 10 BY 20 SHOWING
  !  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
  !  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
  !
  !     DIMENSION TABLEU(20,20)
  !     M = 10
  !     N = 20
  !     LDTABL = 20
  !     IDIGIT = -6
  !     CALL SMOUT(M,N,LDTABL,TABLEU,21H(16H1SIMPLEX TABLEAU),IDIGIT)
  !
  !***
  ! **See also:**  FC, SBOCLS
  !***
  ! **Routines called:**  I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891107  Added comma after 1P edit descriptor in FORMAT
  !           statements.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : I1MACH
  INTEGER :: Idigit, Lda, M, N
  REAL(SP) :: A(Lda,N)
  CHARACTER :: Ifmt*(*)
  INTEGER :: i, j, k1, k2, lout, ndigit
  CHARACTER(3), PARAMETER :: icol = 'COL'
  !* FIRST EXECUTABLE STATEMENT  SMOUT
  lout = I1MACH(2)
  WRITE (lout,Ifmt)
  IF ( M<=0.OR.N<=0.OR.Lda<=0 ) RETURN
  ndigit = Idigit
  IF ( Idigit==0 ) ndigit = 4
  IF ( Idigit<0 ) THEN
    !
    ndigit = -Idigit
    IF ( ndigit<=4 ) THEN
      !
      DO k1 = 1, N, 5
        k2 = MIN(N,k1+4)
        WRITE (lout,99001) (icol,i,i=k1,k2)
        DO i = 1, M
          WRITE (lout,99005) i, (A(i,j),j=k1,k2)
        END DO
      END DO
      RETURN
      !
    ELSEIF ( ndigit<=6 ) THEN
      !
      DO k1 = 1, N, 4
        k2 = MIN(N,k1+3)
        WRITE (lout,99002) (icol,i,i=k1,k2)
        DO i = 1, M
          WRITE (lout,99006) i, (A(i,j),j=k1,k2)
        END DO
      END DO
      RETURN
      !
    ELSEIF ( ndigit>10 ) THEN
      !
      DO k1 = 1, N, 2
        k2 = MIN(N,k1+1)
        WRITE (lout,99004) (icol,i,i=k1,k2)
        DO i = 1, M
          WRITE (lout,99008) i, (A(i,j),j=k1,k2)
        END DO
      END DO
      RETURN
    ELSE
      !
      DO k1 = 1, N, 3
        k2 = MIN(N,k1+2)
        WRITE (lout,99003) (icol,i,i=k1,k2)
        DO i = 1, M
          WRITE (lout,99007) i, (A(i,j),j=k1,k2)
        END DO
      END DO
      RETURN
    END IF
    !
  ELSEIF ( ndigit<=4 ) THEN
    !
    DO k1 = 1, N, 10
      k2 = MIN(N,k1+9)
      WRITE (lout,99001) (icol,i,i=k1,k2)
      DO i = 1, M
        WRITE (lout,99005) i, (A(i,j),j=k1,k2)
      END DO
    END DO
    RETURN
    !
  ELSEIF ( ndigit<=6 ) THEN
    !
    DO k1 = 1, N, 8
      k2 = MIN(N,k1+7)
      WRITE (lout,99002) (icol,i,i=k1,k2)
      DO i = 1, M
        WRITE (lout,99006) i, (A(i,j),j=k1,k2)
      END DO
    END DO
    RETURN
    !
  ELSEIF ( ndigit>10 ) THEN
    !
    DO k1 = 1, N, 5
      k2 = MIN(N,k1+4)
      WRITE (lout,99004) (icol,i,i=k1,k2)
      DO i = 1, M
        WRITE (lout,99008) i, (A(i,j),j=k1,k2)
      END DO
    END DO
    RETURN
  END IF
  !
  DO k1 = 1, N, 6
    k2 = MIN(N,k1+5)
    WRITE (lout,99003) (icol,i,i=k1,k2)
    DO i = 1, M
      WRITE (lout,99007) i, (A(i,j),j=k1,k2)
    END DO
  END DO
  RETURN
  99001 FORMAT (10X,10(4X,A,I4,1X))
  99002 FORMAT (10X,8(5X,A,I4,2X))
  99003 FORMAT (10X,6(7X,A,I4,4X))
  99004 FORMAT (10X,5(9X,A,I4,6X))
  99005 FORMAT (1X,'ROW',I4,2X,1P,10E12.3)
  99006 FORMAT (1X,'ROW',I4,2X,1P,8E14.5)
  99007 FORMAT (1X,'ROW',I4,2X,1P,6E18.9)
  99008 FORMAT (1X,'ROW',I4,2X,1P,5E22.13)
END SUBROUTINE SMOUT
