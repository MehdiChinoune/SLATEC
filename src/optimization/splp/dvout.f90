!** DVOUT
SUBROUTINE DVOUT(N,Dx,Ifmt,Idigit)
  !>
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SVOUT-S, DVOUT-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION VECTOR OUTPUT ROUTINE.
  !
  !  INPUT..
  !
  !  N,DX(*) PRINT THE DOUBLE PRECISION ARRAY DX(I),I=1,...,N, ON
  !          OUTPUT UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT
  !          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
  !          STEP. THE COMPONENTS DX(I) ARE INDEXED, ON OUTPUT,
  !          IN A PLEASANT FORMAT.
  !  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
  !          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT
  !                WRITE(LOUT,IFMT)
  !  IDIGIT  PRINT AT LEAST ABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
  !          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10 OR 14
  !          WHICH WILL PRINT AT LEAST ABS(IDIGIT) NUMBER OF
  !          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
  !          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY DX(*). (THIS
  !          CAN BE USED ON MOST TIME-SHARING TERMINALS). IF
  !          IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN
  !          BE USED ON MOST LINE PRINTERS).
  !
  !  EXAMPLE..
  !
  !  PRINT AN ARRAY CALLED (COSTS OF PURCHASES) OF LENGTH 100 SHOWING
  !  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
  !  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
  !
  !     DOUBLE PRECISION COSTS(100)
  !     N = 100
  !     IDIGIT = -6
  !     CALL DVOUT(N,COSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891107  Added comma after 1P edit descriptor in FORMAT
  !           statements.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR section.  (WRB)
  USE service, ONLY : I1MACH
  INTEGER :: Idigit, N
  REAL(8) :: Dx(N)
  CHARACTER :: Ifmt*(*)
  INTEGER :: i, k1, k2, lout, ndigit
  !* FIRST EXECUTABLE STATEMENT  DVOUT
  lout = I1MACH(2)
  WRITE (lout,Ifmt)
  IF ( N<=0 ) RETURN
  ndigit = Idigit
  IF ( Idigit==0 ) ndigit = 6
  IF ( Idigit<0 ) THEN
    !
    ndigit = -Idigit
    IF ( ndigit<=6 ) THEN
      !
      DO k1 = 1, N, 4
        k2 = MIN(N,k1+3)
        WRITE (lout,99001) k1, k2, (Dx(i),i=k1,k2)
      END DO
      RETURN
      !
    ELSEIF ( ndigit<=14 ) THEN
      !
      DO k1 = 1, N, 2
        k2 = MIN(N,k1+1)
        WRITE (lout,99002) k1, k2, (Dx(i),i=k1,k2)
      END DO
      RETURN
      !
    ELSEIF ( ndigit>20 ) THEN
      !
      DO k1 = 1, N
        k2 = k1
        WRITE (lout,99004) k1, k2, (Dx(i),i=k1,k2)
      END DO
      RETURN
    ELSE
      !
      DO k1 = 1, N, 2
        k2 = MIN(N,k1+1)
        WRITE (lout,99003) k1, k2, (Dx(i),i=k1,k2)
      END DO
      RETURN
    END IF
    !
  ELSEIF ( ndigit<=6 ) THEN
    !
    DO k1 = 1, N, 8
      k2 = MIN(N,k1+7)
      WRITE (lout,99001) k1, k2, (Dx(i),i=k1,k2)
    END DO
    RETURN
    !
  ELSEIF ( ndigit<=14 ) THEN
    !
    DO k1 = 1, N, 5
      k2 = MIN(N,k1+4)
      WRITE (lout,99002) k1, k2, (Dx(i),i=k1,k2)
    END DO
    RETURN
    !
  ELSEIF ( ndigit>20 ) THEN
    !
    DO k1 = 1, N, 3
      k2 = MIN(N,k1+2)
      WRITE (lout,99004) k1, k2, (Dx(i),i=k1,k2)
    END DO
    RETURN
  END IF
  !
  DO k1 = 1, N, 4
    k2 = MIN(N,k1+3)
    WRITE (lout,99003) k1, k2, (Dx(i),i=k1,k2)
  END DO
  RETURN
  99001 FORMAT (1X,I4,' - ',I4,1X,1P,8D14.5)
  99002 FORMAT (1X,I4,' - ',I4,1X,1P,5D22.13)
  99003 FORMAT (1X,I4,' - ',I4,1X,1P,4D28.19)
  99004 FORMAT (1X,I4,' - ',I4,1X,1P,3D36.27)
END SUBROUTINE DVOUT
