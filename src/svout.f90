!DECK SVOUT
SUBROUTINE SVOUT(N,Sx,Ifmt,Idigit)
  IMPLICIT NONE
  INTEGER i, I1MACH, Idigit, j, k1, k2, lout, N, ndigit
  REAL Sx
  !***BEGIN PROLOGUE  SVOUT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SVOUT-S, DVOUT-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     SINGLE PRECISION VECTOR OUTPUT ROUTINE.
  !
  !  INPUT..
  !
  !  N,SX(*) PRINT THE SINGLE PRECISION ARRAY SX(I),I=1,...,N, ON
  !          OUTPUT UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT
  !          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
  !          STEP. THE COMPONENTS SX(I) ARE INDEXED, ON OUTPUT,
  !          IN A PLEASANT FORMAT.
  !  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
  !          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT
  !                WRITE(LOUT,IFMT)
  !  IDIGIT  PRINT AT LEAST ABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
  !          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10 OR 14
  !          WHICH WILL PRINT AT LEAST ABS(IDIGIT) NUMBER OF
  !          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
  !          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY SX(*). (THIS
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
  !     DIMENSION COSTS(100)
  !     N = 100
  !     IDIGIT = -6
  !     CALL SVOUT(N,COSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  I1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891107  Added comma after 1P edit descriptor in FORMAT
  !           statements.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  SVOUT
  DIMENSION Sx(*)
  CHARACTER Ifmt*(*)
  !
  !     GET THE UNIT NUMBER WHERE OUTPUT WILL BE WRITTEN.
  !***FIRST EXECUTABLE STATEMENT  SVOUT
  j = 2
  lout = I1MACH(j)
  WRITE (lout,Ifmt)
  IF ( N<=0 ) RETURN
  ndigit = Idigit
  IF ( Idigit==0 ) ndigit = 4
  IF ( Idigit<0 ) THEN
    !
    ndigit = -Idigit
    IF ( ndigit<=4 ) THEN
      !
      DO k1 = 1, N, 5
        k2 = MIN(N,k1+4)
        WRITE (lout,99001) k1, k2, (Sx(i),i=k1,k2)
      ENDDO
      RETURN
      !
    ELSEIF ( ndigit<=6 ) THEN
      !
      DO k1 = 1, N, 4
        k2 = MIN(N,k1+3)
        WRITE (lout,99002) k1, k2, (Sx(i),i=k1,k2)
      ENDDO
      RETURN
      !
    ELSEIF ( ndigit>10 ) THEN
      !
      DO k1 = 1, N, 2
        k2 = MIN(N,k1+1)
        WRITE (lout,99004) k1, k2, (Sx(i),i=k1,k2)
      ENDDO
      RETURN
    ELSE
      !
      DO k1 = 1, N, 3
        k2 = MIN(N,k1+2)
        WRITE (lout,99003) k1, k2, (Sx(i),i=k1,k2)
      ENDDO
      RETURN
    ENDIF
    !
  ELSEIF ( ndigit<=4 ) THEN
    !
    DO k1 = 1, N, 10
      k2 = MIN(N,k1+9)
      WRITE (lout,99001) k1, k2, (Sx(i),i=k1,k2)
    ENDDO
    RETURN
    !
  ELSEIF ( ndigit<=6 ) THEN
    !
    DO k1 = 1, N, 8
      k2 = MIN(N,k1+7)
      WRITE (lout,99002) k1, k2, (Sx(i),i=k1,k2)
    ENDDO
    RETURN
    !
  ELSEIF ( ndigit>10 ) THEN
    !
    DO k1 = 1, N, 5
      k2 = MIN(N,k1+4)
      WRITE (lout,99004) k1, k2, (Sx(i),i=k1,k2)
    ENDDO
    RETURN
  ENDIF
  !
  DO k1 = 1, N, 6
    k2 = MIN(N,k1+5)
    WRITE (lout,99003) k1, k2, (Sx(i),i=k1,k2)
  ENDDO
  RETURN
  99001 FORMAT (1X,I4,' - ',I4,1P,10E12.3)
  99002 FORMAT (1X,I4,' - ',I4,1X,1P,8E14.5)
  99003 FORMAT (1X,I4,' - ',I4,1X,1P,6E18.9)
  99004 FORMAT (1X,I4,' - ',I4,1X,1P,5E24.13)
END SUBROUTINE SVOUT
