!*==IVOUT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK IVOUT
SUBROUTINE IVOUT(N,Ix,Ifmt,Idigit)
  IMPLICIT NONE
  !*--IVOUT5
  !*** Start of declarations inserted by SPAG
  INTEGER i , I1MACH , Idigit , Ix , j , k1 , k2 , lout , N , ndigit
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  IVOUT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (IVOUT-I)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***DESCRIPTION
  !
  !     INTEGER VECTOR OUTPUT ROUTINE.
  !
  !  INPUT..
  !
  !  N,IX(*) PRINT THE INTEGER ARRAY IX(I),I=1,...,N, ON OUTPUT
  !          UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT
  !          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
  !          STEP. THE COMPONENTS IX(I) ARE INDEXED, ON OUTPUT,
  !          IN A PLEASANT FORMAT.
  !  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
  !          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT
  !                WRITE(LOUT,IFMT)
  !  IDIGIT  PRINT UP TO ABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
  !          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10 OR 14
  !          WHICH WILL PRINT AT LEAST ABS(IDIGIT) NUMBER OF
  !          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
  !          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY IX(*). (THIS
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
  !     DIMENSION ICOSTS(100)
  !     N = 100
  !     IDIGIT = -6
  !     CALL IVOUT(N,ICOSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  I1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910403  Updated AUTHOR section.  (WRB)
  !***END PROLOGUE  IVOUT
  DIMENSION Ix(*)
  CHARACTER Ifmt*(*)
  !
  !     GET THE UNIT NUMBER WHERE OUTPUT WILL BE WRITTEN.
  !***FIRST EXECUTABLE STATEMENT  IVOUT
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
      DO k1 = 1 , N , 10
        k2 = MIN(N,k1+9)
        WRITE (lout,99001) k1 , k2 , (Ix(i),i=k1,k2)
      ENDDO
      RETURN
      !
    ELSEIF ( ndigit<=6 ) THEN
      !
      DO k1 = 1 , N , 7
        k2 = MIN(N,k1+6)
        WRITE (lout,99002) k1 , k2 , (Ix(i),i=k1,k2)
      ENDDO
      RETURN
      !
    ELSEIF ( ndigit>10 ) THEN
      !
      DO k1 = 1 , N , 3
        k2 = MIN(N,k1+2)
        WRITE (lout,99004) k1 , k2 , (Ix(i),i=k1,k2)
      ENDDO
      RETURN
    ELSE
      !
      DO k1 = 1 , N , 5
        k2 = MIN(N,k1+4)
        WRITE (lout,99003) k1 , k2 , (Ix(i),i=k1,k2)
      ENDDO
      RETURN
    ENDIF
    !
  ELSEIF ( ndigit<=4 ) THEN
    !
    DO k1 = 1 , N , 20
      k2 = MIN(N,k1+19)
      WRITE (lout,99001) k1 , k2 , (Ix(i),i=k1,k2)
    ENDDO
    RETURN
    !
  ELSEIF ( ndigit<=6 ) THEN
    !
    DO k1 = 1 , N , 15
      k2 = MIN(N,k1+14)
      WRITE (lout,99002) k1 , k2 , (Ix(i),i=k1,k2)
    ENDDO
    RETURN
    !
  ELSEIF ( ndigit>10 ) THEN
    !
    DO k1 = 1 , N , 7
      k2 = MIN(N,k1+6)
      WRITE (lout,99004) k1 , k2 , (Ix(i),i=k1,k2)
    ENDDO
    RETURN
  ENDIF
  !
  DO k1 = 1 , N , 10
    k2 = MIN(N,k1+9)
    WRITE (lout,99003) k1 , k2 , (Ix(i),i=k1,k2)
  ENDDO
  RETURN
  99001 FORMAT (1X,I4,' - ',I4,20(1X,I5))
  99002 FORMAT (1X,I4,' - ',I4,15(1X,I7))
  99003 FORMAT (1X,I4,' - ',I4,10(1X,I11))
  99004 FORMAT (1X,I4,' - ',I4,7(1X,I15))
END SUBROUTINE IVOUT
