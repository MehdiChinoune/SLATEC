!** SPLPUP
PURE SUBROUTINE SPLPUP(USRMAT,Mrelas,Nvars,Dattrv,Bl,Bu,Ind,Info,Amat,&
    Imat,Sizeup,Asmall,Abig)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SPLPUP-S, DPLPUP-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/.
  !
  !     REVISED 810613-1130
  !     REVISED YYMMDD-HHMM
  !
  !     THIS SUBROUTINE COLLECTS INFORMATION ABOUT THE BOUNDS AND MATRIX
  !     FROM THE USER.  IT IS PART OF THE SPLP( ) PACKAGE.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  PCHNGS, PNNZRS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Corrected references to XERRWV.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls, changed do-it-yourself
  !           DO loops to DO loops.  (RWC)
  !   900602  Get rid of ASSIGNed GOTOs.  (RWC)

  INTERFACE
    PURE SUBROUTINE USRMAT(I,J,Aij,Indcat,Dattrv,Iflag)
      IMPORT SP
      INTEGER, INTENT(OUT) :: I, J, Indcat
      INTEGER, INTENT(INOUT) :: Iflag(4)
      REAL(SP), INTENT(IN) :: Dattrv(:)
      REAL(SP), INTENT(OUT) :: Aij
    END SUBROUTINE USRMAT
  END INTERFACE
  INTEGER, INTENT(IN) :: Mrelas, Nvars
  INTEGER, INTENT(OUT) :: Info
  REAL(SP), INTENT(IN) :: Abig, Asmall
  LOGICAL, INTENT(IN) :: Sizeup
  INTEGER, INTENT(IN) :: Ind(Nvars+Mrelas)
  INTEGER, INTENT(INOUT) :: Imat(:)
  REAL(SP), INTENT(IN) :: Bl(Nvars+Mrelas), Bu(Nvars+Mrelas), Dattrv(:)
  REAL(SP), INTENT(INOUT) :: Amat(:)
  INTEGER :: i, indcat, indexx, iplace, itcnt, itmax, j
  INTEGER :: iflag(10)
  LOGICAL :: first
  REAL(SP) :: aij, amn, amx, xval
  CHARACTER(8) :: xern1, xern2
  CHARACTER(16) :: xern3, xern4
  !
  !* FIRST EXECUTABLE STATEMENT  SPLPUP
  !
  !     CHECK USER-SUPPLIED BOUNDS
  !
  !     CHECK THAT IND(*) VALUES ARE 1,2,3 OR 4.
  !     ALSO CHECK CONSISTENCY OF UPPER AND LOWER BOUNDS.
  !
  DO j = 1, Nvars
    IF( Ind(j)<1 .OR. Ind(j)>4 ) THEN
      WRITE (xern1,'(I8)') j
      Info = -10
      ERROR STOP 'SPLPUP : IN SPLP, INDEPENDENT VARIABLE IS NOT DEFINED.'
      RETURN
    END IF
    !
    IF( Ind(j)==3 ) THEN
      IF( Bl(j)>Bu(j) ) THEN
        WRITE (xern1,'(I8)') j
        WRITE (xern3,'(1PE15.6)') Bl(j)
        WRITE (xern4,'(1PE15.6)') Bu(j)
        ERROR STOP 'SPLPUP : IN SPLP, LOWER BOUND  AND UPPER BOUND = &
          & FOR INDEPENDENT VARIABLE ARE NOT CONSISTENT.'
        RETURN
      END IF
    END IF
  END DO
  !
  DO i = Nvars + 1, Nvars + Mrelas
    IF( Ind(i)<1 .OR. Ind(i)>4 ) THEN
      WRITE (xern1,'(I8)') i - Nvars
      Info = -12
      ERROR STOP 'SPLPUP : IN SPLP, DEPENDENT VARIABLE IS NOT DEFINED.'
      RETURN
    END IF
    !
    IF( Ind(i)==3 ) THEN
      IF( Bl(i)>Bu(i) ) THEN
        WRITE (xern1,'(I8)') i
        WRITE (xern3,'(1PE15.6)') Bl(i)
        WRITE (xern4,'(1PE15.6)') Bu(i)
        Info = -13
        ERROR STOP 'SPLPUP : IN SPLP, LOWER BOUND AND UPPER BOUND &
          &FOR DEPENDANT VARIABLE  ARE NOT CONSISTENT.'
        RETURN
      END IF
    END IF
  END DO
  !
  !     GET UPDATES OR DATA FOR MATRIX FROM THE USER
  !
  !     GET THE ELEMENTS OF THE MATRIX FROM THE USER.  IT WILL BE STORED
  !     BY COLUMNS USING THE SPARSE STORAGE CODES OF RJ HANSON AND
  !     JA WISNIEWSKI.
  !
  iflag(1) = 1
  !
  !     KEEP ACCEPTING ELEMENTS UNTIL THE USER IS FINISHED GIVING THEM.
  !     LIMIT THIS LOOP TO 2*NVARS*MRELAS ITERATIONS.
  !
  itmax = 2*Nvars*Mrelas + 1
  itcnt = 0
  first = .TRUE.
  DO
    !
    !     CHECK ON THE ITERATION COUNT.
    !
    itcnt = itcnt + 1
    IF( itcnt>itmax ) THEN
      Info = -7
      ERROR STOP 'SPLPUP : IN SPLP, MORE THAN 2*NVARS*MRELAS ITERATIONS DEFINING OR UPDATING MATRIX DATA.'
      RETURN
    END IF
    !
    aij = 0._SP
    CALL USRMAT(i,j,aij,indcat,Dattrv,iflag)
    IF( iflag(1)==1 ) THEN
      iflag(1) = 2
      CYCLE
    END IF
    !
    !     CHECK TO SEE THAT THE SUBSCRIPTS I AND J ARE VALID.
    !
    IF( i<1 .OR. i>Mrelas .OR. j<1 .OR. j>Nvars ) THEN
      !
      !        CHECK ON SIZE OF MATRIX DATA
      !        RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
      !
      IF( iflag(1)==3 ) THEN
        IF( Sizeup .AND. ABS(aij)/=0._SP ) THEN
          IF( first ) THEN
            amx = ABS(aij)
            amn = ABS(aij)
            first = .FALSE.
          ELSEIF( ABS(aij)>amx ) THEN
            amx = ABS(aij)
          ELSEIF( ABS(aij)<amn ) THEN
            amn = ABS(aij)
          END IF
        END IF
        EXIT
      END IF
      !
      WRITE (xern1,'(I8)') i
      WRITE (xern2,'(I8)') j
      Info = -8
      ERROR STOP 'SPLPUP : IN SPLP, ROW INDEX OR COLUMN INDEX  IS OUT OF RANGE.'
      RETURN
    END IF
    !
    !     IF INDCAT=0 THEN SET A(I,J)=AIJ.
    !     IF INDCAT=1 THEN ACCUMULATE ELEMENT, A(I,J)=A(I,J)+AIJ.
    !
    IF( indcat==0 ) THEN
      CALL PCHNGS(i,aij,iplace,Amat,Imat,j)
    ELSEIF( indcat==1 ) THEN
      indexx = -(i-1)
      CALL PNNZRS(indexx,xval,iplace,Amat,Imat,j)
      IF( indexx==i ) aij = aij + xval
      CALL PCHNGS(i,aij,iplace,Amat,Imat,j)
    ELSE
      WRITE (xern1,'(I8)') indcat
      ERROR STOP 'SPLPUP : IN SPLP, INDICATION FLAG FOR MATRIX DATA MUST BE EITHER 0 OR 1.'
      Info = -9
      RETURN
    END IF
    !
    !     CHECK ON SIZE OF MATRIX DATA
    !     RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
    !
    IF( Sizeup .AND. ABS(aij)/=0._SP ) THEN
      IF( first ) THEN
        amx = ABS(aij)
        amn = ABS(aij)
        first = .FALSE.
      ELSEIF( ABS(aij)>amx ) THEN
        amx = ABS(aij)
      ELSEIF( ABS(aij)<amn ) THEN
        amn = ABS(aij)
      END IF
    END IF
    IF( iflag(1)==3 ) EXIT
  END DO
  !
  IF( Sizeup .AND. .NOT. first ) THEN
    IF( amn<Asmall .OR. amx>Abig ) THEN
      Info = -22
      ERROR STOP 'SPLPUP : IN SPLP, A MATRIX ELEMENT''S SIZE IS OUT OF THE SPECIFIED RANGE.'
      RETURN
    END IF
  END IF
END SUBROUTINE SPLPUP
