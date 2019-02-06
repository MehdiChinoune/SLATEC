!*==DPLPUP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DPLPUP
      SUBROUTINE DPLPUP(DUSRMT,Mrelas,Nvars,Prgopt,Dattrv,Bl,Bu,Ind,Info,Amat,
     &                  Imat,Sizeup,Asmall,Abig)
      IMPLICIT NONE
!*--DPLPUP6
!*** Start of declarations inserted by SPAG
      INTEGER i , indcat , index , Info , iplace , itcnt , itmax , j , Mrelas , 
     &        Nvars
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DPLPUP
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DSPLP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SPLPUP-S, DPLPUP-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
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
!     FROM THE USER.  IT IS PART OF THE DSPLP( ) PACKAGE.
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  DPCHNG, DPNNZR, XERMSG
!***REVISION HISTORY  (YYMMDD)
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
!***END PROLOGUE  DPLPUP
      DOUBLE PRECISION Abig , aij , Amat(*) , amn , amx , Asmall , Bl(*) , 
     &                 Bu(*) , Dattrv(*) , Prgopt(*) , xval , zero
      INTEGER iflag(10) , Imat(*) , Ind(*)
      LOGICAL Sizeup , first
      CHARACTER*8 xern1 , xern2
      CHARACTER*16 xern3 , xern4
!
!***FIRST EXECUTABLE STATEMENT  DPLPUP
      zero = 0.D0
!
!     CHECK USER-SUPPLIED BOUNDS
!
!     CHECK THAT IND(*) VALUES ARE 1,2,3 OR 4.
!     ALSO CHECK CONSISTENCY OF UPPER AND LOWER BOUNDS.
!
      DO j = 1 , Nvars
        IF ( Ind(j)<1.OR.Ind(j)>4 ) THEN
          WRITE (xern1,'(I8)') j
          CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, INDEPENDENT VARIABLE = '//
     &                xern1//' IS NOT DEFINED.',10,1)
          Info = -10
          RETURN
        ENDIF
!
        IF ( Ind(j)==3 ) THEN
          IF ( Bl(j)>Bu(j) ) THEN
            WRITE (xern1,'(I8)') j
            WRITE (xern3,'(1PE15.6)') Bl(j)
            WRITE (xern4,'(1PE15.6)') Bu(j)
            CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, LOWER BOUND = '//xern3//
     &                  ' AND UPPER BOUND = '//xern4//
     &                  ' FOR INDEPENDENT VARIABLE = '//xern1//
     &                  ' ARE NOT CONSISTENT.',11,1)
            RETURN
          ENDIF
        ENDIF
      ENDDO
!
      DO i = Nvars + 1 , Nvars + Mrelas
        IF ( Ind(i)<1.OR.Ind(i)>4 ) THEN
          WRITE (xern1,'(I8)') i - Nvars
          CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, DEPENDENT VARIABLE = '//
     &                xern1//' IS NOT DEFINED.',12,1)
          Info = -12
          RETURN
        ENDIF
!
        IF ( Ind(i)==3 ) THEN
          IF ( Bl(i)>Bu(i) ) THEN
            WRITE (xern1,'(I8)') i
            WRITE (xern3,'(1PE15.6)') Bl(i)
            WRITE (xern4,'(1PE15.6)') Bu(i)
            CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, LOWER BOUND = '//xern3//
     &                  ' AND UPPER BOUND = '//xern4//
     &                  ' FOR DEPENDANT VARIABLE = '//xern1//
     &                  ' ARE NOT CONSISTENT.',13,1)
            Info = -13
            RETURN
          ENDIF
        ENDIF
      ENDDO
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
        IF ( itcnt>itmax ) THEN
          CALL XERMSG('SLATEC','DPLPUP',
     &                'IN DSPLP, MORE THAN 2*NVARS*MRELAS ITERATIONS DEFINING '
     &                //'OR UPDATING MATRIX DATA.',7,1)
          Info = -7
          RETURN
        ENDIF
!
        aij = zero
        CALL DUSRMT(i,j,aij,indcat,Prgopt,Dattrv,iflag)
        IF ( iflag(1)==1 ) THEN
          iflag(1) = 2
          CYCLE
        ENDIF
!
!     CHECK TO SEE THAT THE SUBSCRIPTS I AND J ARE VALID.
!
        IF ( i<1.OR.i>Mrelas.OR.j<1.OR.j>Nvars ) THEN
!
!        CHECK ON SIZE OF MATRIX DATA
!        RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
!
          IF ( iflag(1)==3 ) THEN
            IF ( Sizeup.AND.ABS(aij)/=zero ) THEN
              IF ( first ) THEN
                amx = ABS(aij)
                amn = ABS(aij)
                first = .FALSE.
              ELSEIF ( ABS(aij)>amx ) THEN
                amx = ABS(aij)
              ELSEIF ( ABS(aij)<amn ) THEN
                amn = ABS(aij)
              ENDIF
            ENDIF
            EXIT
          ENDIF
!
          WRITE (xern1,'(I8)') i
          WRITE (xern2,'(I8)') j
          CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, ROW INDEX = '//xern1//
     &                ' OR COLUMN INDEX = '//xern2//' IS OUT OF RANGE.',8,1)
          Info = -8
          RETURN
        ENDIF
!
!     IF INDCAT=0 THEN SET A(I,J)=AIJ.
!     IF INDCAT=1 THEN ACCUMULATE ELEMENT, A(I,J)=A(I,J)+AIJ.
!
        IF ( indcat==0 ) THEN
          CALL DPCHNG(i,aij,iplace,Amat,Imat,j)
        ELSEIF ( indcat==1 ) THEN
          index = -(i-1)
          CALL DPNNZR(index,xval,iplace,Amat,Imat,j)
          IF ( index==i ) aij = aij + xval
          CALL DPCHNG(i,aij,iplace,Amat,Imat,j)
        ELSE
          WRITE (xern1,'(I8)') indcat
          CALL XERMSG('SLATEC','DPLPUP','IN DSPLP, INDICATION FLAG = '//xern1//
     &                ' FOR MATRIX DATA MUST BE EITHER 0 OR 1.',9,1)
          Info = -9
          RETURN
        ENDIF
!
!     CHECK ON SIZE OF MATRIX DATA
!     RECORD THE LARGEST AND SMALLEST(IN MAGNITUDE) NONZERO ELEMENTS.
!
        IF ( Sizeup.AND.ABS(aij)/=zero ) THEN
          IF ( first ) THEN
            amx = ABS(aij)
            amn = ABS(aij)
            first = .FALSE.
          ELSEIF ( ABS(aij)>amx ) THEN
            amx = ABS(aij)
          ELSEIF ( ABS(aij)<amn ) THEN
            amn = ABS(aij)
          ENDIF
        ENDIF
        IF ( iflag(1)==3 ) EXIT
      ENDDO
!
      IF ( Sizeup.AND..NOT.first ) THEN
        IF ( amn<Asmall.OR.amx>Abig ) THEN
          CALL XERMSG('SLATEC','DPLPUP',
     &                'IN DSPLP, A MATRIX ELEMENT''S SIZE IS OUT OF THE '//
     &                'SPECIFIED RANGE.',22,1)
          Info = -22
          RETURN
        ENDIF
      ENDIF
      END SUBROUTINE DPLPUP
