!*==DUSRMT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DUSRMT
      SUBROUTINE DUSRMT(I,J,Aij,Indcat,Prgopt,Dattrv,Iflag)
      IMPLICIT NONE
!*--DUSRMT5
!*** Start of declarations inserted by SPAG
      INTEGER I , Indcat , J , l
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DUSRMT
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DSPLP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (USRMAT-S, DUSRMT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   The user may supply this code
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DUSRMT
      DOUBLE PRECISION Prgopt(*) , Dattrv(*) , Aij
      INTEGER Iflag(*)
!
!***FIRST EXECUTABLE STATEMENT  DUSRMT
      IF ( Iflag(1)==1 ) THEN
!
!     THIS IS THE INITIALIZATION STEP.  THE VALUES OF IFLAG(K),K=2,3,4,
!     ARE RESPECTIVELY THE COLUMN INDEX, THE ROW INDEX (OR THE NEXT COL.
!     INDEX), AND THE POINTER TO THE MATRIX ENTRY'S VALUE WITHIN
!     DATTRV(*).  ALSO CHECK (DATTRV(1)=0.) SIGNIFYING NO DATA.
        IF ( Dattrv(1)==0.D0 ) THEN
          I = 0
          J = 0
          Iflag(1) = 3
        ELSE
          Iflag(2) = -Dattrv(1)
          Iflag(3) = Dattrv(2)
          Iflag(4) = 3
        ENDIF
!
        RETURN
      ELSE
        J = Iflag(2)
        I = Iflag(3)
        l = Iflag(4)
        IF ( I==0 ) THEN
!
!     SIGNAL THAT ALL OF THE NONZERO ENTRIES HAVE BEEN DEFINED.
          Iflag(1) = 3
          RETURN
        ELSEIF ( I<0 ) THEN
!
!     SIGNAL THAT A SWITCH IS MADE TO A NEW COLUMN.
          J = -I
          I = Dattrv(l)
          l = l + 1
        ENDIF
!
        Aij = Dattrv(l)
!
!     UPDATE THE INDICES AND POINTERS FOR THE NEXT ENTRY.
        Iflag(2) = J
        Iflag(3) = Dattrv(l+1)
        Iflag(4) = l + 2
!
!     INDCAT=0 DENOTES THAT ENTRIES OF THE MATRIX ARE ASSIGNED THE
!     VALUES FROM DATTRV(*).  NO ACCUMULATION IS PERFORMED.
        Indcat = 0
        RETURN
      ENDIF
      END SUBROUTINE DUSRMT
