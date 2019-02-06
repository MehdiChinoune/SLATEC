!*==CSCALE.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CSCALE
      SUBROUTINE CSCALE(A,Nrda,Nrow,Ncol,Cols,Colsav,Rows,Rowsav,Anorm,Scales,
     &                  Iscale,Ic)
      IMPLICIT NONE
!*--CSCALE6
!*** Start of declarations inserted by SPAG
      REAL A , alog2 , Anorm , ascale , Cols , Colsav , cs , p , Rows , Rowsav , 
     &     s , Scales , SDOT , ten20 , ten4
      INTEGER Ic , ip , Iscale , j , k , Ncol , Nrda , Nrow
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CSCALE
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CSCALE-S, DCSCAL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     This routine scales the matrix A by columns when needed
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  CSCALE
      DIMENSION A(Nrda,*) , Cols(*) , Colsav(*) , Scales(*) , Rows(*) , 
     &          Rowsav(*)
!
      SAVE ten4 , ten20
      DATA ten4 , ten20/1.E+4 , 1.E+20/
!
!***FIRST EXECUTABLE STATEMENT  CSCALE
      IF ( Iscale==(-1) ) THEN
!
        IF ( Ic/=0 ) THEN
          DO k = 1 , Ncol
            Cols(k) = SDOT(Nrow,A(1,k),1,A(1,k),1)
          ENDDO
        ENDIF
!
        ascale = Anorm/Ncol
        DO k = 1 , Ncol
          cs = Cols(k)
          IF ( (cs>ten4*ascale).OR.(ten4*cs<ascale) ) GOTO 100
          IF ( (cs<1./ten20).OR.(cs>ten20) ) GOTO 100
        ENDDO
      ENDIF
!
      DO k = 1 , Ncol
        Scales(k) = 1.
      ENDDO
      RETURN
!
 100  alog2 = LOG(2.)
      Anorm = 0.
      DO k = 1 , Ncol
        cs = Cols(k)
        IF ( cs/=0. ) THEN
          p = LOG(cs)/alog2
          ip = -0.5*p
          s = 2.**ip
          Scales(k) = s
          IF ( Ic/=1 ) THEN
            Cols(k) = s*s*Cols(k)
            Anorm = Anorm + Cols(k)
            Colsav(k) = Cols(k)
          ENDIF
          DO j = 1 , Nrow
            A(j,k) = s*A(j,k)
          ENDDO
        ELSE
          Scales(k) = 1.
        ENDIF
      ENDDO
!
      IF ( Ic==0 ) RETURN
!
      DO k = 1 , Nrow
        Rows(k) = SDOT(Ncol,A(k,1),Nrda,A(k,1),Nrda)
        Rowsav(k) = Rows(k)
        Anorm = Anorm + Rows(k)
      ENDDO
      END SUBROUTINE CSCALE
