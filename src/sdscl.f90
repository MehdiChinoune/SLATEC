!*==SDSCL.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SDSCL
      SUBROUTINE SDSCL(Hmax,N,Nq,Rmax,H,Rc,Rh,Yh)
      IMPLICIT NONE
!*--SDSCL5
!***BEGIN PROLOGUE  SDSCL
!***SUBSIDIARY
!***PURPOSE  Subroutine SDSCL rescales the YH array whenever the step
!            size is changed.
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      SINGLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  SDSCL
      INTEGER i , j , N , Nq
      REAL H , Hmax , Rc , Rh , Rmax , r1 , Yh(N,*)
!***FIRST EXECUTABLE STATEMENT  SDSCL
      IF ( H<1.E0 ) THEN
        Rh = MIN(ABS(H)*Rh,ABS(H)*Rmax,Hmax)/ABS(H)
      ELSE
        Rh = MIN(Rh,Rmax,Hmax/ABS(H))
      ENDIF
      r1 = 1.E0
      DO j = 1 , Nq
        r1 = r1*Rh
        DO i = 1 , N
          Yh(i,j+1) = Yh(i,j+1)*r1
        ENDDO
      ENDDO
      H = H*Rh
      Rc = Rc*Rh
      END SUBROUTINE SDSCL
