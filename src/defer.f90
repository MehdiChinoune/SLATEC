!*==DEFER.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DEFER
      SUBROUTINE DEFER(COFX,COFY,Idmn,Usol,Grhs)
      IMPLICIT NONE
!*--DEFER5
!*** Start of declarations inserted by SPAG
      REAL ai , AIT , bi , BIT , ci , CIT , DIT , dj , DLX , DLX4 , DLY , DLY4 , 
     &     ej , fj , Grhs , TDLx3 , TDLy3 , tx , ty , Usol
      REAL uxxx , uxxxx , uyyy , uyyyy , xi , yj
      INTEGER i , Idmn , IS , j , JS , K , KSWx , KSWy , L , MIT , MS , NIT , NS
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DEFER
!***SUBSIDIARY
!***PURPOSE  Subsidiary to SEPELI
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DEFER-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine first approximates the truncation error given by
!     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where
!     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and
!     at the boundaries if periodic (here UXXX,UXXXX are the third
!     and fourth partial derivatives of U with respect to X).
!     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
!     at X=A or X=B if the boundary condition there is mixed.
!     TX=0.0 along specified boundaries.  TY has symmetric form
!     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y).
!     The second order solution in USOL is used to approximate
!     (via second order finite differencing) the truncation error
!     and the result is added to the right hand side in GRHS
!     and then transferred to USOL to be used as a new right
!     hand side when calling BLKTRI for a fourth order solution.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  DX, DY
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  DEFER
!
      COMMON /SPLPCM/ KSWx , KSWy , K , L , AIT , BIT , CIT , DIT , MIT , NIT , 
     &                IS , MS , JS , NS , DLX , DLY , TDLx3 , TDLy3 , DLX4 , 
     &                DLY4
      DIMENSION Grhs(Idmn,*) , Usol(Idmn,*)
      EXTERNAL COFX , COFY
!***FIRST EXECUTABLE STATEMENT  DEFER
      DO j = JS , NS
        yj = CIT + (j-1)*DLY
        CALL COFY(yj,dj,ej,fj)
        DO i = IS , MS
          xi = AIT + (i-1)*DLX
          CALL COFX(xi,ai,bi,ci)
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
!
          CALL DX(Usol,Idmn,i,j,uxxx,uxxxx)
          CALL DY(Usol,Idmn,i,j,uyyy,uyyyy)
          tx = ai*uxxxx/12.0 + bi*uxxx/6.0
          ty = dj*uyyyy/12.0 + ej*uyyy/6.0
!
!     RESET FORM OF TRUNCATION IF AT BOUNDARY WHICH IS NON-PERIODIC
!
          IF ( .NOT.(KSWx==1.OR.(i>1.AND.i<K)) )
     &         tx = ai/3.0*(uxxxx/4.0+uxxx/DLX)
          IF ( .NOT.(KSWy==1.OR.(j>1.AND.j<L)) )
     &         ty = dj/3.0*(uyyyy/4.0+uyyy/DLY)
          Grhs(i,j) = Grhs(i,j) + DLX**2*tx + DLY**2*ty
        ENDDO
      ENDDO
!
!     RESET THE RIGHT HAND SIDE IN USOL
!
      DO i = IS , MS
        DO j = JS , NS
          Usol(i,j) = Grhs(i,j)
        ENDDO
      ENDDO
      END SUBROUTINE DEFER
