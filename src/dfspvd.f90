!*==DFSPVD.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFSPVD
SUBROUTINE DFSPVD(T,K,X,Ileft,Vnikx,Nderiv)
  !***BEGIN PROLOGUE  DFSPVD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DFC
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BSPLVD-S, DFSPVD-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   **** Double Precision Version of BSPLVD ****
  ! Calculates value and deriv.s of all B-splines which do not vanish at X
  !
  !  Fill VNIKX(J,IDERIV), J=IDERIV, ... ,K  with nonzero values of
  !  B-splines of order K+1-IDERIV , IDERIV=NDERIV, ... ,1, by repeated
  !  calls to DFSPVN
  !
  !***SEE ALSO  DFC
  !***ROUTINES CALLED  DFSPVN
  !***REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DFSPVD
  IMPLICIT NONE
  !*--DFSPVD30
  !*** Start of declarations inserted by SPAG
  REAL(8) :: a , diff , fkmd , T , v , Vnikx , X
  INTEGER i , ideriv , idervm , Ileft , ipkmd , j , jlow , jm1 , K , kmd , &
    l , m , Nderiv
  !*** End of declarations inserted by SPAG
  DIMENSION T(*) , Vnikx(K,*)
  DIMENSION a(20,20)
  !***FIRST EXECUTABLE STATEMENT  DFSPVD
  CALL DFSPVN(T,K+1-Nderiv,1,X,Ileft,Vnikx(Nderiv,Nderiv))
  IF ( Nderiv>1 ) THEN
    ideriv = Nderiv
    DO i = 2 , Nderiv
      idervm = ideriv - 1
      DO j = ideriv , K
        Vnikx(j-1,idervm) = Vnikx(j,ideriv)
      ENDDO
      ideriv = idervm
      CALL DFSPVN(T,0,2,X,Ileft,Vnikx(ideriv,ideriv))
    ENDDO
    !
    DO i = 1 , K
      DO j = 1 , K
        a(i,j) = 0.D0
      ENDDO
      a(i,i) = 1.D0
    ENDDO
    kmd = K
    DO m = 2 , Nderiv
      kmd = kmd - 1
      fkmd = kmd
      i = Ileft
      j = K
      DO
        jm1 = j - 1
        ipkmd = i + kmd
        diff = T(ipkmd) - T(i)
        IF ( jm1==0 ) THEN
          IF ( diff/=0. ) a(1,1) = a(1,1)/diff*fkmd
          !
          DO i = 1 , K
            v = 0.D0
            jlow = MAX(i,m)
            DO j = jlow , K
              v = a(i,j)*Vnikx(j,m) + v
            ENDDO
            Vnikx(i,m) = v
          ENDDO
          EXIT
        ELSE
          IF ( diff/=0.D0 ) THEN
            DO l = 1 , j
              a(l,j) = (a(l,j)-a(l,j-1))/diff*fkmd
            ENDDO
          ENDIF
          j = jm1
          i = i - 1
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE DFSPVD
