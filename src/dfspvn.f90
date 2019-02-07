!*==DFSPVN.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFSPVN
SUBROUTINE DFSPVN(T,Jhigh,Index,X,Ileft,Vnikx)
  !***BEGIN PROLOGUE  DFSPVN
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DFC
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BSPLVN-S, DFSPVN-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !  **** Double Precision version of BSPLVN ****
  !
  ! Calculates the value of all possibly nonzero B-splines at *X* of
  !  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*.
  !
  !***SEE ALSO  DFC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DFSPVN
  IMPLICIT NONE
  !*--DFSPVN25
  !*** Start of declarations inserted by SPAG
  REAL(8) :: deltam , deltap , T , vm , vmprev , Vnikx , X
  INTEGER i , Ileft , imjp1 , Index , ipj , j , Jhigh , jp1 , jp1ml , l
  !*** End of declarations inserted by SPAG
  DIMENSION T(*) , Vnikx(*)
  DIMENSION deltam(20) , deltap(20)
  SAVE j , deltam , deltap
  DATA j/1/ , (deltam(i),i=1,20) , (deltap(i),i=1,20)/40*0.0D0/
  !***FIRST EXECUTABLE STATEMENT  DFSPVN
  IF ( Index/=2 ) THEN
    j = 1
    Vnikx(1) = 1.D0
    IF ( j>=Jhigh ) GOTO 99999
  ENDIF
  DO
    !
    ipj = Ileft + j
    deltap(j) = T(ipj) - X
    imjp1 = Ileft - j + 1
    deltam(j) = X - T(imjp1)
    vmprev = 0.D0
    jp1 = j + 1
    DO l = 1 , j
      jp1ml = jp1 - l
      vm = Vnikx(l)/(deltap(l)+deltam(jp1ml))
      Vnikx(l) = vm*deltap(l) + vmprev
      vmprev = vm*deltam(jp1ml)
    ENDDO
    Vnikx(jp1) = vmprev
    j = jp1
    IF ( j>=Jhigh ) EXIT
  ENDDO
  !
  99999 END SUBROUTINE DFSPVN
