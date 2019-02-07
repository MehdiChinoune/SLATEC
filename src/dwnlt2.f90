!*==DWNLT2.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DWNLT2
LOGICAL FUNCTION DWNLT2(Me,Mend,Ir,Factor,Tau,Scale,Wic)
  IMPLICIT NONE
  !*--DWNLT25
  !***BEGIN PROLOGUE  DWNLT2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to WNLIT
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
  !
  !     To test independence of incoming column.
  !
  !     Test the column IC to determine if it is linearly independent
  !     of the columns already in the basis.  In the initial tri. step,
  !     we usually want the heavy weight ALAMDA to be included in the
  !     test for independence.  In this case, the value of FACTOR will
  !     have been set to 1.E0 before this procedure is invoked.
  !     In the potentially rank deficient problem, the value of FACTOR
  !     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
  !     heavy weight from the test for independence.
  !
  !     Write new column as partitioned vector
  !           (A1)  number of components in solution so far = NIV
  !           (A2)  M-NIV components
  !     And compute  SN = inverse weighted length of A1
  !                  RN = inverse weighted length of A2
  !     Call the column independent when RN .GT. TAU*SN
  !
  !***SEE ALSO  DWNLIT
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
  !   900604  DP version created from SP version.  (RWC)
  !***END PROLOGUE  DWNLT2
  DOUBLE PRECISION Factor , Scale(*) , Tau , Wic(*)
  INTEGER Ir , Me , Mend
  !
  DOUBLE PRECISION rn , sn , t
  INTEGER j
  !
  !***FIRST EXECUTABLE STATEMENT  DWNLT2
  sn = 0.E0
  rn = 0.E0
  DO j = 1 , Mend
    t = Scale(j)
    IF ( j<=Me ) t = t/Factor
    t = t*Wic(j)**2
    !
    IF ( j<Ir ) THEN
      sn = sn + t
    ELSE
      rn = rn + t
    ENDIF
  ENDDO
  DWNLT2 = rn>sn*Tau**2
END FUNCTION DWNLT2
