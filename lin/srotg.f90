!DECK SROTG
SUBROUTINE SROTG(Sa,Sb,Sc,Ss)
  IMPLICIT NONE
  REAL r, Sa, Sb, Sc, Ss, u, v
  !***BEGIN PROLOGUE  SROTG
  !***PURPOSE  Construct a plane Givens rotation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B10
  !***TYPE      SINGLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
  !***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
  !             LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !       SA  single precision scalar
  !       SB  single precision scalar
  !
  !     --Output--
  !       SA  single precision result R
  !       SB  single precision result Z
  !       SC  single precision result
  !       SS  single precision result
  !
  !     Construct the Givens transformation
  !
  !         ( SC  SS )
  !     G = (        ),    SC**2 + SS**2 = 1 ,
  !         (-SS  SC )
  !
  !     which zeros the second entry of the 2-vector  (SA,SB)**T.
  !
  !     The quantity R = (+/-)SQRT(SA**2 + SB**2) overwrites SA in
  !     storage.  The value of SB is overwritten by a value Z which
  !     allows SC and SS to be recovered by the following algorithm:
  !
  !           If Z=1  set  SC=0.0  and  SS=1.0
  !           If ABS(Z) .LT. 1  set  SC=SQRT(1-Z**2)  and  SS=Z
  !           If ABS(Z) .GT. 1  set  SC=1/Z  and  SS=SQRT(1-SC**2)
  !
  !     Normally, the subprogram SROT(N,SX,INCX,SY,INCY,SC,SS) will
  !     next be called to apply the transformation to a 2 by N matrix.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SROTG
  !***FIRST EXECUTABLE STATEMENT  SROTG
  IF ( ABS(Sa)>ABS(Sb) ) THEN
    !
    ! *** HERE ABS(SA) .GT. ABS(SB) ***
    !
    u = Sa + Sa
    v = Sb/u
    !
    !     NOTE THAT U AND R HAVE THE SIGN OF SA
    !
    r = SQRT(0.25E0+v**2)*u
    !
    !     NOTE THAT SC IS POSITIVE
    !
    Sc = Sa/r
    Ss = v*(Sc+Sc)
    Sb = Ss
    Sa = r
    RETURN
    !
    ! *** HERE ABS(SA) .LE. ABS(SB) ***
    !
  ELSEIF ( Sb==0.0E0 ) THEN
    !
    ! *** HERE SA = SB = 0.0 ***
    !
    Sc = 1.0E0
    Ss = 0.0E0
    RETURN
  ELSE
    u = Sb + Sb
    v = Sa/u
    !
    !     NOTE THAT U AND R HAVE THE SIGN OF SB
    !     (R IS IMMEDIATELY STORED IN SA)
    !
    Sa = SQRT(0.25E0+v**2)*u
    !
    !     NOTE THAT SS IS POSITIVE
    !
    Ss = Sb/Sa
    Sc = v*(Ss+Ss)
    IF ( Sc/=0.0E0 ) THEN
      Sb = 1.0E0/Sc
      RETURN
    ENDIF
  ENDIF
  Sb = 1.0E0
  RETURN
END SUBROUTINE SROTG
