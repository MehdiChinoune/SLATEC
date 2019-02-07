!*==SROTM.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SROTM
SUBROUTINE SROTM(N,Sx,Incx,Sy,Incy,Sparam)
  IMPLICIT NONE
  !*--SROTM5
  !*** Start of declarations inserted by SPAG
  INTEGER i , Incx , Incy , kx , ky , N , nsteps
  REAL sflag , sh11 , sh12 , sh21 , sh22 , Sparam , Sx , Sy , two , w , z , &
    zero
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SROTM
  !***PURPOSE  Apply a modified Givens transformation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A8
  !***TYPE      SINGLE PRECISION (SROTM-S, DROTM-D)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
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
  !        N  number of elements in input vector(s)
  !       SX  single precision vector with N elements
  !     INCX  storage spacing between elements of SX
  !       SY  single precision vector with N elements
  !     INCY  storage spacing between elements of SY
  !   SPARAM  5-element vector. SPARAM(1) is SFLAG described below.
  !           Locations 2-5 of SPARAM contain elements of the
  !           transformation matrix H described below.
  !
  !     --Output--
  !       SX  rotated vector (unchanged if N .LE. 0)
  !       SY  rotated vector (unchanged if N .LE. 0)
  !
  !     Apply the modified Givens transformation, H, to the 2 by N matrix
  !     (SX**T)
  !     (SY**T) , where **T indicates transpose.  The elements of SX are
  !     in SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
  !     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
  !
  !     With SPARAM(1)=SFLAG, H has one of the following forms:
  !
  !     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
  !
  !       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
  !     H=(          )    (          )    (          )    (          )
  !       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
  !
  !     See SROTMG for a description of data storage in SPARAM.
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
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SROTM
  DIMENSION Sx(*) , Sy(*) , Sparam(5)
  SAVE zero , two
  DATA zero , two/0.0E0 , 2.0E0/
  !***FIRST EXECUTABLE STATEMENT  SROTM
  sflag = Sparam(1)
  IF ( N>0.AND.(sflag+two/=zero) ) THEN
    IF ( Incx/=Incy.OR.Incx<=0 ) THEN
      kx = 1
      ky = 1
      IF ( Incx<0 ) kx = 1 + (1-N)*Incx
      IF ( Incy<0 ) ky = 1 + (1-N)*Incy
      !
      IF ( sflag<0 ) THEN
        sh11 = Sparam(2)
        sh12 = Sparam(4)
        sh21 = Sparam(3)
        sh22 = Sparam(5)
        DO i = 1 , N
          w = Sx(kx)
          z = Sy(ky)
          Sx(kx) = w*sh11 + z*sh12
          Sy(ky) = w*sh21 + z*sh22
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ELSEIF ( sflag==0 ) THEN
        sh12 = Sparam(4)
        sh21 = Sparam(3)
        DO i = 1 , N
          w = Sx(kx)
          z = Sy(ky)
          Sx(kx) = w + z*sh12
          Sy(ky) = w*sh21 + z
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ELSE
        sh11 = Sparam(2)
        sh22 = Sparam(5)
        DO i = 1 , N
          w = Sx(kx)
          z = Sy(ky)
          Sx(kx) = w*sh11 + z
          Sy(ky) = -w + sh22*z
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ENDIF
    ELSE
      !
      nsteps = N*Incx
      IF ( sflag<0 ) THEN
        sh11 = Sparam(2)
        sh12 = Sparam(4)
        sh21 = Sparam(3)
        sh22 = Sparam(5)
        DO i = 1 , nsteps , Incx
          w = Sx(i)
          z = Sy(i)
          Sx(i) = w*sh11 + z*sh12
          Sy(i) = w*sh21 + z*sh22
        ENDDO
      ELSEIF ( sflag==0 ) THEN
        sh12 = Sparam(4)
        sh21 = Sparam(3)
        DO i = 1 , nsteps , Incx
          w = Sx(i)
          z = Sy(i)
          Sx(i) = w + z*sh12
          Sy(i) = w*sh21 + z
        ENDDO
      ELSE
        sh11 = Sparam(2)
        sh22 = Sparam(5)
        DO i = 1 , nsteps , Incx
          w = Sx(i)
          z = Sy(i)
          Sx(i) = w*sh11 + z
          Sy(i) = -w + sh22*z
        ENDDO
      ENDIF
    ENDIF
  ENDIF
END SUBROUTINE SROTM
