!** DROTM
SUBROUTINE DROTM(N,Dx,Incx,Dy,Incy,Dparam)
  IMPLICIT NONE
  !>
  !***
  !  Apply a modified Givens transformation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A8
  !***
  ! **Type:**      DOUBLE PRECISION (SROTM-S, DROTM-D)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
  !***
  ! **Author:**  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***
  ! **Description:**
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below.
  !           Locations 2-5 of SPARAM contain elements of the
  !           transformation matrix H described below.
  !
  !     --Output--
  !       DX  rotated vector (unchanged if N .LE. 0)
  !       DY  rotated vector (unchanged if N .LE. 0)
  !
  !     Apply the modified Givens transformation, H, to the 2 by N matrix
  !     (DX**T)
  !     (DY**T), where **T indicates transpose.  The elements of DX are
  !     in DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
  !     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
  !
  !     With DPARAM(1)=DFLAG, H has one of the following forms:
  !
  !     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
  !
  !       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
  !     H=(          )    (          )    (          )    (          )
  !       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
  !
  !     See DROTMG for a description of data storage in DPARAM.
  !
  !***
  ! **References:**  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER i, Incx, Incy, kx, ky, N, nsteps
  REAL(8) :: dflag, dh12, dh22, Dx, two, z, dh11, dh21, &
    Dparam, Dy, w, zero
  DIMENSION Dx(*), Dy(*), Dparam(5)
  SAVE zero, two
  DATA zero, two/0.0D0, 2.0D0/
  !* FIRST EXECUTABLE STATEMENT  DROTM
  dflag = Dparam(1)
  IF ( N>0.AND.(dflag+two/=zero) ) THEN
    IF ( Incx/=Incy.OR.Incx<=0 ) THEN
      kx = 1
      ky = 1
      IF ( Incx<0 ) kx = 1 + (1-N)*Incx
      IF ( Incy<0 ) ky = 1 + (1-N)*Incy
      !
      IF ( dflag<0 ) THEN
        dh11 = Dparam(2)
        dh12 = Dparam(4)
        dh21 = Dparam(3)
        dh22 = Dparam(5)
        DO i = 1, N
          w = Dx(kx)
          z = Dy(ky)
          Dx(kx) = w*dh11 + z*dh12
          Dy(ky) = w*dh21 + z*dh22
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ELSEIF ( dflag==0 ) THEN
        dh12 = Dparam(4)
        dh21 = Dparam(3)
        DO i = 1, N
          w = Dx(kx)
          z = Dy(ky)
          Dx(kx) = w + z*dh12
          Dy(ky) = w*dh21 + z
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ELSE
        dh11 = Dparam(2)
        dh22 = Dparam(5)
        DO i = 1, N
          w = Dx(kx)
          z = Dy(ky)
          Dx(kx) = w*dh11 + z
          Dy(ky) = -w + dh22*z
          kx = kx + Incx
          ky = ky + Incy
        ENDDO
      ENDIF
    ELSE
      !
      nsteps = N*Incx
      IF ( dflag<0 ) THEN
        dh11 = Dparam(2)
        dh12 = Dparam(4)
        dh21 = Dparam(3)
        dh22 = Dparam(5)
        DO i = 1, nsteps, Incx
          w = Dx(i)
          z = Dy(i)
          Dx(i) = w*dh11 + z*dh12
          Dy(i) = w*dh21 + z*dh22
        ENDDO
      ELSEIF ( dflag==0 ) THEN
        dh12 = Dparam(4)
        dh21 = Dparam(3)
        DO i = 1, nsteps, Incx
          w = Dx(i)
          z = Dy(i)
          Dx(i) = w + z*dh12
          Dy(i) = w*dh21 + z
        ENDDO
      ELSE
        dh11 = Dparam(2)
        dh22 = Dparam(5)
        DO i = 1, nsteps, Incx
          w = Dx(i)
          z = Dy(i)
          Dx(i) = w*dh11 + z
          Dy(i) = -w + dh22*z
        ENDDO
      ENDIF
    ENDIF
  ENDIF
END SUBROUTINE DROTM
