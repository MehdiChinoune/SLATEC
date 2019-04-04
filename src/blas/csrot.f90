!** CSROT
SUBROUTINE CSROT(N,Cx,Incx,Cy,Incy,C,S)
  IMPLICIT NONE
  !>
  !***
  !  Apply a plane Givens rotation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B10
  !***
  ! **Type:**      COMPLEX (SROT-S, DROT-D, CSROT-C)
  !***
  ! **Keywords:**  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
  !             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
  !***
  ! **Author:**  Dongarra, J., (ANL)
  !***
  ! **Description:**
  !
  !     CSROT applies the complex Givens rotation
  !
  !          (X)   ( C S)(X)
  !          (Y) = (-S C)(Y)
  !
  !     N times where for I = 0,...,N-1
  !
  !          X = CX(LX+I*INCX)
  !          Y = CY(LY+I*INCY),
  !
  !     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !     defined in a similar way using INCY.
  !
  !     Argument Description
  !
  !        N      (integer)  number of elements in each vector
  !
  !        CX     (complex array)  beginning of one vector
  !
  !        INCX   (integer)  memory spacing of successive elements
  !               of vector CX
  !
  !        CY     (complex array)  beginning of the other vector
  !
  !        INCY   (integer)  memory spacing of successive elements
  !               of vector CY
  !
  !        C      (real)  cosine term of the rotation
  !
  !        S      (real)  sine term of the rotation.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  COMPLEX Cx(*), Cy(*), ctemp
  REAL C, S
  INTEGER i, Incx, Incy, ix, iy, N
  !* FIRST EXECUTABLE STATEMENT  CSROT
  IF ( N<=0 ) RETURN
  IF ( Incx==1.AND.Incy==1 ) THEN
    !
    !     Code for both increments equal to 1.
    !
    DO i = 1, N
      ctemp = C*Cx(i) + S*Cy(i)
      Cy(i) = C*Cy(i) - S*Cx(i)
      Cx(i) = ctemp
    END DO
    RETURN
  END IF
  !
  !     Code for unequal increments or equal increments not equal to 1.
  !
  ix = 1
  iy = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  IF ( Incy<0 ) iy = (-N+1)*Incy + 1
  DO i = 1, N
    ctemp = C*Cx(ix) + S*Cy(iy)
    Cy(iy) = C*Cy(iy) - S*Cx(ix)
    Cx(ix) = ctemp
    ix = ix + Incx
    iy = iy + Incy
  END DO
  RETURN
END SUBROUTINE CSROT
