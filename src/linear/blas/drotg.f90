!** DROTG
SUBROUTINE DROTG(Da,Db,Dc,Ds)
  !>
  !***
  !  Construct a plane Givens rotation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B10
  !***
  ! **Type:**      DOUBLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
  !***
  ! **Keywords:**  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
  !             LINEAR ALGEBRA, VECTOR
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
  !       DA  double precision scalar
  !       DB  double precision scalar
  !
  !     --Output--
  !       DA  double precision result R
  !       DB  double precision result Z
  !       DC  double precision result
  !       DS  double precision result
  !
  !     Construct the Givens transformation
  !
  !         ( DC  DS )
  !     G = (        ),    DC**2 + DS**2 = 1 ,
  !         (-DS  DC )
  !
  !     which zeros the second entry of the 2-vector  (DA,DB)**T .
  !
  !     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
  !     storage.  The value of DB is overwritten by a value Z which
  !     allows DC and DS to be recovered by the following algorithm.
  !
  !           If Z=1  set  DC=0.0  and  DS=1.0
  !           If ABS(Z) .LT. 1  set  DC=SQRT(1-Z**2)  and  DS=Z
  !           If ABS(Z) .GT. 1  set  DC=1/Z  and  DS=SQRT(1-DC**2)
  !
  !     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will
  !     next be called to apply the transformation to a 2 by N matrix.
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL(8) :: Da, Db, Dc, Ds, u, v, r
  !* FIRST EXECUTABLE STATEMENT  DROTG
  IF ( ABS(Da)>ABS(Db) ) THEN
    !
    !- ** HERE ABS(DA) .GT. ABS(DB) ***
    !
    u = Da + Da
    v = Db/u
    !
    !     NOTE THAT U AND R HAVE THE SIGN OF DA
    !
    r = SQRT(0.25D0+v**2)*u
    !
    !     NOTE THAT DC IS POSITIVE
    !
    Dc = Da/r
    Ds = v*(Dc+Dc)
    Db = Ds
    Da = r
    RETURN
    !
    !- ** HERE ABS(DA) .LE. ABS(DB) ***
    !
  ELSEIF ( Db==0.0D0 ) THEN
    !
    !- ** HERE DA = DB = 0.0 ***
    !
    Dc = 1.0D0
    Ds = 0.0D0
    RETURN
  ELSE
    u = Db + Db
    v = Da/u
    !
    !     NOTE THAT U AND R HAVE THE SIGN OF DB
    !     (R IS IMMEDIATELY STORED IN DA)
    !
    Da = SQRT(0.25D0+v**2)*u
    !
    !     NOTE THAT DS IS POSITIVE
    !
    Ds = Db/Da
    Dc = v*(Ds+Ds)
    IF ( Dc/=0.0D0 ) THEN
      Db = 1.0D0/Dc
      RETURN
    END IF
  END IF
  Db = 1.0D0
  RETURN
END SUBROUTINE DROTG
