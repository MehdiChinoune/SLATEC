!** DCSCAL
SUBROUTINE DCSCAL(A,Nrda,Nrow,Ncol,Cols,Colsav,Rows,Rowsav,Anorm,Scales,Iscale,Ic)
  !> Subsidiary to DBVSUP and DSUDS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (CSCALE-S, DCSCAL-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !     This routine scales the matrix A by columns when needed.
  !
  !***
  ! **See also:**  DBVSUP, DSUDS
  !***
  ! **Routines called:**  DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  INTEGER :: Ic, Iscale, Ncol, Nrda, Nrow
  REAL(DP) :: Anorm, A(Nrda,Ncol), Cols(Ncol), Colsav(Ncol), Rows(Nrow), Rowsav(Nrow), &
    Scales(Ncol)
  INTEGER :: ip, j, k
  REAL(DP) :: alog2, ascale, cs, p, s
  REAL(DP), PARAMETER :: ten4 = 1.E4_DP, ten20 = 1.E20_DP
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 130
  !        BEGIN BLOCK PERMITTING ...EXITS TO 60
  !* FIRST EXECUTABLE STATEMENT  DCSCAL
  IF( Iscale==(-1) ) THEN
    !
    IF( Ic/=0 ) THEN
      DO k = 1, Ncol
        Cols(k) = NORM2(A(1:Nrow,k))**2
      END DO
    END IF
    !
    ascale = Anorm/Ncol
    DO k = 1, Ncol
      cs = Cols(k)
      !        .........EXIT
      IF( (cs>ten4*ascale) .OR. (ten4*cs<ascale) ) GOTO 100
      !        .........EXIT
      IF( (cs<1._DP/ten20) .OR. (cs>ten20) ) GOTO 100
    END DO
  END IF
  !
  DO k = 1, Ncol
    Scales(k) = 1._DP
  END DO
  !     ......EXIT
  RETURN
  !
  100  alog2 = LOG(2._DP)
  Anorm = 0._DP
  DO k = 1, Ncol
    cs = Cols(k)
    IF( cs/=0._DP ) THEN
      p = LOG(cs)/alog2
      ip = INT( -0.5_DP*p )
      s = 2._DP**ip
      Scales(k) = s
      IF( Ic/=1 ) THEN
        Cols(k) = s*s*Cols(k)
        Anorm = Anorm + Cols(k)
        Colsav(k) = Cols(k)
      END IF
      DO j = 1, Nrow
        A(j,k) = s*A(j,k)
      END DO
    ELSE
      Scales(k) = 1._DP
    END IF
  END DO
  !
  !     ...EXIT
  IF( Ic/=0 ) THEN
    !
    DO k = 1, Nrow
      Rows(k) = NORM2(A(k,1:Ncol))**2
      Rowsav(k) = Rows(k)
      Anorm = Anorm + Rows(k)
    END DO
  END IF
  RETURN
END SUBROUTINE DCSCAL
