!*==DCSCAL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DCSCAL
SUBROUTINE DCSCAL(A,Nrda,Nrow,Ncol,Cols,Colsav,Rows,Rowsav,Anorm,Scales,&
    Iscale,Ic)
  IMPLICIT NONE
  !*--DCSCAL6
  !***BEGIN PROLOGUE  DCSCAL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP and DSUDS
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CSCALE-S, DCSCAL-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     This routine scales the matrix A by columns when needed.
  !
  !***SEE ALSO  DBVSUP, DSUDS
  !***ROUTINES CALLED  DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DCSCAL
  REAL(8) :: DDOT
  INTEGER Ic , ip , Iscale , j , k , Ncol , Nrda , Nrow
  REAL(8) :: A(Nrda,*) , alog2 , Anorm , ascale , Cols(*) , Colsav(*)&
    , cs , p , Rows(*) , Rowsav(*) , s , Scales(*) , ten20 , &
    ten4
  !
  SAVE ten4 , ten20
  DATA ten4 , ten20/1.0D4 , 1.0D20/
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 130
  !        BEGIN BLOCK PERMITTING ...EXITS TO 60
  !***FIRST EXECUTABLE STATEMENT  DCSCAL
  IF ( Iscale==(-1) ) THEN
    !
    IF ( Ic/=0 ) THEN
      DO k = 1 , Ncol
        Cols(k) = DDOT(Nrow,A(1,k),1,A(1,k),1)
      ENDDO
    ENDIF
    !
    ascale = Anorm/Ncol
    DO k = 1 , Ncol
      cs = Cols(k)
      !        .........EXIT
      IF ( (cs>ten4*ascale).OR.(ten4*cs<ascale) ) GOTO 100
      !        .........EXIT
      IF ( (cs<1.0D0/ten20).OR.(cs>ten20) ) GOTO 100
    ENDDO
  ENDIF
  !
  DO k = 1 , Ncol
    Scales(k) = 1.0D0
  ENDDO
  !     ......EXIT
  GOTO 99999
  !
  100  alog2 = LOG(2.0D0)
  Anorm = 0.0D0
  DO k = 1 , Ncol
    cs = Cols(k)
    IF ( cs/=0.0D0 ) THEN
      p = LOG(cs)/alog2
      ip = -0.5D0*p
      s = 2.0D0**ip
      Scales(k) = s
      IF ( Ic/=1 ) THEN
        Cols(k) = s*s*Cols(k)
        Anorm = Anorm + Cols(k)
        Colsav(k) = Cols(k)
      ENDIF
      DO j = 1 , Nrow
        A(j,k) = s*A(j,k)
      ENDDO
    ELSE
      Scales(k) = 1.0D0
    ENDIF
  ENDDO
  !
  !     ...EXIT
  IF ( Ic/=0 ) THEN
    !
    DO k = 1 , Nrow
      Rows(k) = DDOT(Ncol,A(k,1),Nrda,A(k,1),Nrda)
      Rowsav(k) = Rows(k)
      Anorm = Anorm + Rows(k)
    ENDDO
  ENDIF
  99999 END SUBROUTINE DCSCAL
