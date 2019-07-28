!** CDSCL
PURE SUBROUTINE CDSCL(Hmax,N,Nq,Rmax,H,Rc,Rh,Yh)
  !> Subroutine CDSCL rescales the YH array whenever the stepw size is changed.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      COMPLEX (SDSCL-S, DDSCL-D, CDSCL-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  INTEGER, INTENT(IN) :: N, Nq
  REAL(SP), INTENT(IN) :: Hmax, Rmax
  REAL(SP), INTENT(INOUT) :: H, Rc, Rh
  COMPLEX(SP), INTENT(INOUT) :: Yh(N,Nq+1)
  !
  INTEGER :: i, j
  REAL(SP) :: r1
  !* FIRST EXECUTABLE STATEMENT  CDSCL
  IF( H<1._SP ) THEN
    Rh = MIN(ABS(H)*Rh,ABS(H)*Rmax,Hmax)/ABS(H)
  ELSE
    Rh = MIN(Rh,Rmax,Hmax/ABS(H))
  END IF
  r1 = 1._SP
  DO j = 1, Nq
    r1 = r1*Rh
    DO i = 1, N
      Yh(i,j+1) = Yh(i,j+1)*r1
    END DO
  END DO
  H = H*Rh
  Rc = Rc*Rh
  !
END SUBROUTINE CDSCL