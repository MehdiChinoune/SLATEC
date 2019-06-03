!** DDSCL
SUBROUTINE DDSCL(Hmax,N,Nq,Rmax,H,Rc,Rh,Yh)
  !>
  !  Subroutine DDSCL rescales the YH array whenever the step
  !            size is changed.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      DOUBLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C)
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

  INTEGER :: N, Nq
  REAL(DP) :: H, Hmax, Rc, Rh, Rmax, Yh(N,Nq+1)
  INTEGER :: i, j
  REAL(DP) :: r1
  !* FIRST EXECUTABLE STATEMENT  DDSCL
  IF ( H<1.D0 ) THEN
    Rh = MIN(ABS(H)*Rh,ABS(H)*Rmax,Hmax)/ABS(H)
  ELSE
    Rh = MIN(Rh,Rmax,Hmax/ABS(H))
  END IF
  r1 = 1.D0
  DO j = 1, Nq
    r1 = r1*Rh
    DO i = 1, N
      Yh(i,j+1) = Yh(i,j+1)*r1
    END DO
  END DO
  H = H*Rh
  Rc = Rc*Rh
END SUBROUTINE DDSCL
