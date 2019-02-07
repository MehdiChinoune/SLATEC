!*==POS3D1.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK POS3D1
SUBROUTINE POS3D1(Lp,L,Mp,M,N,A,B,C,Ldimf,Mdimf,F,Xrt,Yrt,T,D,Wx,Wy,C1,C2,&
    Bb)
  IMPLICIT NONE
  !*--POS3D16
  !*** Start of declarations inserted by SPAG
  REAL A, B, Bb, C, C1, C2, D, di, dj, dum, dx, dy, F, pi, &
    PIMACH, scalx, scaly, T, Wx, Wy
  REAL Xrt, Yrt
  INTEGER i, ifwrd, j, k, L, Ldimf, Lp, lr, lrdel, M, Mdimf, Mp, &
    mr, mrdel, N, nr
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  POS3D1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to POIS3D
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (POS3D1-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  POIS3D
  !***ROUTINES CALLED  COSQB, COSQF, COSQI, COST, COSTI, PIMACH, RFFTB,
  !                    RFFTF, RFFTI, SINQB, SINQF, SINQI, SINT, SINTI,
  !                    TRIDQ
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900308  Changed call to TRID to call to TRIDQ.  (WRB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  POS3D1
  DIMENSION A(*), B(*), C(*), F(Ldimf,Mdimf,*), Xrt(*), Yrt(*), T(*), &
    D(*), Wx(*), Wy(*), Bb(*)
  !***FIRST EXECUTABLE STATEMENT  POS3D1
  pi = PIMACH(dum)
  lr = L
  mr = M
  nr = N
  !
  !     GENERATE TRANSFORM ROOTS
  !
  lrdel = ((Lp-1)*(Lp-3)*(Lp-5))/3
  scalx = lr + lrdel
  dx = pi/(2.*scalx)
  SELECT CASE (Lp)
    CASE (1)
      Xrt(1) = 0.
      Xrt(lr) = -4.*C1
      DO i = 3, lr, 2
        Xrt(i-1) = -4.*C1*(SIN((i-1)*dx))**2
        Xrt(i) = Xrt(i-1)
      ENDDO
      CALL RFFTI(lr,Wx)
      GOTO 100
    CASE (2)
      di = 0.0
    CASE (4)
      di = 1.0
    CASE DEFAULT
      di = 0.5
      scalx = 2.*scalx
  END SELECT
  DO i = 1, lr
    Xrt(i) = -4.*C1*(SIN((i-di)*dx))**2
  ENDDO
  scalx = 2.*scalx
  SELECT CASE (Lp)
    CASE (1)
    CASE (3)
      CALL SINQI(lr,Wx)
    CASE (4)
      CALL COSTI(lr,Wx)
    CASE (5)
      CALL COSQI(lr,Wx)
    CASE DEFAULT
      CALL SINTI(lr,Wx)
  END SELECT
  100  mrdel = ((Mp-1)*(Mp-3)*(Mp-5))/3
  scaly = mr + mrdel
  dy = pi/(2.*scaly)
  SELECT CASE (Mp)
    CASE (1)
      Yrt(1) = 0.
      Yrt(mr) = -4.*C2
      DO j = 3, mr, 2
        Yrt(j-1) = -4.*C2*(SIN((j-1)*dy))**2
        Yrt(j) = Yrt(j-1)
      ENDDO
      CALL RFFTI(mr,Wy)
      GOTO 200
    CASE (2)
      dj = 0.0
    CASE (4)
      dj = 1.0
    CASE DEFAULT
      dj = 0.5
      scaly = 2.*scaly
  END SELECT
  DO j = 1, mr
    Yrt(j) = -4.*C2*(SIN((j-dj)*dy))**2
  ENDDO
  scaly = 2.*scaly
  SELECT CASE (Mp)
    CASE (1)
    CASE (3)
      CALL SINQI(mr,Wy)
    CASE (4)
      CALL COSTI(mr,Wy)
    CASE (5)
      CALL COSQI(mr,Wy)
    CASE DEFAULT
      CALL SINTI(mr,Wy)
  END SELECT
  200 CONTINUE
  IFwrd = 1
  DO
    !
    !     TRANSFORM X
    !
    DO j = 1, mr
      DO k = 1, nr
        DO i = 1, lr
          T(i) = F(i,j,k)
        ENDDO
        SELECT CASE (Lp)
          CASE (2)
            CALL SINT(lr,T,Wx)
          CASE (3)
            IF ( ifwrd==2 ) THEN
              CALL SINQB(lr,T,Wx)
            ELSE
              CALL SINQF(lr,T,Wx)
            ENDIF
          CASE (4)
            CALL COST(lr,T,Wx)
          CASE (5)
            IF ( ifwrd==2 ) THEN
              CALL COSQB(lr,T,Wx)
            ELSE
              CALL COSQF(lr,T,Wx)
            ENDIF
          CASE DEFAULT
            IF ( ifwrd==2 ) THEN
              CALL RFFTB(lr,T,Wx)
            ELSE
              CALL RFFTF(lr,T,Wx)
            ENDIF
        END SELECT
        DO i = 1, lr
          F(i,j,k) = T(i)
        ENDDO
      ENDDO
    ENDDO
    IF ( ifwrd==2 ) THEN
      DO i = 1, lr
        DO j = 1, mr
          DO k = 1, nr
            F(i,j,k) = F(i,j,k)/(scalx*scaly)
          ENDDO
        ENDDO
      ENDDO
      EXIT
    ELSE
      DO
        !
        !     TRANSFORM Y
        !
        DO i = 1, lr
          DO k = 1, nr
            DO j = 1, mr
              T(j) = F(i,j,k)
            ENDDO
            SELECT CASE (Mp)
              CASE (2)
                CALL SINT(mr,T,Wy)
              CASE (3)
                IF ( ifwrd==2 ) THEN
                  CALL SINQB(mr,T,Wy)
                ELSE
                  CALL SINQF(mr,T,Wy)
                ENDIF
              CASE (4)
                CALL COST(mr,T,Wy)
              CASE (5)
                IF ( ifwrd==2 ) THEN
                  CALL COSQB(mr,T,Wy)
                ELSE
                  CALL COSQF(mr,T,Wy)
                ENDIF
              CASE DEFAULT
                IF ( ifwrd==2 ) THEN
                  CALL RFFTB(mr,T,Wy)
                ELSE
                  CALL RFFTF(mr,T,Wy)
                ENDIF
            END SELECT
            DO j = 1, mr
              F(i,j,k) = T(j)
            ENDDO
          ENDDO
        ENDDO
        IF ( ifwrd==2 ) EXIT
        !
        !     SOLVE TRIDIAGONAL SYSTEMS IN Z
        !
        DO i = 1, lr
          DO j = 1, mr
            DO k = 1, nr
              Bb(k) = B(k) + Xrt(i) + Yrt(j)
              T(k) = F(i,j,k)
            ENDDO
            CALL TRIDQ(nr,A,Bb,C,T,D)
            DO k = 1, nr
              F(i,j,k) = T(k)
            ENDDO
          ENDDO
        ENDDO
        ifwrd = 2
      ENDDO
    ENDIF
  ENDDO
END SUBROUTINE POS3D1
