!** CDCOR
PURE SUBROUTINE CDCOR(Dfdy,El,FA,H,Ierror,Impl,Ipvt,Matdim,Miter,Ml,Mu,N,Nde,&
    Nq,T,USERS,Y,Yh,Ywt,Evalfa,Save1,Save2,A,D,Jstate)
  !> Subroutine CDCOR computes corrections to the Y array.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      COMPLEX (SDCOR-S, DDCOR-D, CDCOR-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
  !
  !  In the case of functional iteration, update Y directly from the
  !  result of the last call to F.
  !  In the case of the chord method, compute the corrector error and
  !  solve the linear system with that as right hand side and DFDY as
  !  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4, or 5.
  !
  !***
  ! **Routines called:**  CGBSL, CGESL, SCNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  USE blas, ONLY : SCNRM2
  USE lapack, ONLY : CGBTRS, CGETRS
  !
  INTERFACE
    PURE SUBROUTINE USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      IMPORT SP
      INTEGER, INTENT(IN) :: Impl, N, Nde, iflag
      REAL(SP), INTENT(IN) :: T, H, El
      COMPLEX(SP), INTENT(IN) :: Y(N), Yh(N,13), Ywt(N)
      COMPLEX(SP), INTENT(INOUT) :: Save1(N), Save2(N)
    END SUBROUTINE USERS
    PURE SUBROUTINE FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IMPORT SP
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu, Nde
      REAL(SP), INTENT(IN) :: T
      COMPLEX(SP), INTENT(IN) :: Y(N)
      COMPLEX(SP), INTENT(INOUT) :: A(:,:)
    END SUBROUTINE FA
  END INTERFACE
  INTEGER, INTENT(IN) :: Ierror, Impl, Matdim, Miter, Ml, Mu, N, Nde, Nq
  INTEGER, INTENT(OUT) :: Jstate
  INTEGER, INTENT(IN) :: Ipvt(N)
  REAL(SP), INTENT(IN) :: H, T, El(13,12)
  REAL(SP), INTENT(OUT) :: D
  COMPLEX(SP), INTENT(IN) :: Dfdy(Matdim,N), Y(N), Yh(N,13), Ywt(N)
  COMPLEX(SP), INTENT(INOUT) :: A(Matdim,N), Save1(N), Save2(N)
  LOGICAL, INTENT(INOUT) :: Evalfa
  !
  INTEGER :: i, iflag, j, mw, info
  !* FIRST EXECUTABLE STATEMENT  CDCOR
  IF( Miter==0 ) THEN
    IF( Ierror==1 .OR. Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = (H*Save2(i)-Yh(i,2)-Save1(i))/Ywt(i)
      END DO
    ELSE
      DO i = 1, N
        Save1(i) = (H*Save2(i)-Yh(i,2)-Save1(i))/MAX(ABS(Y(i)),ABS(Ywt(i)))
      END DO
    END IF
    D = SCNRM2(N,Save1,1)/SQRT(REAL(N,SP))
    DO i = 1, N
      Save1(i) = H*Save2(i) - Yh(i,2)
    END DO
  ELSEIF( Miter==1 .OR. Miter==2 ) THEN
    IF( Impl==0 ) THEN
      DO i = 1, N
        Save2(i) = H*Save2(i) - Yh(i,2) - Save1(i)
      END DO
    ELSEIF( Impl==1 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      END DO
      DO j = 1, N
        DO i = 1, N
          Save2(i) = Save2(i) - A(i,j)*(Yh(j,2)+Save1(j))
        END DO
      END DO
    ELSEIF( Impl==2 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i) - A(i,1)*(Yh(i,2)+Save1(i))
      END DO
    ELSEIF( Impl==3 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      END DO
      DO j = 1, Nde
        DO i = 1, Nde
          Save2(i) = Save2(i) - A(i,j)*(Yh(j,2)+Save1(j))
        END DO
      END DO
    END IF
    CALL CGETRS('N',N,1,Dfdy,Matdim,Ipvt,Save2,N,info)
    IF( Ierror==1 .OR. Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      END DO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      END DO
    END IF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N,SP))
  ELSEIF( Miter==4 .OR. Miter==5 ) THEN
    IF( Impl==0 ) THEN
      DO i = 1, N
        Save2(i) = H*Save2(i) - Yh(i,2) - Save1(i)
      END DO
    ELSEIF( Impl==1 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      END DO
      mw = Ml + 1 + Mu
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Save2(i+j-mw) = Save2(i+j-mw) - A(i,j)*(Yh(j,2)+Save1(j))
        END DO
      END DO
    ELSEIF( Impl==2 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i) - A(i,1)*(Yh(i,2)+Save1(i))
      END DO
    ELSEIF( Impl==3 ) THEN
      IF( Evalfa ) THEN
        CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
      ELSE
        Evalfa = .TRUE.
      END IF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      END DO
      mw = Ml + 1 + Mu
      DO j = 1, Nde
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+Nde-j,mw+Ml)
          Save2(i+j-mw) = Save2(i+j-mw) - A(i,j)*(Yh(j,2)+Save1(j))
        END DO
      END DO
    END IF
    CALL CGBTRS('N',N,Ml,Mu,1,Dfdy,Matdim,Ipvt,Save2,N,info)
    IF( Ierror==1 .OR. Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      END DO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      END DO
    END IF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N,SP))
  ELSEIF( Miter==3 ) THEN
    iflag = 2
    CALL USERS(Y,Yh(1,2),Ywt,Save1,Save2,T,H,El(1,Nq),Impl,N,Nde,iflag)
    IF( N==0 ) THEN
      Jstate = 10
      RETURN
    END IF
    IF( Ierror==1 .OR. Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      END DO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      END DO
    END IF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N,SP))
  END IF
  !
END SUBROUTINE CDCOR