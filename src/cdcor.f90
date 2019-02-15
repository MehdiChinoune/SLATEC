!DECK CDCOR
SUBROUTINE CDCOR(Dfdy,El,FA,H,Ierror,Impl,Ipvt,Matdim,Miter,Ml,Mu,N,Nde,&
    Nq,T,USERS,Y,Yh,Ywt,Evalfa,Save1,Save2,A,D,Jstate)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CDCOR
  !***SUBSIDIARY
  !***PURPOSE  Subroutine CDCOR computes corrections to the Y array.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      COMPLEX (SDCOR-S, DDCOR-D, CDCOR-C)
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***DESCRIPTION
  !
  !  In the case of functional iteration, update Y directly from the
  !  result of the last call to F.
  !  In the case of the chord method, compute the corrector error and
  !  solve the linear system with that as right hand side and DFDY as
  !  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
  !  or 5.
  !
  !***ROUTINES CALLED  CGBSL, CGESL, SCNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  CDCOR
  INTEGER i, Ierror, iflag, Impl, j, Jstate, Matdim, Miter, Ml, &
    Mu, mw, N, Nde, Nq
  COMPLEX A(Matdim,*), Dfdy(Matdim,*), Save1(*), Save2(*), Y(*), &
    Yh(N,*), Ywt(*)
  REAL D, El(13,12), H, SCNRM2, T
  INTEGER Ipvt(*)
  LOGICAL Evalfa
  !***FIRST EXECUTABLE STATEMENT  CDCOR
  IF ( Miter==0 ) THEN
    IF ( Ierror==1.OR.Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = (H*Save2(i)-Yh(i,2)-Save1(i))/Ywt(i)
      ENDDO
    ELSE
      DO i = 1, N
        Save1(i) = (H*Save2(i)-Yh(i,2)-Save1(i))/MAX(ABS(Y(i)),ABS(Ywt(i)))
      ENDDO
    ENDIF
    D = SCNRM2(N,Save1,1)/SQRT(REAL(N))
    DO i = 1, N
      Save1(i) = H*Save2(i) - Yh(i,2)
    ENDDO
  ELSEIF ( Miter==1.OR.Miter==2 ) THEN
    IF ( Impl==0 ) THEN
      DO i = 1, N
        Save2(i) = H*Save2(i) - Yh(i,2) - Save1(i)
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      ENDDO
      DO j = 1, N
        DO i = 1, N
          Save2(i) = Save2(i) - A(i,j)*(Yh(j,2)+Save1(j))
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i) - A(i,1)*(Yh(i,2)+Save1(i))
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      ENDDO
      DO j = 1, Nde
        DO i = 1, Nde
          Save2(i) = Save2(i) - A(i,j)*(Yh(j,2)+Save1(j))
        ENDDO
      ENDDO
    ENDIF
    CALL CGESL(Dfdy,Matdim,N,Ipvt,Save2,0)
    IF ( Ierror==1.OR.Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      ENDDO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      ENDDO
    ENDIF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N))
  ELSEIF ( Miter==4.OR.Miter==5 ) THEN
    IF ( Impl==0 ) THEN
      DO i = 1, N
        Save2(i) = H*Save2(i) - Yh(i,2) - Save1(i)
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      ENDDO
      mw = Ml + 1 + Mu
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Save2(i+j-mw) = Save2(i+j-mw) - A(i,j)*(Yh(j,2)+Save1(j))
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i) - A(i,1)*(Yh(i,2)+Save1(i))
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      IF ( Evalfa ) THEN
        CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
        IF ( N==0 ) THEN
          Jstate = 9
          RETURN
        ENDIF
      ELSE
        Evalfa = .TRUE.
      ENDIF
      DO i = 1, N
        Save2(i) = H*Save2(i)
      ENDDO
      mw = Ml + 1 + Mu
      DO j = 1, Nde
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+Nde-j,mw+Ml)
          Save2(i+j-mw) = Save2(i+j-mw) - A(i,j)*(Yh(j,2)+Save1(j))
        ENDDO
      ENDDO
    ENDIF
    CALL CGBSL(Dfdy,Matdim,N,Ml,Mu,Ipvt,Save2,0)
    IF ( Ierror==1.OR.Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      ENDDO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      ENDDO
    ENDIF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N))
  ELSEIF ( Miter==3 ) THEN
    iflag = 2
    CALL USERS(Y,Yh(1,2),Ywt,Save1,Save2,T,H,El(1,Nq),Impl,N,Nde,iflag)
    IF ( N==0 ) THEN
      Jstate = 10
      RETURN
    ENDIF
    IF ( Ierror==1.OR.Ierror==5 ) THEN
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/Ywt(i)
      ENDDO
    ELSE
      DO i = 1, N
        Save1(i) = Save1(i) + Save2(i)
        Save2(i) = Save2(i)/MAX(ABS(Y(i)),ABS(Ywt(i)))
      ENDDO
    ENDIF
    D = SCNRM2(N,Save2,1)/SQRT(REAL(N))
  ENDIF
END SUBROUTINE CDCOR
