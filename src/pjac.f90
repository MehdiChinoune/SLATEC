!*==PJAC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PJAC
      SUBROUTINE PJAC(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wm,Iwm,F,JAC,Rpar,Ipar)
      IMPLICIT NONE
!*--PJAC5
!*** Start of declarations inserted by SPAG
      INTEGER Ipar
      REAL Rpar
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  PJAC
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DEBDF
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PJAC-S, DPJAC-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   PJAC sets up the iteration matrix (involving the Jacobian) for the
!   integration package DEBDF.
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  SGBFA, SGEFA, VNWRMS
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  PJAC
!
!LLL. OPTIMIZE
      INTEGER Neq , Nyh , Iwm , i , i1 , i2 , IER , ii , IOWnd , IOWns , j , 
     &        j1 , jj , JSTart , KFLag , L , lenp , MAXord , mba , mband , 
     &        meb1 , meband , METh , MITer , ml , ml3 , mu , N , NFE , NJE , 
     &        NQ , NQU , NST
      EXTERNAL F , JAC
      REAL Y , Yh , Ewt , Ftem , Savf , Wm , ROWnd , ROWns , EL0 , H , HMIn , 
     &     HMXi , HU , TN , UROund , con , di , fac , hl0 , r , r0 , srur , yi , 
     &     yj , yjj , VNWRMS
      DIMENSION Y(*) , Yh(Nyh,*) , Ewt(*) , Ftem(*) , Savf(*) , Wm(*) , Iwm(*) , 
     &          Rpar(*) , Ipar(*)
      COMMON /DEBDF1/ ROWnd , ROWns(210) , EL0 , H , HMIn , HMXi , HU , TN , 
     &                UROund , IOWnd(14) , IOWns(6) , IER , JSTart , KFLag , L , 
     &                METh , MITer , MAXord , N , NQ , NST , NFE , NJE , NQU
!-----------------------------------------------------------------------
! PJAC IS CALLED BY STOD  TO COMPUTE AND PROCESS THE MATRIX
! P = I - H*EL(1)*J , WHERE J IS AN APPROXIMATION TO THE JACOBIAN.
! HERE J IS COMPUTED BY THE USER-SUPPLIED ROUTINE JAC IF
! MITER = 1 OR 4, OR BY FINITE DIFFERENCING IF MITER = 2, 3, OR 5.
! IF MITER = 3, A DIAGONAL APPROXIMATION TO J IS USED.
! J IS STORED IN WM AND REPLACED BY P.  IF MITER .NE. 3, P IS THEN
! SUBJECTED TO LU DECOMPOSITION IN PREPARATION FOR LATER SOLUTION
! OF LINEAR SYSTEMS WITH P AS COEFFICIENT MATRIX. THIS IS DONE
! BY SGEFA IF MITER = 1 OR 2, AND BY SGBFA IF MITER = 4 OR 5.
!
! IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
! WITH PJAC USES THE FOLLOWING..
! Y    = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
! FTEM = WORK ARRAY OF LENGTH N (ACOR IN STOD ).
! SAVF = ARRAY CONTAINING F EVALUATED AT PREDICTED Y.
! WM   = REAL WORK SPACE FOR MATRICES.  ON OUTPUT IT CONTAINS THE
!        INVERSE DIAGONAL MATRIX IF MITER = 3 AND THE LU DECOMPOSITION
!        OF P IF MITER IS 1, 2 , 4, OR 5.
!        STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
!        WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
!        WM(1) = SQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.
!        WM(2) = H*EL0, SAVED FOR LATER USE IF MITER = 3.
! IWM  = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING AT
!        IWM(21), IF MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS THE
!        BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) IF MITER IS 4 OR 5.
! EL0  = EL(1) (INPUT).
! IER  = OUTPUT ERROR FLAG,  = 0 IF NO TROUBLE, .NE. 0 IF
!        P MATRIX FOUND TO BE SINGULAR.
! THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, TN, UROUND,
! MITER, N, NFE, AND NJE.
!-----------------------------------------------------------------------
!***FIRST EXECUTABLE STATEMENT  PJAC
      NJE = NJE + 1
      hl0 = H*EL0
      SELECT CASE (MITer)
      CASE (2)
! IF MITER = 2, MAKE N CALLS TO F TO APPROXIMATE J. --------------------
        fac = VNWRMS(N,Savf,Ewt)
        r0 = 1000.0E0*ABS(H)*UROund*N*fac
        IF ( r0==0.0E0 ) r0 = 1.0E0
        srur = Wm(1)
        j1 = 2
        DO j = 1 , N
          yj = Y(j)
          r = MAX(srur*ABS(yj),r0*Ewt(j))
          Y(j) = Y(j) + r
          fac = -hl0/r
          CALL F(TN,Y,Ftem,Rpar,Ipar)
          DO i = 1 , N
            Wm(i+j1) = (Ftem(i)-Savf(i))*fac
          ENDDO
          Y(j) = yj
          j1 = j1 + N
        ENDDO
        NFE = NFE + N
      CASE (3)
! IF MITER = 3, CONSTRUCT A DIAGONAL APPROXIMATION TO J AND P. ---------
        Wm(2) = hl0
        IER = 0
        r = EL0*0.1E0
        DO i = 1 , N
          Y(i) = Y(i) + r*(H*Savf(i)-Yh(i,2))
        ENDDO
        CALL F(TN,Y,Wm(3),Rpar,Ipar)
        NFE = NFE + 1
        DO i = 1 , N
          r0 = H*Savf(i) - Yh(i,2)
          di = 0.1E0*r0 - H*(Wm(i+2)-Savf(i))
          Wm(i+2) = 1.0E0
          IF ( ABS(r0)>=UROund*Ewt(i) ) THEN
            IF ( ABS(di)==0.0E0 ) GOTO 100
            Wm(i+2) = 0.1E0*r0/di
          ENDIF
        ENDDO
        RETURN
      CASE (4)
! IF MITER = 4, CALL JAC AND MULTIPLY BY SCALAR. -----------------------
        ml = Iwm(1)
        mu = Iwm(2)
        ml3 = 3
        mband = ml + mu + 1
        meband = mband + ml
        lenp = meband*N
        DO i = 1 , lenp
          Wm(i+2) = 0.0E0
        ENDDO
        CALL JAC(TN,Y,Wm(ml3),meband,Rpar,Ipar)
        con = -hl0
        DO i = 1 , lenp
          Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 200
      CASE (5)
! IF MITER = 5, MAKE MBAND CALLS TO F TO APPROXIMATE J. ----------------
        ml = Iwm(1)
        mu = Iwm(2)
        mband = ml + mu + 1
        mba = MIN(mband,N)
        meband = mband + ml
        meb1 = meband - 1
        srur = Wm(1)
        fac = VNWRMS(N,Savf,Ewt)
        r0 = 1000.0E0*ABS(H)*UROund*N*fac
        IF ( r0==0.0E0 ) r0 = 1.0E0
        DO j = 1 , mba
          DO i = j , N , mband
            yi = Y(i)
            r = MAX(srur*ABS(yi),r0*Ewt(i))
            Y(i) = Y(i) + r
          ENDDO
          CALL F(TN,Y,Ftem,Rpar,Ipar)
          DO jj = j , N , mband
            Y(jj) = Yh(jj,1)
            yjj = Y(jj)
            r = MAX(srur*ABS(yjj),r0*Ewt(jj))
            fac = -hl0/r
            i1 = MAX(jj-mu,1)
            i2 = MIN(jj+ml,N)
            ii = jj*meb1 - ml + 2
            DO i = i1 , i2
              Wm(ii+i) = (Ftem(i)-Savf(i))*fac
            ENDDO
          ENDDO
        ENDDO
        NFE = NFE + mba
        GOTO 200
      CASE DEFAULT
! IF MITER = 1, CALL JAC AND MULTIPLY BY SCALAR. -----------------------
        lenp = N*N
        DO i = 1 , lenp
          Wm(i+2) = 0.0E0
        ENDDO
        CALL JAC(TN,Y,Wm(3),N,Rpar,Ipar)
        con = -hl0
        DO i = 1 , lenp
          Wm(i+2) = Wm(i+2)*con
        ENDDO
      END SELECT
! ADD IDENTITY MATRIX. -------------------------------------------------
      j = 3
      DO i = 1 , N
        Wm(j) = Wm(j) + 1.0E0
        j = j + (N+1)
      ENDDO
! DO LU DECOMPOSITION ON P. --------------------------------------------
      CALL SGEFA(Wm(3),N,N,Iwm(21),IER)
      RETURN
 100  IER = -1
      RETURN
! ADD IDENTITY MATRIX. -------------------------------------------------
 200  ii = mband + 2
      DO i = 1 , N
        Wm(ii) = Wm(ii) + 1.0E0
        ii = ii + meband
      ENDDO
! DO LU DECOMPOSITION OF P. --------------------------------------------
      CALL SGBFA(Wm(3),meband,N,ml,mu,Iwm(21),IER)
!----------------------- END OF SUBROUTINE PJAC -----------------------
      END SUBROUTINE PJAC
