!** BVPOR
SUBROUTINE BVPOR(Y,Nrowy,Ncomp,Xpts,Nxpts,A,Nrowa,Alpha,Nic,B,Nrowb,Beta,&
    Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,W,Niv,Yhp,U,V,Coef,S,Stowa,G,Work,Iwork,Nfcc)
  !>
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BVPOR-S, DBVPOR-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !     INPUT to BVPOR    (items not defined in BVSUP comments)
  !- *********************************************************************
  !
  !     NOPG = 0 -- Orthonormalization points not pre-assigned
  !          = 1 -- Orthonormalization points pre-assigned
  !
  !     MXNON = Maximum number of orthogonalizations allowed.
  !
  !     NDISK = 0 -- IN-CORE storage
  !           = 1 -- DISK storage.  Value of NTAPE in data statement
  !                  is set to 13.  If another value is desired,
  !                  the data statement must be changed.
  !
  !     INTEG = Type of integrator and associated test to be used
  !             to determine when to orthonormalize.
  !
  !             1 -- Use GRAM-SCHMIDT test and DERKF
  !             2 -- Use GRAM-SCHMIDT test and DEABM
  !
  !     TOL = Tolerance for allowable error in orthogonalization test.
  !
  !     NPS = 0 Normalize particular solution to unit length at each
  !             point of orthonormalization.
  !         = 1 Do not normalize particular solution.
  !
  !     NTP = Must be .GE. NFC*(NFC+1)/2.
  !
  !
  !     NFCC = 2*NFC for special treatment of a complex valued problem
  !
  !     ICOCO = 0 Skip final computations (superposition coefficients
  !               and ,hence, boundary problem solution)
  !           = 1 Calculate superposition coefficients and obtain
  !               solution to the boundary value problem
  !
  !- *********************************************************************
  !     OUTPUT from BVPOR
  !- *********************************************************************
  !
  !     Y(NROWY,NXPTS) = Solution at specified output points.
  !
  !     MXNON = Number of orthonormalizations performed by BVPOR.
  !
  !     Z(MXNON+1) = Locations of orthonormalizations performed by BVPOR.
  !
  !     NIV = Number of independent vectors returned from MGSBV. Normally
  !        this parameter will be meaningful only when MGSBV returns with
  !           MFLAG = 2.
  !
  !- *********************************************************************
  !
  !     The following variables are in the argument list because of
  !     variable dimensioning. In general, they contain no information of
  !     use to the user.  The amount of storage set aside by the user must
  !     be greater than or equal to that indicated by the dimension
  !     statements.   For the DISK storage mode, NON = 0 and KPTS = 1,
  !     while for the IN-CORE storage mode, NON = MXNON and KPTS = NXPTS.
  !
  !     P(NTP,NON+1)
  !     IP(NFCC,NON+1)
  !     YHP(NCOMP,NFC+1)  plus an additional column of the length  NEQIVP
  !     U(NCOMP,NFC,KPTS)
  !     V(NCOMP,KPTS)
  !     W(NFCC,NON+1)
  !     COEF(NFCC)
  !     S(NFC+1)
  !     STOWA(NCOMP*(NFC+1)+NEQIVP+1)
  !     G(NCOMP)
  !     WORK(KKKWS)
  !     IWORK(LLLIWS)
  !
  !- *********************************************************************
  !     Subroutines used by BVPOR
  !         LSSUDS -- Solves an underdetermined system of linear
  !                   equations.  This routine is used to get a full
  !                   set of initial conditions for integration.
  !                   Called by BVPOR
  !
  !         SVECS -- Obtains starting vectors for special treatment
  !                  of complex valued problems, called by BVPOR
  !
  !         RKFAB -- Routine which conducts integration using DERKF or
  !                   DEABM
  !
  !         STWAY -- Storage for backup capability, called by
  !                   BVPOR and REORT
  !
  !         STOR1 -- Storage at output points, called by BVPOR,
  !                  RKFAB, REORT and STWAY.
  !
  !         SDOT -- Single precision vector inner product routine,
  !                   called by BVPOR, SCOEF, LSSUDS, MGSBV,
  !                   BKSOL, REORT and PRVEC.
  !         ** NOTE **
  !         A considerable improvement in speed can be achieved if a
  !         machine language version is used for SDOT.
  !
  !         SCOEF -- Computes the superposition constants from the
  !                  boundary conditions at Xfinal.
  !
  !         BKSOL -- Solves an upper triangular set of linear equations.
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  BKSOL, LSSUDS, RKFAB, SCOEF, SDOT, STOR1, STWAY,
  !                    SVECS
  !***
  ! COMMON BLOCKS    ML15TO, ML18JR, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY : C, INHomo, IVP, PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, KNSwot, &
    KOP, LOTjp, NSWot, AE, RE, TOL, NDIsk, NTApe, NEQ, NEQivp, NUMort, ICOco
  INTEGER :: Mxnon, Ncomp, Nfc, Nfcc, Nic, Niv, Nrowa, Nrowb, Nrowy, Ntp, Iflag, &
    Nxpts
  INTEGER :: Ip(Nfcc,Mxnon+1), Iwork(*)
  REAL :: A(Nrowa,Ncomp), Alpha(:), B(Nrowb,Ncomp), Beta(Nfc), Coef(Nfcc), &
    G(Ncomp), P(Ntp,Mxnon+1), S(Nfc+1), Stowa(:), U(Ncomp,Nfc,Nxpts), &
    V(Ncomp,Nxpts), W(Nfcc,Mxnon+1), Work(*), Xpts(Nxpts), Y(Nrowy,Nxpts), &
    Yhp(Ncomp,Nfc+1), Z(Mxnon+1)
  INTEGER kod, kpts, kwc, kwd, kws, kwt, l, m, n, ncomp2, ndw, &
    nfcp1, nfcp2, nn, non, i, i1, i2, ic, ira, isflg, j, k
  !* FIRST EXECUTABLE STATEMENT  BVPOR
  nfcp1 = Nfc + 1
  NUMort = 0
  C = 1.0
  !
  !- *********************************************************************
  !     CALCULATE INITIAL CONDITIONS WHICH SATISFY
  !                   A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
  !     WHEN NFC .NE. NFCC LSSUDS DEFINES VALUES YHP IN A MATRIX OF SIZE
  !     (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE ALLOCATION INTO
  !     THE U ARRAY. HOWEVER, THIS IS OKAY SINCE PLENTY OF SPACE IS
  !     AVAILABLE IN U AND IT HAS NOT YET BEEN USED.
  !
  ndw = Nrowa*Ncomp
  kws = ndw + Nic + 1
  kwd = kws + Nic
  kwt = kwd + Nic
  kwc = kwt + Nic
  Iflag = 0
  CALL LSSUDS(A,Yhp(1,Nfcc+1),Alpha,Nic,Ncomp,Nrowa,Yhp,Ncomp,Iflag,1,ira,0,&
    Work(1),Work(ndw+1),Iwork,Work(kws),Work(kwd),Work(kwt),isflg,Work(kwc))
  IF ( Iflag==1 ) THEN
    IF ( Nfc/=Nfcc ) CALL SVECS(Ncomp,Nfc,Yhp,Work,Iwork,INHomo,Iflag)
    IF ( Iflag==1 ) THEN
      !
      !- *********************************************************************
      !     DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE INTEGRATED,
      !     INITIALIZE VARIABLES FOR AUXILIARY INITIAL VALUE PROBLEM AND
      !     STORE INITIAL CONDITIONS.
      !
      NEQ = Ncomp*Nfc
      IF ( INHomo==1 ) NEQ = NEQ + Ncomp
      IVP = 0
      IF ( NEQivp/=0 ) THEN
        IVP = NEQ
        NEQ = NEQ + NEQivp
        nfcp2 = nfcp1
        IF ( INHomo==1 ) nfcp2 = nfcp1 + 1
        DO k = 1, NEQivp
          Yhp(k,nfcp2) = Alpha(Nic+k)
        END DO
      END IF
      CALL STOR1(U(:,1,1),Yhp(:,1),V(:,1),Yhp(:,nfcp1),0,NDIsk,NTApe)
      !
      !- *********************************************************************
      !     SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND
      !     SAVE INITIAL CONDITIONS IN CASE A RESTART IS NECESSARY.
      !
      NSWot = 1
      KNSwot = 0
      LOTjp = 1
      TND = LOG10(10.*TOL)
      PWCnd = LOG10(SQRT(TOL))
      X = XBEg
      PX = X
      XOT = XENd
      XOP = X
      KOP = 1
      CALL STWAY(U(:,1,1),V(:,1),Yhp(:,1),0,Stowa)
      !
      !- *********************************************************************
      !- ******* FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS **********
      !- *********************************************************************
      !
      CALL RKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,W,&
        S,Stowa,G,Work,Iwork,Nfcc)
      IF ( Iflag==0.AND.ICOco/=0 ) THEN
        !
        !- *********************************************************************
        !- *************** BACKWARD SWEEP TO OBTAIN SOLUTION *******************
        !- *********************************************************************
        !
        !     CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
        !
        !   FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO READ  U  AND  V
        !   AT THE LAST OUTPUT POINT, SINCE THE LOCAL COPY OF EACH STILL EXISTS.
        !
        kod = 1
        IF ( NDIsk==0 ) kod = Nxpts
        i1 = 1 + Nfcc*Nfcc
        i2 = i1 + Nfcc
        CALL SCOEF(U(1,1,kod),V(1,kod),Ncomp,Nrowb,Nfc,B,Beta,Coef,&
          INHomo,RE,AE,Work,Work(i1),Work(i2),Iwork,Iflag,Nfcc)
        !
        !- *********************************************************************
        !     CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING BACKWARDS.
        !     AS WE RECUR BACKWARDS FROM XFINAL TO XINITIAL WE MUST CALCULATE
        !     NEW SUPERPOSITION COEFFICIENTS EACH TIME WE CROSS A POINT OF
        !     ORTHONORMALIZATION.
        !
        k = NUMort
        ncomp2 = Ncomp/2
        ic = 1
        IF ( Nfc/=Nfcc ) ic = 2
        DO j = 1, Nxpts
          kpts = Nxpts - j + 1
          kod = kpts
          IF ( NDIsk==1 ) kod = 1
          DO WHILE ( k/=0 )
            IF ( XENd>XBEg.AND.Xpts(kpts)>=Z(k) ) EXIT
            IF ( XENd<XBEg.AND.Xpts(kpts)<=Z(k) ) EXIT
            non = k
            IF ( NDIsk/=0 ) THEN
              non = 1
              BACKSPACE NTApe
              READ (NTApe) (Ip(i,1),i=1,Nfcc), (P(i,1),i=1,Ntp)
              BACKSPACE NTApe
            END IF
            IF ( INHomo==1 ) THEN
              IF ( NDIsk/=0 ) THEN
                BACKSPACE NTApe
                READ (NTApe) (W(i,1),i=1,Nfcc)
                BACKSPACE NTApe
              END IF
              DO n = 1, Nfcc
                Coef(n) = Coef(n) - W(n,non)
              END DO
            END IF
            CALL BKSOL(Nfcc,P(1,non),Coef)
            DO m = 1, Nfcc
              Work(m) = Coef(m)
            END DO
            DO m = 1, Nfcc
              l = Ip(m,non)
              Coef(l) = Work(m)
            END DO
            k = k - 1
          END DO
          IF ( NDIsk/=0 ) THEN
            BACKSPACE NTApe
            READ (NTApe) (V(i,1),i=1,Ncomp), ((U(i,m,1),i=1,Ncomp),m=1,Nfc)
            BACKSPACE NTApe
          END IF
          DO n = 1, Ncomp
            Y(n,kpts) = V(n,kod) + DOT_PRODUCT(U(n,1:Nfc,kod),Coef(1:ic*Nfc:ic))
          END DO
          IF ( Nfc/=Nfcc ) THEN
            DO n = 1, ncomp2
              nn = ncomp2 + n
              Y(n,kpts) = Y(n,kpts) - DOT_PRODUCT(U(nn,1:Nfc,kod),Coef(2:2*Nfc:2))
              Y(nn,kpts) = Y(nn,kpts) + DOT_PRODUCT(U(n,1:Nfc,kod),Coef(2:2*Nfc:2))
            END DO
          END IF
        END DO
      END IF
    ELSE
      Iflag = -5
    END IF
  ELSE
    Iflag = -4
  END IF
  !
  !- *********************************************************************
  !
  Mxnon = NUMort
END SUBROUTINE BVPOR
