!*==DBVPOR.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBVPOR
SUBROUTINE DBVPOR(Y,Nrowy,Ncomp,Xpts,Nxpts,A,Nrowa,Alpha,Nic,B,Nrowb,Beta,&
    Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,W,Niv,Yhp,U,V,Coef,S,Stowa,G,&
    Work,Iwork,Nfcc)
  IMPLICIT NONE
  !*--DBVPOR7
  !***BEGIN PROLOGUE  DBVPOR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BVPOR-S, DBVPOR-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !     INPUT to DBVPOR    (items not defined in DBVSUP comments)
  ! **********************************************************************
  !
  !     NOPG = 0 -- orthonormalization points not pre-assigned
  !          = 1 -- orthonormalization points pre-assigned
  !
  !     MXNON = maximum number of orthogonalizations allowed.
  !
  !     NDISK = 0 -- in-core storage
  !           = 1 -- disk storage.  Value of NTAPE in data statement
  !                  is set to 13.  If another value is desired,
  !                  the data statement must be changed.
  !
  !     INTEG = type of integrator and associated test to be used
  !             to determine when to orthonormalize.
  !
  !             1 -- use GRAM-SCHMIDT test and DDERKF
  !             2 -- use GRAM-SCHMIDT test and DDEABM
  !
  !     TOL = tolerance for allowable error in orthogonalization test.
  !
  !     NPS = 0 normalize particular solution to unit length at each
  !             point of orthonormalization.
  !         = 1 do not normalize particular solution.
  !
  !     NTP = must be .GE. NFC*(NFC+1)/2.
  !
  !     NFCC = 2*NFC for special treatment of a COMPLEX*16 valued problem
  !
  !     ICOCO = 0 skip final computations (superposition coefficients
  !               and, hence, boundary problem solution)
  !           = 1 calculate superposition coefficients and obtain
  !               solution to the boundary value problem
  !
  ! **********************************************************************
  !     OUTPUT from DBVPOR
  ! **********************************************************************
  !
  !     Y(NROWY,NXPTS) = solution at specified output points.
  !
  !     MXNON = number of orthonormalizations performed by DBVPOR.
  !
  !     Z(MXNON+1) = locations of orthonormalizations performed by DBVPOR.
  !
  !     NIV = number of independent vectors returned from DMGSBV. Normally
  !           this parameter will be meaningful only when DMGSBV returns
  !           with MFLAG = 2.
  !
  ! **********************************************************************
  !
  !     The following variables are in the argument list because of
  !     variable dimensioning.  In general, they contain no information of
  !     use to the user.  The amount of storage set aside by the user must
  !     be greater than or equal to that indicated by the dimension
  !     statements.  For the disk storage mode, NON = 0 and KPTS = 1,
  !     while for the in-core storage mode, NON = MXNON and KPTS = NXPTS.
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
  ! **********************************************************************
  !     SUBROUTINES used by DBVPOR
  !         DLSSUD -- solves an underdetermined system of linear
  !                   equations.  This routine is used to get a full
  !                   set of initial conditions for integration.
  !                   Called by DBVPOR.
  !
  !         DVECS -- obtains starting vectors for special treatment
  !                   of COMPLEX*16 valued problems, called by DBVPOR.
  !
  !         DRKFAB -- routine which conducts integration using DDERKF or
  !                   DDEABM.
  !
  !         DSTWAY -- storage for backup capability, called by
  !                   DBVPOR and DREORT.
  !
  !         DSTOR1 -- storage at output points, called by DBVPOR,
  !                   DRKFAB, DREORT and DSTWAY.
  !
  !         DDOT -- single precision vector inner product routine,
  !                   called by DBVPOR, DCOEF, DLSSUD, DMGSBV,
  !                   DBKSOL, DREORT and DPRVEC.
  !         ** NOTE **
  !         a considerable improvement in speed can be achieved if a
  !         machine language version is used for DDOT.
  !
  !         DCOEF -- computes the superposition constants from the
  !                   boundary conditions at XFINAL.
  !
  !         DBKSOL -- solves an upper triangular set of linear equations.
  !
  ! **********************************************************************
  !
  !***SEE ALSO  DBVSUP
  !***ROUTINES CALLED  DBKSOL, DCOEF, DDOT, DLSSUD, DRKFAB, DSTOR1,
  !                    DSTWAY, DVECS
  !***COMMON BLOCKS    DML15T, DML18J, DML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DBVPOR
  !
  REAL(8) :: DDOT
  INTEGER i, i1, i2, ic, ICOco, Iflag, IGOfx, INDpvt, INFo, &
    INHomo, INTeg, ira, isflg, ISTkop, IVP, j, k, KNSwot, &
    kod, KOP, kpts, kwc, kwd, kws, kwt, l, LOTjp, m, &
    MNSwot, Mxnon, MXNond, n, Ncomp, ncomp2, NCOmpd, NDIsk, &
    ndw, NEQ, NEQivp, Nfc, Nfcc, NFCcd, NFCd, nfcp1, nfcp2, &
    Nic, NICd, Niv, nn, non, NOPg, NPS, Nrowa, Nrowb, Nrowy, &
    NSWot, NTApe, Ntp, NTPd, NUMort, Nxpts, NXPtsd, Ip(Nfcc,*)&
    , Iwork(*)
  REAL(8) :: A(Nrowa,*), AE, Alpha(*), B(Nrowb,*), Beta(*), C, &
    Coef(*), G(*), P(Ntp,*), PWCnd, PX, RE, S(*), &
    Stowa(*), TND, TOL, U(Ncomp,Nfc,*), V(Ncomp,*), &
    W(Nfcc,*), Work(*), X, XBEg, XENd, XOP, XOT, &
    Xpts(*), XSAv, Y(Nrowy,*), Yhp(Ncomp,*), Z(*)
  !
  !     ******************************************************************
  !
  COMMON /DML8SZ/ C, XSAv, IGOfx, INHomo, IVP, NCOmpd, NFCd
  COMMON /DML15T/ PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, INFo(15)&
    , ISTkop, KNSwot, KOP, LOTjp, MNSwot, NSWot
  COMMON /DML18J/ AE, RE, TOL, NXPtsd, NICd, NOPg, MXNond, NDIsk, &
    NTApe, NEQ, INDpvt, INTeg, NPS, NTPd, NEQivp, &
    NUMort, NFCcd, ICOco
  !
  !      *****************************************************************
  !
  !***FIRST EXECUTABLE STATEMENT  DBVPOR
  nfcp1 = Nfc + 1
  NUMort = 0
  C = 1.0D0
  !
  !     ******************************************************************
  !         CALCULATE INITIAL CONDITIONS WHICH SATISFY
  !                       A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
  !         WHEN NFC .NE. NFCC DLSSUD DEFINES VALUES YHP IN A MATRIX OF
  !         SIZE (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE
  !         ALLOCATION INTO THE U ARRAY. HOWEVER, THIS IS OKAY SINCE
  !         PLENTY OF SPACE IS AVAILABLE IN U AND IT HAS NOT YET BEEN
  !         USED.
  !
  ndw = Nrowa*Ncomp
  kws = ndw + Nic + 1
  kwd = kws + Nic
  kwt = kwd + Nic
  kwc = kwt + Nic
  Iflag = 0
  CALL DLSSUD(A,Yhp(1,Nfcc+1),Alpha,Nic,Ncomp,Nrowa,Yhp,Ncomp,Iflag,1,ira,0,&
    Work(1),Work(ndw+1),Iwork,Work(kws),Work(kwd),Work(kwt),isflg,&
    Work(kwc))
  IF ( Iflag==1 ) THEN
    IF ( Nfc/=Nfcc ) CALL DVECS(Ncomp,Nfc,Yhp,Work,Iwork,INHomo,Iflag)
    IF ( Iflag==1 ) THEN
      !
      !           ************************************************************
      !               DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE
      !               INTEGRATED, INITIALIZE VARIABLES FOR AUXILIARY INITIAL
      !               VALUE PROBLEM AND STORE INITIAL CONDITIONS.
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
        ENDDO
      ENDIF
      CALL DSTOR1(U,Yhp,V,Yhp(1,nfcp1),0,NDIsk,NTApe)
      !
      !           ************************************************************
      !               SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE
      !               AND SAVE INITIAL CONDITIONS IN CASE A RESTART IS
      !               NECESSARY.
      !
      NSWot = 1
      KNSwot = 0
      LOTjp = 1
      TND = LOG10(10.0D0*TOL)
      PWCnd = LOG10(SQRT(TOL))
      X = XBEg
      PX = X
      XOT = XENd
      XOP = X
      KOP = 1
      CALL DSTWAY(U,V,Yhp,0,Stowa)
      !
      !           ************************************************************
      !           ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS
      !           **********
      !           ************************************************************
      !
      CALL DRKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,W,&
        S,Stowa,G,Work,Iwork,Nfcc)
      IF ( Iflag==0.AND.ICOco/=0 ) THEN
        !
        !              *********************************************************
        !              **************** BACKWARD SWEEP TO OBTAIN SOLUTION
        !              *******************
        !              *********************************************************
        !
        !                  CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
        !
        !                FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO
        !                READ  U  AND  V AT THE LAST OUTPUT POINT, SINCE THE
        !                LOCAL COPY OF EACH STILL EXISTS.
        !
        kod = 1
        IF ( NDIsk==0 ) kod = Nxpts
        i1 = 1 + Nfcc*Nfcc
        i2 = i1 + Nfcc
        CALL DCOEF(U(1,1,kod),V(1,kod),Ncomp,Nrowb,Nfc,Nic,B,Beta,Coef,&
          INHomo,RE,AE,Work,Work(i1),Work(i2),Iwork,Iflag,Nfcc)
        !
        !              *********************************************************
        !                  CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING
        !                  BACKWARDS.  AS WE RECUR BACKWARDS FROM XFINAL TO
        !                  XINITIAL WE MUST CALCULATE NEW SUPERPOSITION
        !                  COEFFICIENTS EACH TIME WE CROSS A POINT OF
        !                  ORTHONORMALIZATION.
        !
        k = NUMort
        ncomp2 = Ncomp/2
        ic = 1
        IF ( Nfc/=Nfcc ) ic = 2
        DO j = 1, Nxpts
          kpts = Nxpts - j + 1
          kod = kpts
          IF ( NDIsk==1 ) kod = 1
          !                 ...EXIT
          DO WHILE ( k/=0 )
            !                 ...EXIT
            IF ( XENd>XBEg.AND.Xpts(kpts)>=Z(k) ) EXIT
            !                 ...EXIT
            IF ( XENd<XBEg.AND.Xpts(kpts)<=Z(k) ) EXIT
            non = k
            IF ( NDIsk/=0 ) THEN
              non = 1
              BACKSPACE NTApe
              READ (NTApe) (Ip(i,1),i=1,Nfcc), (P(i,1),i=1,Ntp)
              BACKSPACE NTApe
            ENDIF
            IF ( INHomo==1 ) THEN
              IF ( NDIsk/=0 ) THEN
                BACKSPACE NTApe
                READ (NTApe) (W(i,1),i=1,Nfcc)
                BACKSPACE NTApe
              ENDIF
              DO n = 1, Nfcc
                Coef(n) = Coef(n) - W(n,non)
              ENDDO
            ENDIF
            CALL DBKSOL(Nfcc,P(1,non),Coef)
            DO m = 1, Nfcc
              Work(m) = Coef(m)
            ENDDO
            DO m = 1, Nfcc
              l = Ip(m,non)
              Coef(l) = Work(m)
            ENDDO
            k = k - 1
          ENDDO
          IF ( NDIsk/=0 ) THEN
            BACKSPACE NTApe
            READ (NTApe) (V(i,1),i=1,Ncomp), ((U(i,m,1),i=1,Ncomp),m=1,Nfc)
            BACKSPACE NTApe
          ENDIF
          DO n = 1, Ncomp
            Y(n,kpts) = V(n,kod) + DDOT(Nfc,U(n,1,kod),Ncomp,Coef,ic)
          ENDDO
          IF ( Nfc/=Nfcc ) THEN
            DO n = 1, ncomp2
              nn = ncomp2 + n
              Y(n,kpts) = Y(n,kpts) - DDOT(Nfc,U(nn,1,kod),Ncomp,Coef(2),2)
              Y(nn,kpts) = Y(nn,kpts) + DDOT(Nfc,U(n,1,kod),Ncomp,Coef(2),2)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ELSE
      Iflag = -5
    ENDIF
  ELSE
    Iflag = -4
  ENDIF
  !
  !     ******************************************************************
  !
  Mxnon = NUMort
END SUBROUTINE DBVPOR
