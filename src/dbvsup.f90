!DECK DBVSUP
SUBROUTINE DBVSUP(Y,Nrowy,Ncomp,Xpts,Nxpts,A,Nrowa,Alpha,Nic,B,Nrowb,Beta,&
    Nfc,Igofx,Re,Ae,Iflag,Work,Ndw,Iwork,Ndiw,Neqivp)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DBVSUP
  !***PURPOSE  Solve a linear two-point boundary value problem using
  !            superposition coupled with an orthonormalization procedure
  !            and a variable-step integration scheme.
  !***LIBRARY   SLATEC
  !***CATEGORY  I1B1
  !***TYPE      DOUBLE PRECISION (BVSUP-S, DBVSUP-D)
  !***KEYWORDS  ORTHONORMALIZATION, SHOOTING,
  !             TWO-POINT BOUNDARY VALUE PROBLEM
  !***AUTHOR  Scott, M. R., (SNLA)
  !           Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !
  !     Subroutine DBVSUP solves a linear two-point boundary-value problem
  !     of the form
  !                        DY/DX = MATRIX(X,U)*Y(X) + G(X,U)
  !                A*Y(XINITIAL) = ALPHA,  B*Y(XFINAL) = BETA
  !
  !     coupled with the solution of the initial value problem
  !
  !                        DU/DX = F(X,U)
  !                      U(XINITIAL) = ETA
  !
  ! **********************************************************************
  !     ABSTRACT
  !        The method of solution uses superposition coupled with an
  !     orthonormalization procedure and a variable-step integration
  !     scheme.  Each time the superposition solutions start to
  !     lose their numerical linear independence, the vectors are
  !     reorthonormalized before integration proceeds.  The underlying
  !     principle of the algorithm is then to piece together the
  !     intermediate (orthogonalized) solutions, defined on the various
  !     subintervals, to obtain the desired solutions.
  !
  ! **********************************************************************
  !     INPUT to DBVSUP
  ! **********************************************************************
  !
  !     NROWY = actual row dimension of Y in calling program.
  !             NROWY must be .GE. NCOMP
  !
  !     NCOMP = number of components per solution vector.
  !             NCOMP is equal to number of original differential
  !             equations.  NCOMP = NIC + NFC.
  !
  !     XPTS = desired output points for solution. They must be monotonic.
  !            XINITIAL = XPTS(1)
  !            XFINAL = XPTS(NXPTS)
  !
  !     NXPTS = number of output points.
  !
  !     A(NROWA,NCOMP) = boundary condition matrix at XINITIAL
  !                      must be contained in (NIC,NCOMP) sub-matrix.
  !
  !     NROWA = actual row dimension of A in calling program,
  !             NROWA must be .GE. NIC.
  !
  !     ALPHA(NIC+NEQIVP) = boundary conditions at XINITIAL.
  !                         If NEQIVP .GT. 0 (see below), the boundary
  !                         conditions at XINITIAL for the initial value
  !                         equations must be stored starting in
  !                         position (NIC + 1) of ALPHA.
  !                         Thus,  ALPHA(NIC+K) = ETA(K).
  !
  !     NIC = number of boundary conditions at XINITIAL.
  !
  !     B(NROWB,NCOMP) = boundary condition matrix at XFINAL.
  !                      Must be contained in (NFC,NCOMP) sub-matrix.
  !
  !     NROWB = actual row dimension of B in calling program,
  !             NROWB must be .GE. NFC.
  !
  !     BETA(NFC) = boundary conditions at XFINAL.
  !
  !     NFC = number of boundary conditions at XFINAL.
  !
  !     IGOFX =0 -- The inhomogeneous term G(X) is identically zero.
  !           =1 -- The inhomogeneous term G(X) is not identically zero.
  !                 (if IGOFX=1, then Subroutine DGVEC (or DUVEC) must be
  !                  supplied).
  !
  !     RE = relative error tolerance used by the integrator.
  !          (see one of the integrators)
  !
  !     AE = absolute error tolerance used by the integrator.
  !          (see one of the integrators)
  ! **NOTE-  RE and AE should not both be zero.
  !
  !     IFLAG = a status parameter used principally for output.
  !             However, for efficient solution of problems which
  !             are originally defined as COMPLEX*16 valued (but
  !             converted to double precision systems to use this code),
  !             the user must set IFLAG=13 on input. See the comment
  !             below for more information on solving such problems.
  !
  !     WORK(NDW) = floating point array used for internal storage.
  !
  !     NDW = actual dimension of work array allocated by user.
  !           An estimate for NDW can be computed from the following
  !            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of
  !                                           orthonormalizations/8)
  !           For the disk or tape storage mode,
  !            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130
  !  However, when the ADAMS integrator is to be used, the estimates are
  !            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of
  !                                           orthonormalizations/8)
  !    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130  , respectively.
  !
  !     IWORK(NDIW) = integer array used for internal storage.
  !
  !     NDIW = actual dimension of IWORK array allocated by user.
  !            An estimate for NDIW can be computed from the following
  !            NDIW = 68 + NCOMP * (1 + expected number of
  !                                            orthonormalizations)
  ! **NOTE --  the amount of storage required is problem dependent and may
  !            be difficult to predict in advance.  Experience has shown
  !            that for most problems 20 or fewer orthonormalizations
  !            should suffice. If the problem cannot be completed with the
  !            allotted storage, then a message will be printed which
  !            estimates the amount of storage necessary. In any case, the
  !            user can examine the IWORK array for the actual storage
  !            requirements, as described in the output information below.
  !
  !     NEQIVP = number of auxiliary initial value equations being added
  !              to the boundary value problem.
  ! **NOTE -- Occasionally the coefficients  matrix  and/or  G  may be
  !           functions which depend on the independent variable  X  and
  !           on  U, the solution of an auxiliary initial value problem.
  !           In order to avoid the difficulties associated with
  !           interpolation, the auxiliary equations may be solved
  !           simultaneously with the given boundary value problem.
  !           This initial value problem may be linear or nonlinear.
  !                 See SAND77-1328 for an example.
  !
  !
  !     The user must supply subroutines DFMAT, DGVEC, DUIVP and DUVEC,
  !     when needed (they must be so named), to evaluate the derivatives
  !     as follows
  !
  !        A. DFMAT must be supplied.
  !
  !              SUBROUTINE DFMAT(X,Y,YP)
  !              X = independent variable (input to DFMAT)
  !              Y = dependent variable vector (input to DFMAT)
  !              YP = DY/DX = derivative vector (output from DFMAT)
  !
  !            Compute the derivatives for the homogeneous problem
  !              YP(I) = DY(I)/DX = MATRIX(X) * Y(I) , I = 1,...,NCOMP
  !
  !            When (NEQIVP .GT. 0) and  matrix  is dependent on  U  as
  !            well as on  X, the following common statement must be
  !            included in DFMAT
  !                    COMMON /DMLIVP/ NOFST
  !            for convenience, the  U  vector is stored at the bottom
  !            of the  Y  array.  Thus, during any call to DFMAT,
  !            U(I) is referenced by  Y(NOFST + I).
  !
  !
  !            Subroutine DBVDER calls DFMAT NFC times to evaluate the
  !            homogeneous equations and, if necessary, it calls DFMAT
  !            once in evaluating the particular solution. since X remains
  !            unchanged in this sequence of calls it is possible to
  !            realize considerable computational savings for complicated
  !            and expensive evaluations of the matrix entries. To do this
  !            the user merely passes a variable, say XS, via common where
  !            XS is defined in the main program to be any value except
  !            the initial X. Then the non-constant elements of matrix(x)
  !            appearing in the differential equations need only be
  !            computed if X is unequal to XS, whereupon XS is reset to X.
  !
  !
  !        B. If  NEQIVP .GT. 0,  DUIVP must also be supplied.
  !
  !              SUBROUTINE DUIVP(X,U,UP)
  !              X = independent variable (input to DUIVP)
  !              U = dependent variable vector (input to DUIVP)
  !              UP = DU/DX = derivative vector (output from DUIVP)
  !
  !            Compute the derivatives for the auxiliary initial value eqs
  !              UP(I) = DU(I)/DX, I = 1,...,NEQIVP.
  !
  !            Subroutine DBVDER calls DUIVP once to evaluate the
  !            derivatives for the auxiliary initial value equations.
  !
  !
  !        C. If  NEQIVP = 0  and  IGOFX = 1,  DGVEC must be supplied.
  !
  !              SUBROUTINE DGVEC(X,G)
  !              X = independent variable (input to DGVEC)
  !              G = vector of inhomogeneous terms G(X) (output from
  !              DGVEC)
  !
  !            Compute the inhomogeneous terms G(X)
  !                G(I) = G(X) values for I = 1,...,NCOMP.
  !
  !            Subroutine DBVDER calls DGVEC in evaluating the particular
  !            solution provided G(X) is not identically zero. Thus, when
  !            IGOFX=0, the user need not write a DGVEC subroutine. Also,
  !            the user does not have to bother with the computational
  !            savings scheme for DGVEC as this is automatically achieved
  !            via the DBVDER subroutine.
  !
  !
  !        D. If  NEQIVP .GT. 0  and  IGOFX = 1,  DUVEC must be supplied.
  !
  !             SUBROUTINE DUVEC(X,U,G)
  !             X = independent variable (input to DUVEC)
  !             U = dependent variable vector from the auxiliary initial
  !                 value problem    (input to DUVEC)
  !             G = array of inhomogeneous terms G(X,U)(output from DUVEC)
  !
  !            Compute the inhomogeneous terms G(X,U)
  !                G(I) = G(X,U) values for I = 1,...,NCOMP.
  !
  !            Subroutine DBVDER calls DUVEC in evaluating the particular
  !            solution provided G(X,U) is not identically zero.  Thus,
  !            when IGOFX=0, the user need not write a DUVEC subroutine.
  !
  !
  !
  !     The following is optional input to DBVSUP to give user more
  !     flexibility in use of code.  See SAND75-0198, SAND77-1328,
  !     SAND77-1690, SAND78-0522, and SAND78-1501 for more information.
  !
  ! ****CAUTION -- The user must zero out IWORK(1),...,IWORK(15)
  !                prior to calling DBVSUP. These locations define
  !                optional input and must be zero unless set to special
  !                values by the user as described below.
  !
  !     IWORK(1) -- number of orthonormalization points.
  !                 A value need be set only if IWORK(11) = 1
  !
  !     IWORK(9) -- integrator and orthonormalization parameter
  !                 (default value is 1)
  !                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test.
  !                 2 = ADAMS code using GRAM-SCHMIDT test.
  !
  !     IWORK(11) -- orthonormalization points parameter
  !                  (default value is 0)
  !                  0 - orthonormalization points not pre-assigned.
  !                  1 - orthonormalization points pre-assigned in
  !                      the first IWORK(1) positions of work.
  !
  !     IWORK(12) -- storage parameter
  !                  (default value is 0)
  !                  0 - all storage in core.
  !                  LUN - homogeneous and inhomogeneous solutions at
  !                      output points and orthonormalization information
  !                      are stored on disk.  The logical unit number to
  !                      be used for disk I/O (NTAPE) is set to IWORK(12).
  !
  !     WORK(1),... -- pre-assigned orthonormalization points, stored
  !                    monotonically, corresponding to the direction
  !                    of integration.
  !
  !
  !
  !                 ******************************************************
  !                 *** COMPLEX*16 VALUED PROBLEM ***
  !                 ******************************************************
  ! **NOTE***
  !       Suppose the original boundary value problem is NC equations
  !     of the form
  !                   DW/DX = MAT(X,U)*W(X) + H(X,U)
  !                 R*W(XINITIAL)=GAMMA, S*W(XFINAL)=DELTA
  !     where all variables are COMPLEX*16 valued. The DBVSUP code can be
  !     used by converting to a double precision system of size 2*NC. To
  !     solve the larger dimensioned problem efficiently, the user must
  !     initialize IFLAG=13 on input and order the vector components
  !     according to Y(1)=DOUBLE PRECISION(W(1)),...,Y(NC)=DOUBLE
  !     PRECISION(W(NC)),Y(NC+1)=IMAG(W(1)),...., Y(2*NC)=IMAG(W(NC)).
  !     Then define
  !                        ...............................................
  !                        . DOUBLE PRECISION(MAT)    -IMAG(MAT) .
  !            MATRIX  =   .                         .
  !                        . IMAG(MAT)     DOUBLE PRECISION(MAT) .
  !                        ...............................................
  !
  !     The matrices A,B and vectors G,ALPHA,BETA must be defined
  !     similarly. Further details can be found in SAND78-1501.
  !
  !
  ! **********************************************************************
  !     OUTPUT from DBVSUP
  ! **********************************************************************
  !
  !     Y(NROWY,NXPTS) = solution at specified output points.
  !
  !     IFLAG Output Values
  !            =-5 algorithm ,for obtaining starting vectors for the
  !                special COMPLEX*16 problem structure, was unable to
  !                obtain the initial vectors satisfying the necessary
  !                independence criteria.
  !            =-4 rank of boundary condition matrix A is less than NIC,
  !                as determined by DLSSUD.
  !            =-2 invalid input parameters.
  !            =-1 insufficient number of storage locations allocated for
  !                WORK or IWORK.
  !
  !            =0 indicates successful solution.
  !
  !            =1 a computed solution is returned but uniqueness of the
  !               solution of the boundary-value problem is questionable.
  !               For an eigenvalue problem, this should be treated as a
  !               successful execution since this is the expected mode
  !               of return.
  !            =2 a computed solution is returned but the existence of the
  !               solution to the boundary-value problem is questionable.
  !            =3 a nontrivial solution approximation is returned although
  !               the boundary condition matrix B*Y(XFINAL) is found to be
  !               nonsingular (to the desired accuracy level) while the
  !               right hand side vector is zero. To eliminate this type
  !               of return, the accuracy of the eigenvalue parameter
  !               must be improved.
  !            ***NOTE-We attempt to diagnose the correct problem behavior
  !               and report possible difficulties by the appropriate
  !               error flag.  However, the user should probably resolve
  !               the problem using smaller error tolerances and/or
  !               perturbations in the boundary conditions or other
  !               parameters. This will often reveal the correct
  !               interpretation for the problem posed.
  !
  !            =13 maximum number of orthonormalizations attained before
  !                reaching XFINAL.
  !            =20-flag from integrator (DDERKF or DDEABM) values can
  !                range from 21 to 25.
  !            =30 solution vectors form a dependent set.
  !
  !     WORK(1),...,WORK(IWORK(1)) = orthonormalization points
  !                                  determined by DBVPOR.
  !
  !     IWORK(1) = number of orthonormalizations performed by DBVPOR.
  !
  !     IWORK(2) = maximum number of orthonormalizations allowed as
  !                calculated from storage allocated by user.
  !
  !     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   give information about
  !                actual storage requirements for WORK and IWORK
  !                arrays.  In particular,
  !                       required storage for  work array is
  !        IWORK(3) + IWORK(4)*(expected number of orthonormalizations)
  !
  !                       required storage for IWORK array is
  !        IWORK(5) + IWORK(6)*(expected number of orthonormalizations)
  !
  !     IWORK(8) = final value of exponent parameter used in tolerance
  !                test for orthonormalization.
  !
  !     IWORK(16) = number of independent vectors returned from DMGSBV.
  !                It is only of interest when IFLAG=30 is obtained.
  !
  !     IWORK(17) = numerically estimated rank of the boundary
  !                 condition matrix defined from B*Y(XFINAL)
  !
  ! **********************************************************************
  !
  !     Necessary machine constants are defined in the Function
  !     Routine D1MACH. The user must make sure that the values
  !     set in D1MACH are relevant to the computer being used.
  !
  ! **********************************************************************
  ! **********************************************************************
  !
  !***REFERENCES  M. R. Scott and H. A. Watts, SUPORT - a computer code
  !                 for two-point boundary-value problems via
  !                 orthonormalization, SIAM Journal of Numerical
  !                 Analysis 14, (1977), pp. 40-70.
  !               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
  !                 of SUPORT, a linear boundary value problem solver
  !                 Part I - pre-assigning orthonormalization points,
  !                 auxiliary initial value problem, disk or tape storage,
  !                 Report SAND77-1328, Sandia Laboratories, Albuquerque,
  !                 New Mexico, 1977.
  !               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
  !                 of SUPORT, a linear boundary value problem solver
  !                 Part II - inclusion of an Adams integrator, Report
  !                 SAND77-1690, Sandia Laboratories, Albuquerque,
  !                 New Mexico, 1977.
  !               M. E. Lord and H. A. Watts, Modifications of SUPORT,
  !                 a linear boundary value problem solver Part III -
  !                 orthonormalization improvements, Report SAND78-0522,
  !                 Sandia Laboratories, Albuquerque, New Mexico, 1978.
  !               H. A. Watts, M. R. Scott and M. E. Lord, Computational
  !                 solution of complex*16 valued boundary problems,
  !                 Report SAND78-1501, Sandia Laboratories,
  !                 Albuquerque, New Mexico, 1978.
  !***ROUTINES CALLED  DEXBVP, DMACON, XERMSG
  !***COMMON BLOCKS    DML15T, DML17B, DML18J, DML5MC, DML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   890921  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Convert XERRWV calls to XERMSG calls, remove some extraneous
  !           comments.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DBVSUP
  ! **********************************************************************
  !
  INTEGER ICOco, Iflag, Igofx, IGOfxd, INDpvt, INFo, INHomo, INTeg ,&
    is, ISTkop, IVP, Iwork(*), j, k, K1, K10, K11, K2, K3 ,&
    K4, K5, K6, K7, K8, K9, kkkcoe, kkkcof, kkkg, KKKint ,&
    kkks, kkksto, kkksud, kkksvc, kkku, kkkv, kkkws, kkkyhp ,&
    KKKzpw, KNSwot, KOP, kpts, L1, L2, lllcof, LLLint, lllip ,&
    llliws, lllsud, lllsvc, LOTjp, LPAr, MNSwot, MXNon ,&
    mxnoni, mxnonr, Ncomp, NCOmpd, ndeq, NDIsk, Ndiw, Ndw ,&
    NEEdiw, NEEdw, NEQ, NEQivd, Neqivp, Nfc, NFCc, NFCd, Nic ,&
    NICd, nitemp, non, NOPg, NPS, Nrowa, Nrowb, Nrowy ,&
    nrtemp, NSWot, NTApe, NTP, NUMort, Nxpts, NXPtsd, nxptsm
  REAL(8) :: A(Nrowa,*), Ae, AED, Alpha(*), B(Nrowb,*), Beta(*) ,&
    C, EPS, FOUru, PWCnd, PX, Re, RED, SQOvfl, SRU ,&
    TND, TOL, TWOu, URO, Work(Ndw), X, XBEg, XENd ,&
    XOP, XOT, Xpts(*), XSAv, Y(Nrowy,*)
  CHARACTER(8) :: xern1, xern2, xern3, xern4
  !
  !     ******************************************************************
  !         THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE
  !         DBVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN
  !         THE CALLING PROGRAM.
  !
  COMMON /DML8SZ/ C, XSAv, IGOfxd, INHomo, IVP, NCOmpd, NFCd
  !
  !     ******************************************************************
  !         THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE
  !         ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE
  !
  COMMON /DML18J/ AED, RED, TOL, NXPtsd, NICd, NOPg, MXNon, NDIsk ,&
    NTApe, NEQ, INDpvt, INTeg, NPS, NTP, NEQivd ,&
    NUMort, NFCc, ICOco
  COMMON /DML17B/ KKKzpw, NEEdw, NEEdiw, K1, K2, K3, K4, K5, K6 ,&
    K7, K8, K9, K10, K11, L1, L2, KKKint, LLLint
  !
  !     ******************************************************************
  !         THIS COMMON BLOCK IS USED IN SUBROUTINES DBVSUP,DBVPOR,DRKFAB,
  !         DREORT, AND DSTWAY. IT CONTAINS INFORMATION NECESSARY
  !         FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP
  !         RESTARTING CAPABILITY.
  !
  COMMON /DML15T/ PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, INFo(15)&
    , ISTkop, KNSwot, KOP, LOTjp, MNSwot, NSWot
  !
  !     ******************************************************************
  !         THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS
  !         USED BY THE CODE
  !
  COMMON /DML5MC/ URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  !
  !      *****************************************************************
  !          SET UP MACHINE DEPENDENT CONSTANTS.
  !
  !***FIRST EXECUTABLE STATEMENT  DBVSUP
  CALL DMACON
  !
  !                       ************************************************
  !                           TEST FOR INVALID INPUT
  !
  IF ( Nrowy>=Ncomp ) THEN
    IF ( Ncomp==Nic+Nfc ) THEN
      IF ( Nxpts>=2 ) THEN
        IF ( Nic>0 ) THEN
          IF ( Nrowa>=Nic ) THEN
            IF ( Nfc>0 ) THEN
              IF ( Nrowb>=Nfc ) THEN
                IF ( Igofx>=0.AND.Igofx<=1 ) THEN
                  IF ( Re>=0.0D0 ) THEN
                    IF ( Ae>=0.0D0 ) THEN
                      IF ( Re/=0.0D0.OR.Ae/=0.0D0 ) THEN
                        !                          BEGIN BLOCK PERMITTING ...EXITS TO 70
                        is = 1
                        IF ( Xpts(Nxpts)<Xpts(1) ) is = 2
                        nxptsm = Nxpts - 1
                        DO k = 1, nxptsm
                          IF ( is==2 ) THEN
                            !                          .........EXIT
                            IF ( Xpts(k)<=Xpts(k+1) ) GOTO 100
                            !                          .........EXIT
                          ELSEIF ( Xpts(k+1)<=Xpts(k) ) THEN
                            GOTO 100
                          ENDIF
                        ENDDO
                        !
                        !                             ******************************************
                        !                                 CHECK FOR DISK STORAGE
                        !
                        kpts = Nxpts
                        NDIsk = 0
                        IF ( Iwork(12)/=0 ) THEN
                          NTApe = Iwork(12)
                          kpts = 1
                          NDIsk = 1
                        ENDIF
                        !
                        !                             ******************************************
                        !                                 SET INTEG PARAMETER ACCORDING TO
                        !                                 CHOICE OF INTEGRATOR.
                        !
                        INTeg = 1
                        IF ( Iwork(9)==2 ) INTeg = 2
                        !
                        !                             ******************************************
                        !                                 COMPUTE INHOMO
                        !
                        !                 ............EXIT
                        IF ( Igofx==1 ) GOTO 300
                        DO j = 1, Nic
                          !                 ...............EXIT
                          IF ( Alpha(j)/=0.0D0 ) GOTO 300
                        ENDDO
                        DO j = 1, Nfc
                          !                    ............EXIT
                          IF ( Beta(j)/=0.0D0 ) GOTO 200
                        ENDDO
                        INHomo = 3
                        !              ...............EXIT
                        GOTO 400
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  100 CONTINUE
  IFlag = -2
  !     ..................EXIT
  RETURN
  200  INHomo = 2
  !              ......EXIT
  GOTO 400
  300  INHomo = 1
  !
  !              *********************************************************
  !                  TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN
  !                  SOLVING A COMPLEX*16 VALUED PROBLEM,WE INTRODUCE
  !                  NFCC=NFC WHILE CHANGING THE INTERNAL VALUE OF NFC
  !
  400  NFCc = Nfc
  IF ( Iflag==13 ) Nfc = Nfc/2
  !
  !              *********************************************************
  !                  DETERMINE NECESSARY STORAGE REQUIREMENTS
  !
  !              FOR BASIC ARRAYS IN DBVPOR
  kkkyhp = Ncomp*(Nfc+1) + Neqivp
  kkku = Ncomp*Nfc*kpts
  kkkv = Ncomp*kpts
  kkkcoe = NFCc
  kkks = Nfc + 1
  kkksto = Ncomp*(Nfc+1) + Neqivp + 1
  kkkg = Ncomp
  !
  !              FOR ORTHONORMALIZATION RELATED MATTERS
  NTP = (NFCc*(NFCc+1))/2
  KKKzpw = 1 + NTP + NFCc
  lllip = NFCc
  !
  !              FOR ADDITIONAL REQUIRED WORK SPACE
  !                (DLSSUD)
  kkksud = 4*Nic + (Nrowa+1)*Ncomp
  lllsud = Nic
  !              (DVECS)
  kkksvc = 1 + 4*NFCc + 2*NFCc**2
  lllsvc = 2*NFCc
  !
  ndeq = Ncomp*Nfc + Neqivp
  IF ( INHomo==1 ) ndeq = ndeq + Ncomp
  IF ( INTeg==2 ) THEN
    !              (DDEABM)
    KKKint = 130 + 21*ndeq
    LLLint = 51
  ELSE
    !              (DDERKF)
    KKKint = 33 + 7*ndeq
    LLLint = 34
  ENDIF
  !
  !              (COEF)
  kkkcof = 5*NFCc + NFCc**2
  lllcof = 3 + NFCc
  !
  kkkws = MAX(kkksud,kkksvc,KKKint,kkkcof)
  llliws = MAX(lllsud,lllsvc,LLLint,lllcof)
  !
  NEEdw = kkkyhp + kkku + kkkv + kkkcoe + kkks + kkksto + kkkg + KKKzpw +&
    kkkws
  NEEdiw = 17 + lllip + llliws
  !              *********************************************************
  !                  COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS
  !                  WITH THE ALLOTTED STORAGE
  !
  Iwork(3) = NEEdw
  Iwork(4) = KKKzpw
  Iwork(5) = NEEdiw
  Iwork(6) = lllip
  nrtemp = Ndw - NEEdw
  nitemp = Ndiw - NEEdiw
  !           ...EXIT
  IF ( nrtemp>=0 ) THEN
    !           ...EXIT
    IF ( nitemp>=0 ) THEN
      !
      IF ( NDIsk==0 ) THEN
        !
        mxnonr = nrtemp/KKKzpw
        mxnoni = nitemp/lllip
        MXNon = MIN(mxnonr,mxnoni)
        non = MXNon
      ELSE
        non = 0
        MXNon = nrtemp
      ENDIF
      !
      Iwork(2) = MXNon
      !
      !              *********************************************************
      !                  CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS
      !
      NOPg = 0
      !        ......EXIT
      IF ( Iwork(11)/=1 ) GOTO 500
      IF ( MXNon>=Iwork(1) ) THEN
        NOPg = 1
        MXNon = Iwork(1)
        Work(MXNon+1) = 2.0D0*Xpts(Nxpts) - Xpts(1)
        !        .........EXIT
        GOTO 500
      ENDIF
    ENDIF
  ENDIF
  !
  Iflag = -1
  IF ( NDIsk/=1 ) THEN
    WRITE (xern1,'(I8)') NEEdw
    WRITE (xern2,'(I8)') KKKzpw
    WRITE (xern3,'(I8)') NEEdiw
    WRITE (xern4,'(I8)') lllip
    CALL XERMSG('SLATEC','DBVSUP','REQUIRED STORAGE FOR WORK ARRAY IS '//&
      xern1//' + '//xern2//&
      '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS) $$'//&
      'REQUIRED STORAGE FOR IWORK ARRAY IS '//xern3//' + '//&
      xern4//'*(EXPECTED NUMBER OF ORTHONORMALIZATIONS)',1,0)
  ELSE
    WRITE (xern1,'(I8)') NEEdw
    WRITE (xern2,'(I8)') NEEdiw
    CALL XERMSG('SLATEC','DBVSUP','REQUIRED STORAGE FOR WORK ARRAY IS '//&
      xern1//' + NUMBER OF ORTHONOMALIZATIONS. $$'//&
      'REQUIRED STORAGE FOR IWORK ARRAY IS '//xern2,1,0)
  ENDIF
  RETURN
  !
  !        ***************************************************************
  !            ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS
  !
  !         (Z)
  500  K1 = 1 + (MXNon+1)
  !        (P)
  K2 = K1 + NTP*(non+1)
  !        (W)
  K3 = K2 + NFCc*(non+1)
  !        (YHP)
  K4 = K3 + kkkyhp
  !        (U)
  K5 = K4 + kkku
  !        (V)
  K6 = K5 + kkkv
  !        (COEF)
  K7 = K6 + kkkcoe
  !        (S)
  K8 = K7 + kkks
  !        (STOWA)
  K9 = K8 + kkksto
  !        (G)
  K10 = K9 + kkkg
  K11 = K10 + kkkws
  !                  REQUIRED ADDITIONAL DOUBLE PRECISION WORK SPACE
  !                  STARTS AT WORK(K10) AND EXTENDS TO WORK(K11-1)
  !
  !           FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL
  !           INPUT AND OUTPUT ITEMS
  !        (IP)
  L1 = 18 + NFCc*(non+1)
  L2 = L1 + llliws
  !                   REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1)
  !                   AND EXTENDS TO IWORK(L2-1)
  !
  !        ***************************************************************
  !            SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION
  !
  NPS = 0
  IF ( Iwork(10)==1 ) NPS = 1
  !
  !        ***************************************************************
  !            SET PIVOTING PARAMETER
  !
  INDpvt = 0
  IF ( Iwork(15)==1 ) INDpvt = 1
  !
  !        ***************************************************************
  !            SET OTHER COMMON BLOCK PARAMETERS
  !
  NFCd = Nfc
  NCOmpd = Ncomp
  IGOfxd = Igofx
  NXPtsd = Nxpts
  NICd = Nic
  RED = Re
  AED = Ae
  NEQivd = Neqivp
  MNSwot = 20
  IF ( Iwork(13)==-1 ) MNSwot = MAX(1,Iwork(14))
  XBEg = Xpts(1)
  XENd = Xpts(Nxpts)
  XSAv = XENd
  ICOco = 1
  IF ( INHomo==3.AND.NOPg==1 ) Work(MXNon+1) = XENd
  !
  !        ***************************************************************
  !
  CALL DEXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  Nfc = NFCc
  Iwork(17) = Iwork(L1)
  RETURN
END SUBROUTINE DBVSUP
