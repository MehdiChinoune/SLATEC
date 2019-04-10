!** WNLSM
SUBROUTINE WNLSM(W,Mdw,Mme,Ma,N,L,Prgopt,X,Rnorm,Mode,Ipivot,Itype,Wd,H,&
    Scalee,Z,Temp,D)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to WNNLS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (WNLSM-S, DWNLSM-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to WNNLS.
  !     The documentation for WNNLS has complete usage instructions.
  !
  !     In addition to the parameters discussed in the prologue to
  !     subroutine WNNLS, the following work arrays are used in
  !     subroutine WNLSM  (they are passed through the calling
  !     sequence from WNNLS for purposes of variable dimensioning).
  !     Their contents will in general be of no interest to the user.
  !
  !         IPIVOT(*)
  !            An array of length N.  Upon completion it contains the
  !         pivoting information for the cols of W(*,*).
  !
  !         ITYPE(*)
  !            An array of length M which is used to keep track
  !         of the classification of the equations.  ITYPE(I)=0
  !         denotes equation I as an equality constraint.
  !         ITYPE(I)=1 denotes equation I as a least squares
  !         equation.
  !
  !         WD(*)
  !            An array of length N.  Upon completion it contains the
  !         dual solution vector.
  !
  !         H(*)
  !            An array of length N.  Upon completion it contains the
  !         pivot scalars of the Householder transformations performed
  !         in the case KRANK.LT.L.
  !
  !         SCALE(*)
  !            An array of length M which is used by the subroutine
  !         to store the diagonal matrix of weights.
  !         These are used to apply the modified Givens
  !         transformations.
  !
  !         Z(*),TEMP(*)
  !            Working arrays of length N.
  !
  !         D(*)
  !            An array of length N that contains the
  !         column scaling for the matrix (E).
  !                                       (A)
  !
  !***
  ! **See also:**  WNNLS
  !***
  ! **Routines called:**  H12, ISAMAX, R1MACH, SASUM, SAXPY, SCOPY, SNRM2,
  !                    SROTM, SROTMG, SSCAL, SSWAP, WNLIT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and revised.  (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Fixed an error message.  (RWC)

  INTEGER Ipivot(*), Itype(*), L, Ma, Mdw, Mme, Mode, N
  REAL D(*), H(*), Prgopt(*), Rnorm, Scalee(*), Temp(*), W(Mdw,*), Wd(*), X(*), Z(*)
  !
  INTEGER, EXTERNAL :: ISAMAX
  REAL, EXTERNAL :: R1MACH, SASUM, SNRM2
  EXTERNAL :: H12, SAXPY, SCOPY, SROTM, SROTMG, SSCAL, SSWAP, WNLIT, XERMSG
  !
  REAL alamda, alpha, alsq, amax, blowup, bnorm, dope(3), eanorm, fac, sm, &
    sparam(5), t, tau, wmax, z2, zz
  INTEGER i, idope(3), imax, isol, itemp, iter, itmax, iwmax, j, jcon, jp, key, &
    krank, l1, last, link, m, me, next, niv, nlink, nopt, nsoln, ntimes
  LOGICAL done, feasbl, hitcon, pos
  !
  REAL, SAVE :: srelpr
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  WNLSM
  !
  !     Initialize variables.
  !     SRELPR is the precision for the particular machine
  !     being used.  This logic avoids resetting it every entry.
  !
  IF ( first ) srelpr = R1MACH(4)
  first = .FALSE.
  !
  !     Set the nominal tolerance used in the code.
  !
  tau = SQRT(srelpr)
  !
  m = Ma + Mme
  me = Mme
  Mode = 2
  !
  !     To process option vector
  !
  fac = 1.E-4
  !
  !     Set the nominal blow up factor used in the code.
  !
  blowup = tau
  !
  !     The nominal column scaling used in the code is
  !     the identity scaling.
  !
  D(1:N) = 1.E0
  !
  !     Define bound for number of options to change.
  !
  nopt = 1000
  !
  !     Define bound for positive value of LINK.
  !
  nlink = 100000
  ntimes = 0
  last = 1
  link = INT( Prgopt(1) )
  IF ( link<=0.OR.link>nlink ) THEN
    CALL XERMSG('SLATEC','WNLSM','WNNLS, THE OPTION VECTOR IS UNDEFINED',3,1)
    RETURN
  END IF
  DO
    !
    IF ( link>1 ) THEN
      ntimes = ntimes + 1
      IF ( ntimes>nopt ) THEN
        CALL XERMSG('SLATEC','WNLSM',&
          'WNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.',3, 1)
        RETURN
      END IF
      !
      key = INT( Prgopt(last+1) )
      IF ( key==6.AND.Prgopt(last+2)/=0.E0 ) THEN
        DO j = 1, N
          t = SNRM2(m,W(1,j),1)
          IF ( t/=0.E0 ) t = 1.E0/t
          D(j) = t
        END DO
      END IF
      !
      IF ( key==7 ) CALL SCOPY(N,Prgopt(last+2),1,D,1)
      IF ( key==8 ) tau = MAX(srelpr,Prgopt(last+2))
      IF ( key==9 ) blowup = MAX(srelpr,Prgopt(last+2))
      !
      next = INT( Prgopt(link) )
      IF ( next<=0.OR.next>nlink ) THEN
        CALL XERMSG('SLATEC','WNLSM',&
          'WNNLS, THE OPTION VECTOR IS UNDEFINED',3,1)
        RETURN
      END IF
      !
      last = link
      link = next
      CYCLE
    END IF
    !
    DO j = 1, N
      CALL SSCAL(m,D(j),W(1,j),1)
    END DO
    !
    !     Process option vector
    !
    done = .FALSE.
    iter = 0
    itmax = 3*(N-L)
    Mode = 0
    nsoln = L
    l1 = MIN(m,L)
    !
    !     Compute scalee factor to apply to equality constraint equations.
    !
    DO j = 1, N
      Wd(j) = SASUM(m,W(1,j),1)
    END DO
    !
    imax = ISAMAX(N,Wd,1)
    eanorm = Wd(imax)
    bnorm = SASUM(m,W(1,N+1),1)
    alamda = eanorm/(srelpr*fac)
    !
    !     Define scaling diagonal matrix for modified Givens usage and
    !     classify equation types.
    !
    alsq = alamda**2
    DO i = 1, m
      !
      !        When equation I is heavily weighted ITYPE(I)=0,
      !        else ITYPE(I)=1.
      !
      IF ( i<=me ) THEN
        t = alsq
        itemp = 0
      ELSE
        t = 1.E0
        itemp = 1
      END IF
      Scalee(i) = t
      Itype(i) = itemp
    END DO
    !
    !     Set the solution vector X(*) to zero and the column interchange
    !     matrix to the identity.
    !
    X(1:N) = 0.E0
    DO i = 1, N
      Ipivot(i) = i
    END DO
    !
    !     Perform initial triangularization in the submatrix
    !     corresponding to the unconstrained variables.
    !     Set first L components of dual vector to zero because
    !     these correspond to the unconstrained variables.
    !
    Wd(1:L) = 0.E0
    !
    !     The arrays IDOPE(*) and DOPE(*) are used to pass
    !     information to WNLIT().  This was done to avoid
    !     a long calling sequence or the use of COMMON.
    !
    idope(1) = me
    idope(2) = nsoln
    idope(3) = l1
    !
    dope(1) = alsq
    dope(2) = eanorm
    dope(3) = tau
    CALL WNLIT(W,Mdw,m,N,L,Ipivot,Itype,H,Scalee,Rnorm,idope,dope,done)
    me = idope(1)
    krank = idope(2)
    niv = idope(3)
    EXIT
  END DO
  !
  !     Perform WNNLS algorithm using the following steps.
  !
  !     Until(DONE)
  !        compute search direction and feasible point
  !        when (HITCON) add constraints
  !        else perform multiplier test and drop a constraint
  !        fin
  !     Compute-Final-Solution
  !
  !     To compute search direction and feasible point,
  !     solve the triangular system of currently non-active
  !     variables and store the solution in Z(*).
  !
  !     To solve system
  !     Copy right hand side into TEMP vector to use overwriting method.
  !
  DO WHILE ( .NOT.(done) )
    isol = L + 1
    IF ( nsoln>=isol ) THEN
      CALL SCOPY(niv,W(1,N+1),1,Temp,1)
      DO j = nsoln, isol, -1
        IF ( j>krank ) THEN
          i = niv - nsoln + j
        ELSE
          i = j
        END IF
        !
        IF ( j>krank.AND.j<=L ) THEN
          Z(j) = 0.E0
        ELSE
          Z(j) = Temp(i)/W(i,j)
          CALL SAXPY(i-1,-Z(j),W(1,j),1,Temp,1)
        END IF
      END DO
    END IF
    !
    !     Increment iteration counter and check against maximum number
    !     of iterations.
    !
    iter = iter + 1
    IF ( iter>itmax ) THEN
      Mode = 1
      done = .TRUE.
    END IF
    !
    !     Check to see if any constraints have become active.
    !     If so, calculate an interpolation factor so that all
    !     active constraints are removed from the basis.
    !
    alpha = 2.E0
    hitcon = .FALSE.
    DO j = L + 1, nsoln
      zz = Z(j)
      IF ( zz<=0.E0 ) THEN
        t = X(j)/(X(j)-zz)
        IF ( t<alpha ) THEN
          alpha = t
          jcon = j
        END IF
        hitcon = .TRUE.
      END IF
    END DO
    !
    !     Compute search direction and feasible point
    !
    IF ( hitcon ) THEN
      !
      !        To add constraints, use computed ALPHA to interpolate between
      !        last feasible solution X(*) and current unconstrained (and
      !        infeasible) solution Z(*).
      !
      DO j = L + 1, nsoln
        X(j) = X(j) + alpha*(Z(j)-X(j))
      END DO
      feasbl = .FALSE.
      !
      !        Remove column JCON and shift columns JCON+1 through N to the
      !        left.  Swap column JCON into the N th position.  This achieves
      !        upper Hessenberg form for the nonactive constraints and
      !        leaves an upper Hessenberg matrix to retriangularize.
      !
      20 CONTINUE
      DO i = 1, m
        t = W(i,jcon)
        CALL SCOPY(N-jcon,W(i,jcon+1),Mdw,W(i,jcon),Mdw)
        W(i,N) = t
      END DO
      !
      !        Update permuted index vector to reflect this shift and swap.
      !
      itemp = Ipivot(jcon)
      DO i = jcon, N - 1
        Ipivot(i) = Ipivot(i+1)
      END DO
      Ipivot(N) = itemp
      !
      !        Similarly permute X(*) vector.
      !
      CALL SCOPY(N-jcon,X(jcon+1),1,X(jcon),1)
      X(N) = 0.E0
      nsoln = nsoln - 1
      niv = niv - 1
      !
      !        Retriangularize upper Hessenberg matrix after adding
      !        constraints.
      !
      i = krank + jcon - L
      DO j = jcon, nsoln
        IF ( Itype(i)==0.AND.Itype(i+1)==0 ) THEN
          !
          !              Zero IP1 to I in column J
          !
          IF ( W(i+1,j)/=0.E0 ) THEN
            CALL SROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0.E0
            CALL SROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF ( Itype(i)==1.AND.Itype(i+1)==1 ) THEN
          !
          !              Zero IP1 to I in column J
          !
          IF ( W(i+1,j)/=0.E0 ) THEN
            CALL SROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0.E0
            CALL SROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF ( Itype(i)==1.AND.Itype(i+1)==0 ) THEN
          CALL SSWAP(N+1,W(i,1),Mdw,W(i+1,1),Mdw)
          CALL SSWAP(1,Scalee(i),1,Scalee(i+1),1)
          itemp = Itype(i+1)
          Itype(i+1) = Itype(i)
          Itype(i) = itemp
          !
          !              Swapped row was formerly a pivot element, so it will
          !              be large enough to perform elimination.
          !              Zero IP1 to I in column J.
          !
          IF ( W(i+1,j)/=0.E0 ) THEN
            CALL SROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0.E0
            CALL SROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF ( Itype(i)==0.AND.Itype(i+1)==1 ) THEN
          IF ( Scalee(i)*W(i,j)**2/alsq<=(tau*eanorm)**2 ) THEN
            CALL SSWAP(N+1,W(i,1),Mdw,W(i+1,1),Mdw)
            CALL SSWAP(1,Scalee(i),1,Scalee(i+1),1)
            itemp = Itype(i+1)
            Itype(i+1) = Itype(i)
            Itype(i) = itemp
            W(i+1,j) = 0.E0
            !
            !                 Zero IP1 to I in column J
            !
          ELSEIF ( W(i+1,j)/=0.E0 ) THEN
            CALL SROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0.E0
            CALL SROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        END IF
        i = i + 1
      END DO
      !
      !        See if the remaining coefficients in the solution set are
      !        feasible.  They should be because of the way ALPHA was
      !        determined.  If any are infeasible, it is due to roundoff
      !        error.  Any that are non-positive will be set to zero and
      !        removed from the solution set.
      !
      DO jcon = L + 1, nsoln
        IF ( X(jcon)<=0.E0 ) GOTO 40
      END DO
      feasbl = .TRUE.
      40 IF ( .NOT.feasbl ) GOTO 20
    ELSE
      !
      !        To perform multiplier test and drop a constraint.
      !
      CALL SCOPY(nsoln,Z,1,X,1)
      IF ( nsoln<N ) X(nsoln+1:N) = 0.E0
      !
      !        Reclassify least squares equations as equalities as necessary.
      !
      i = niv + 1
      DO
        IF ( i<=me ) THEN
          IF ( Itype(i)==0 ) THEN
            i = i + 1
          ELSE
            CALL SSWAP(N+1,W(i,1),Mdw,W(me,1),Mdw)
            CALL SSWAP(1,Scalee(i),1,Scalee(me),1)
            itemp = Itype(i)
            Itype(i) = Itype(me)
            Itype(me) = itemp
            me = me - 1
          END IF
          CYCLE
        END IF
        !
        !        Form inner product vector WD(*) of dual coefficients.
        !
        DO j = nsoln + 1, N
          sm = 0.E0
          DO i = nsoln + 1, m
            sm = sm + Scalee(i)*W(i,j)*W(i,N+1)
          END DO
          Wd(j) = sm
        END DO
        EXIT
      END DO
      DO
        !
        !        Find J such that WD(J)=WMAX is maximum.  This determines
        !        that the incoming column J will reduce the residual vector
        !        and be positive.
        !
        wmax = 0.E0
        iwmax = nsoln + 1
        DO j = nsoln + 1, N
          IF ( Wd(j)>wmax ) THEN
            wmax = Wd(j)
            iwmax = j
          END IF
        END DO
        IF ( wmax<=0.E0 ) GOTO 100
        !
        !        Set dual coefficients to zero for incoming column.
        !
        Wd(iwmax) = 0.E0
        !
        !        WMAX .GT. 0.E0, so okay to move column IWMAX to solution set.
        !        Perform transformation to retriangularize, and test for near
        !        linear dependence.
        !
        !        Swap column IWMAX into NSOLN-th position to maintain upper
        !        Hessenberg form of adjacent columns, and add new column to
        !        triangular decomposition.
        !
        nsoln = nsoln + 1
        niv = niv + 1
        IF ( nsoln/=iwmax ) THEN
          CALL SSWAP(m,W(1,nsoln),1,W(1,iwmax),1)
          Wd(iwmax) = Wd(nsoln)
          Wd(nsoln) = 0.E0
          itemp = Ipivot(nsoln)
          Ipivot(nsoln) = Ipivot(iwmax)
          Ipivot(iwmax) = itemp
        END IF
        !
        !        Reduce column NSOLN so that the matrix of nonactive constraints
        !        variables is triangular.
        !
        DO j = m, niv + 1, -1
          jp = j - 1
          !
          !           When operating near the ME line, test to see if the pivot
          !           element is near zero.  If so, use the largest element above
          !           it as the pivot.  This is to maintain the sharp interface
          !           between weighted and non-weighted rows in all cases.
          !
          IF ( j==me+1 ) THEN
            imax = me
            amax = Scalee(me)*W(me,nsoln)**2
            DO jp = j - 1, niv, -1
              t = Scalee(jp)*W(jp,nsoln)**2
              IF ( t>amax ) THEN
                imax = jp
                amax = t
              END IF
            END DO
            jp = imax
          END IF
          !
          IF ( W(j,nsoln)/=0.E0 ) THEN
            CALL SROTMG(Scalee(jp),Scalee(j),W(jp,nsoln),W(j,nsoln),sparam)
            W(j,nsoln) = 0.E0
            CALL SROTM(N+1-nsoln,W(jp,nsoln+1),Mdw,W(j,nsoln+1),Mdw,sparam)
          END IF
        END DO
        !
        !        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
        !        this is nonpositive or too large.  If this was true or if the
        !        pivot term was zero, reject the column as dependent.
        !
        IF ( W(niv,nsoln)/=0.E0 ) THEN
          isol = niv
          z2 = W(isol,N+1)/W(isol,nsoln)
          Z(nsoln) = z2
          pos = z2>0.E0
          IF ( z2*eanorm>=bnorm.AND.pos )&
            pos = .NOT.(blowup*z2*eanorm>=bnorm)
          !
          !           Try to add row ME+1 as an additional equality constraint.
          !           Check size of proposed new solution component.
          !           Reject it if it is too large.
          !
        ELSEIF ( niv<=me.AND.W(me+1,nsoln)/=0.E0 ) THEN
          isol = me + 1
          IF ( pos ) THEN
            !
            !              Swap rows ME+1 and NIV, and scale factors for these rows.
            !
            CALL SSWAP(N+1,W(me+1,1),Mdw,W(niv,1),Mdw)
            CALL SSWAP(1,Scalee(me+1),1,Scalee(niv),1)
            itemp = Itype(me+1)
            Itype(me+1) = Itype(niv)
            Itype(niv) = itemp
            me = me + 1
          END IF
        ELSE
          pos = .FALSE.
        END IF
        !
        IF ( .NOT.pos ) THEN
          nsoln = nsoln - 1
          niv = niv - 1
        END IF
        IF ( pos.OR.done ) EXIT
      END DO
    END IF
  END DO
  !
  !     Else perform multiplier test and drop a constraint.  To compute
  !     final solution.  Solve system, store results in X(*).
  !
  !     Copy right hand side into TEMP vector to use overwriting method.
  !
  100  isol = 1
  IF ( nsoln>=isol ) THEN
    CALL SCOPY(niv,W(1,N+1),1,Temp,1)
    DO j = nsoln, isol, -1
      IF ( j>krank ) THEN
        i = niv - nsoln + j
      ELSE
        i = j
      END IF
      !
      IF ( j>krank.AND.j<=L ) THEN
        Z(j) = 0.E0
      ELSE
        Z(j) = Temp(i)/W(i,j)
        CALL SAXPY(i-1,-Z(j),W(1,j),1,Temp,1)
      END IF
    END DO
  END IF
  !
  !     Solve system.
  !
  CALL SCOPY(nsoln,Z,1,X,1)
  !
  !     Apply Householder transformations to X(*) if KRANK.LT.L
  !
  IF ( krank<L ) THEN
    DO i = 1, krank
      CALL H12(2,i,krank+1,L,W(i,1),Mdw,H(i),X,1,1,1)
    END DO
  END IF
  !
  !     Fill in trailing zeroes for constrained variables not in solution.
  !
  IF ( nsoln<N ) X(nsoln+1:N) = 0.E0
  !
  !     Permute solution vector to natural order.
  !
  DO i = 1, N
    j = i
    DO WHILE ( Ipivot(j)/=i )
      j = j + 1
    END DO
    !
    Ipivot(j) = Ipivot(i)
    Ipivot(i) = j
    CALL SSWAP(1,X(j),1,X(i),1)
  END DO
  !
  !     Rescale the solution using the column scaling.
  !
  DO j = 1, N
    X(j) = X(j)*D(j)
  END DO
  !
  DO i = nsoln + 1, m
    t = W(i,N+1)
    IF ( i<=me ) t = t/alamda
    t = (Scalee(i)*t)*t
    Rnorm = Rnorm + t
  END DO
  !
  Rnorm = SQRT(Rnorm)
END SUBROUTINE WNLSM
