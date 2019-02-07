!*==DWNLIT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DWNLIT
SUBROUTINE DWNLIT(W,Mdw,M,N,L,Ipivot,Itype,H,Scale,Rnorm,Idope,Dope,Done)
  IMPLICIT NONE
  !*--DWNLIT5
  !***BEGIN PROLOGUE  DWNLIT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DWNNLS
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
  !
  !     This is a companion subprogram to DWNNLS( ).
  !     The documentation for DWNNLS( ) has complete usage instructions.
  !
  !     Note  The M by (N+1) matrix W(, ) contains the rt. hand side
  !           B as the (N+1)st col.
  !
  !     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
  !     col interchanges.
  !
  !***SEE ALSO  DWNNLS
  !***ROUTINES CALLED  DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1,
  !                    DWNLT2, DWNLT3, IDAMAX
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and revised.  (WRB & RWC)
  !   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900604  DP version created from SP version. .  (RWC)
  !***END PROLOGUE  DWNLIT
  INTEGER Idope(*), Ipivot(*), Itype(*), L, M, Mdw, N
  REAL(8) :: Dope(*), H(*), Rnorm, Scale(*), W(Mdw,*)
  LOGICAL Done
  !
  EXTERNAL DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1, DWNLT2, &
    DWNLT3, IDAMAX
  INTEGER IDAMAX
  LOGICAL DWNLT2
  !
  REAL(8) :: alsq, amax, eanorm, factor, hbar, rn, sparam(5), &
    t, tau
  INTEGER i, i1, imax, ir, j, j1, jj, jp, krank, l1, lb, lend, &
    me, mend, niv, nsoln
  LOGICAL indep, recalc
  !
  !***FIRST EXECUTABLE STATEMENT  DWNLIT
  me = Idope(1)
  nsoln = Idope(2)
  l1 = Idope(3)
  !
  alsq = Dope(1)
  eanorm = Dope(2)
  tau = Dope(3)
  !
  lb = MIN(M-1,L)
  recalc = .TRUE.
  Rnorm = 0.D0
  krank = 0
  !
  !     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
  !     included in the test for column independence.
  !
  factor = 1.D0
  lend = L
  DO i = 1, lb
    !
    !        Set IR to point to the I-th row.
    !
    ir = i
    mend = M
    CALL DWNLT1(i,lend,M,ir,Mdw,recalc,imax,hbar,H,Scale,W)
    !
    !        Update column SS and find pivot column.
    !
    CALL DWNLT3(i,imax,M,Mdw,Ipivot,H,W)
    DO
      !
      !        Perform column interchange.
      !        Test independence of incoming column.
      !
      IF ( DWNLT2(me,mend,ir,factor,tau,Scale,W(1,i)) ) THEN
        !
        !           Eliminate I-th column below diagonal using modified Givens
        !           transformations applied to (A B).
        !
        !           When operating near the ME line, use the largest element
        !           above it as the pivot.
        !
        DO j = M, i + 1, -1
          jp = j - 1
          IF ( j==me+1 ) THEN
            imax = me
            amax = Scale(me)*W(me,i)**2
            DO jp = j - 1, i, -1
              t = Scale(jp)*W(jp,i)**2
              IF ( t>amax ) THEN
                imax = jp
                amax = t
              ENDIF
            ENDDO
            jp = imax
          ENDIF
          !
          IF ( W(j,i)/=0.D0 ) THEN
            CALL DROTMG(Scale(jp),Scale(j),W(jp,i),W(j,i),sparam)
            W(j,i) = 0.D0
            CALL DROTM(N+1-i,W(jp,i+1),Mdw,W(j,i+1),Mdw,sparam)
          ENDIF
        ENDDO
      ELSEIF ( lend>i ) THEN
        !
        !           Column I is dependent.  Swap with column LEND.
        !           Perform column interchange,
        !           and find column in remaining set with largest SS.
        !
        CALL DWNLT3(i,lend,M,Mdw,Ipivot,H,W)
        lend = lend - 1
        imax = IDAMAX(lend-i+1,H(i),1) + i - 1
        hbar = H(imax)
        CYCLE
      ELSE
        krank = i - 1
        GOTO 100
      ENDIF
      EXIT
    ENDDO
  ENDDO
  krank = l1
  !
  100 CONTINUE
  IF ( krank<me ) THEN
    factor = alsq
    DO i = krank + 1, me
      CALL DCOPY(L,0.D0,0,W(i,1),Mdw)
    ENDDO
    !
    !        Determine the rank of the remaining equality constraint
    !        equations by eliminating within the block of constrained
    !        variables.  Remove any redundant constraints.
    !
    recalc = .TRUE.
    lb = MIN(L+me-krank,N)
    DO i = L + 1, lb
      ir = krank + i - L
      lend = N
      mend = me
      CALL DWNLT1(i,lend,me,ir,Mdw,recalc,imax,hbar,H,Scale,W)
      !
      !           Update col ss and find pivot col
      !
      CALL DWNLT3(i,imax,M,Mdw,Ipivot,H,W)
      !
      !           Perform column interchange
      !           Eliminate elements in the I-th col.
      !
      DO j = me, ir + 1, -1
        IF ( W(j,i)/=0.D0 ) THEN
          CALL DROTMG(Scale(j-1),Scale(j),W(j-1,i),W(j,i),sparam)
          W(j,i) = 0.D0
          CALL DROTM(N+1-i,W(j-1,i+1),Mdw,W(j,i+1),Mdw,sparam)
        ENDIF
      ENDDO
      !
      !           I=column being eliminated.
      !           Test independence of incoming column.
      !           Remove any redundant or dependent equality constraints.
      !
      IF ( .NOT.DWNLT2(me,mend,ir,factor,tau,Scale,W(1,i)) ) THEN
        jj = ir
        DO ir = jj, me
          CALL DCOPY(N,0.D0,0,W(ir,1),Mdw)
          Rnorm = Rnorm + (Scale(ir)*W(ir,N+1)/alsq)*W(ir,N+1)
          W(ir,N+1) = 0.D0
          Scale(ir) = 1.D0
          !
          !                 Reclassify the zeroed row as a least squares equation.
          !
          Itype(ir) = 1
        ENDDO
        !
        !              Reduce ME to reflect any discovered dependent equality
        !              constraints.
        !
        me = jj - 1
        EXIT
      ENDIF
    ENDDO
  ENDIF
  !
  !     Try to determine the variables KRANK+1 through L1 from the
  !     least squares equations.  Continue the triangularization with
  !     pivot element W(ME+1,I).
  !
  IF ( krank<l1 ) THEN
    recalc = .TRUE.
    !
    !        Set FACTOR=ALSQ to remove effect of heavy weight from
    !        test for column independence.
    !
    factor = alsq
    DO i = krank + 1, l1
      !
      !           Set IR to point to the ME+1-st row.
      !
      ir = me + 1
      lend = L
      mend = M
      CALL DWNLT1(i,L,M,ir,Mdw,recalc,imax,hbar,H,Scale,W)
      !
      !           Update column SS and find pivot column.
      !
      CALL DWNLT3(i,imax,M,Mdw,Ipivot,H,W)
      !
      !           Perform column interchange.
      !           Eliminate I-th column below the IR-th element.
      !
      DO j = M, ir + 1, -1
        IF ( W(j,i)/=0.D0 ) THEN
          CALL DROTMG(Scale(j-1),Scale(j),W(j-1,i),W(j,i),sparam)
          W(j,i) = 0.D0
          CALL DROTM(N+1-i,W(j-1,i+1),Mdw,W(j,i+1),Mdw,sparam)
        ENDIF
      ENDDO
      !
      !           Test if new pivot element is near zero.
      !           If so, the column is dependent.
      !           Then check row norm test to be classified as independent.
      !
      t = Scale(ir)*W(ir,i)**2
      indep = t>(tau*eanorm)**2
      IF ( indep ) THEN
        rn = 0.D0
        DO i1 = ir, M
          DO j1 = i + 1, N
            rn = MAX(rn,Scale(i1)*W(i1,j1)**2)
          ENDDO
        ENDDO
        indep = t>rn*tau**2
      ENDIF
      !
      !           If independent, swap the IR-th and KRANK+1-th rows to
      !           maintain the triangular form.  Update the rank indicator
      !           KRANK and the equality constraint pointer ME.
      !
      IF ( .NOT.indep ) EXIT
      CALL DSWAP(N+1,W(krank+1,1),Mdw,W(ir,1),Mdw)
      CALL DSWAP(1,Scale(krank+1),1,Scale(ir),1)
      !
      !           Reclassify the least square equation as an equality
      !           constraint and rescale it.
      !
      Itype(ir) = 0
      t = SQRT(Scale(krank+1))
      CALL DSCAL(N+1,t,W(krank+1,1),Mdw)
      Scale(krank+1) = alsq
      me = me + 1
      krank = krank + 1
    ENDDO
  ENDIF
  !
  !     If pseudorank is less than L, apply Householder transformation.
  !     from right.
  !
  IF ( krank<L ) THEN
    DO j = krank, 1, -1
      CALL DH12(1,j,krank+1,L,W(j,1),Mdw,H(j),W,Mdw,1,j-1)
    ENDDO
  ENDIF
  !
  niv = krank + nsoln - L
  IF ( L==N ) Done = .TRUE.
  !
  !     End of initial triangularization.
  !
  Idope(1) = me
  Idope(2) = krank
  Idope(3) = niv
END SUBROUTINE DWNLIT
