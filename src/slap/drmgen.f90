!** DRMGEN
SUBROUTINE DRMGEN(Neltmx,Factor,Ierr,N,Nelt,Isym,Ia,Ja,A,F,Dsum,Itmp,Idiag)
  IMPLICIT NONE
  !>
  !***
  !  This routine generates a "Random" symmetric or non-symmetric matrix
  !  of size N for use in the SLAP Quick Checks.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Type:**      DOUBLE PRECISION (SRMGEN-S, DRMGEN-D)
  !***
  ! **Author:**  Seager, Mark K., (LLNL)
  !             seager@llnl.gov
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-300
  !             Livermore, CA 94550
  !             (510) 423-3141
  !***
  ! **Description:**
  !
  !- Usage:
  !       INTEGER NELTMX, IERR, N, NELT, ISYM,
  !       INTEGER IA(NELTMX), JA(NELTMX), ITMP(N), IDIAG(N)
  !       DOUBLE PRECISION FACTOR, A(NELTMX), F(N), SOLN(N), DSUM(N)
  !
  !       CALL DRMGEN( NELTMX, FACTOR, IERR, N, NELT, ISYM,
  !      $     IA, JA, A, F, SOLN, DSUM, ITMP, IDIAG )
  !
  !- Arguments:
  !
  ! NELTMX :IN       Integer.
  !         Maximum number of non-zeros that can be created by this
  !         routine for storage in the IA, JA, A arrays,  see below.
  ! FACTOR :IN       Double Precision.
  !         Non-zeros in the upper triangle are set to FACTOR times
  !         the corresponding entry in the lower triangle when a non-
  !         symmetric matrix is requested (See ISYM, below).
  ! IERR   :OUT      Integer.
  !         Return error flag.
  !             IERR = 0 => everything went OK.
  !                  = 1 => Ran out of space trying to create matrix.
  !                         Set NELTMX to something larger and retry.
  ! N      :IN       Integer.
  !         Size of the linear system to generate (number of unknowns).
  ! NELT   :OUT      Integer.
  !         Number of non-zeros stored in the IA, JA, A arrays, see below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate the type of matrix to generate:
  !             ISYM = 0 => Non-Symmetric Matrix (See FACTOR, above).
  !                  = 1 => Symmetric Matrix.
  ! IA     :OUT      Integer IA(NELTMX).
  !         Stores the row indices for the non-zeros.
  ! JA     :OUT      Integer JA(NELTMX).
  !         Stores the column indices for the non-zeros.
  ! A      :OUT      Double Precision A(NELTMX).
  !         Stores the values of the non-zeros.
  ! F      :OUT      Double Precision F(N).
  !         The right hand side of the linear system.  Obtained by
  !         multiplying the matrix times SOLN, see below.
  ! SOLN   :OUT      Double Precision SOLN(N).
  !         The true solution to the linear system.  Each component is
  !         chosen at random (0.0<SOLN(I)<1.0, I=1,N)
  ! DSUM   :WORK     Double Precision DSUM(N).
  ! ITMP   :WORK     Integer ITMP(N).
  ! IDIAG  :WORK     Integer IDIAG(N).
  !
  !- Description
  !         The matrix is generated by choosing a random number of
  !         entries for each column and then chosing negative random
  !         numbers for each off diagonal.   The diagonal elements
  !         are chosen to be positive and large enough so the matrix
  !         is slightly diagonally dominant.  The lower triangle of
  !         the matrix is generated and if isym.eq.0 (all matrix elements
  !         stored) the upper triangle elements are chosen so that they
  !         are FACTOR times the corresponding lower triangular element.
  !
  !***
  ! **Routines called:**  ISMPL, RAND

  !* REVISION HISTORY  (YYMMDD)
  !   881120  DATE WRITTEN
  !   890919  Replaced DMPL with ISMPL.  (MKS)
  !   890920  Minor changes to reduce single/double differences.  (FNF)
  !   920511  Added complete declaration section.  (WRB)

  !     .. Scalar Arguments ..
  REAL(8) :: Factor
  INTEGER Ierr, Isym, N, Nelt, Neltmx
  !     .. Array Arguments ..
  REAL(8) :: A(Neltmx), Dsum(N), F(N)
  INTEGER Ia(Neltmx), Idiag(N), Itmp(N), Ja(Neltmx)
  !     .. Arrays in Common ..
  REAL(8) :: SOLn(25)
  !     .. Local Scalars ..
  REAL dummy
  INTEGER i, icol, inum, irow, iseed, k, nl
  !     .. External Functions ..
  REAL, EXTERNAL :: RAND
  !     .. External Subroutines ..
  EXTERNAL :: ISMPL
  !     .. Intrinsic Functions ..
  INTRINSIC INT
  !     .. Common blocks ..
  COMMON /DSLBLK/ SOLn
  !* FIRST EXECUTABLE STATEMENT  DRMGEN
  !
  !     Start by setting the random number generator seed.  This is done
  !     for reproducibility in debugging.
  !
  !     Remove the seed setting call for production testing.
  !
  !     Note:  Double precision version did not work properly with
  !            certain compilers with literal arguments to RAND.
  !
  dummy = 16381.0
  iseed = INT( RAND(dummy) )
  Ierr = 0
  DO i = 1, N
    Idiag(i) = 0
    Dsum(i) = -1.0D0
  ENDDO
  dummy = 0.0
  Nelt = 0
  !
  !     Set the matrix elements.
  !     Loop over the columns.
  !
  !VD$ NOCONCUR
  DO icol = 1, N
    nl = N + 1 - icol
    !
    !         To keep things sparse divide by two, three or four or ...
    !
    inum = (INT(RAND(dummy)*nl)+1)/3
    CALL ISMPL(nl,inum,Itmp)
    !
    !         Set up this column (and row, if non-symmetric structure).
    !VD$ NOVECTOR
    !VD$ NOCONCUR
    DO irow = 1, inum
      Nelt = Nelt + 1
      IF ( Nelt>Neltmx ) THEN
        Ierr = 1
        RETURN
      ENDIF
      Ia(Nelt) = N + 1 - Itmp(irow)
      Ja(Nelt) = icol
      IF ( Ia(Nelt)==icol ) THEN
        Idiag(icol) = Nelt
      ELSE
        A(Nelt) = -RAND(dummy)
        Dsum(icol) = Dsum(icol) + A(Nelt)
        IF ( Isym==0 ) THEN
          !
          !         Copy this element into upper triangle.
          !
          Nelt = Nelt + 1
          IF ( Nelt>Neltmx ) THEN
            Ierr = 1
            RETURN
          ENDIF
          Ia(Nelt) = icol
          Ja(Nelt) = Ia(Nelt-1)
          A(Nelt) = A(Nelt-1)*Factor
          Dsum(Ja(Nelt)) = Dsum(Ja(Nelt)) + A(Nelt)
        ELSE
          Dsum(Ia(Nelt)) = Dsum(Ia(Nelt)) + A(Nelt)
        ENDIF
      ENDIF
    ENDDO
    IF ( Idiag(icol)==0 ) THEN
      !
      !           Add a diagonal to the column.
      !
      Nelt = Nelt + 1
      IF ( Nelt>Neltmx ) THEN
        Ierr = 1
        RETURN
      ENDIF
      Idiag(icol) = Nelt
      A(Nelt) = 0.0D0
      Ia(Nelt) = icol
      Ja(Nelt) = icol
    ENDIF
  ENDDO
  !
  !         Clean up the diagonals.
  !
  !VD$ NODEPCHK
  !LLL. OPTION ASSERT (NOHAZARD)
  !DIR$ IVDEP
  DO i = 1, N
    A(Idiag(i)) = -1.0001D0*Dsum(i)
  ENDDO
  !
  !         Set a random solution and determine the right-hand side.
  !
  !VD$ NOVECTOR
  !VD$ NOCONCUR
  DO i = 1, N
    SOLn(i) = RAND(dummy)
    F(i) = 0.0D0
  ENDDO
  !
  !VD$ NOVECTOR
  !VD$ NOCONCUR
  DO k = 1, Nelt
    F(Ia(k)) = F(Ia(k)) + A(k)*SOLn(Ja(k))
    IF ( Isym/=0.AND.Ia(k)/=Ja(k) ) F(Ja(k)) = F(Ja(k)) + A(k)*SOLn(Ia(k))
  ENDDO
  !------------- LAST LINE OF DRMGEN FOLLOWS ----------------------------
END SUBROUTINE DRMGEN