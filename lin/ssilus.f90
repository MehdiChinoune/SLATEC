!*==SSILUS.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SSILUS
SUBROUTINE SSILUS(N,Nelt,Ia,Ja,A,Isym,Nl,Il,Jl,L,Dinv,Nu,Iu,Ju,U,Nrow,&
    Ncol)
  IMPLICIT NONE
  !*--SSILUS6
  !***BEGIN PROLOGUE  SSILUS
  !***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
  !            Routine to generate the incomplete LDU decomposition of a
  !            matrix.  The unit lower triangular factor L is stored by
  !            rows and the unit upper triangular factor U is stored by
  !            columns.  The inverse of the diagonal matrix D is stored.
  !            No fill in is allowed.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2E
  !***TYPE      SINGLE PRECISION (SSILUS-S, DSILUS-D)
  !***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION,
  !             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
  !***AUTHOR  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  ! *Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
  !     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)
  !     INTEGER NROW(N), NCOL(N)
  !     REAL    A(NELT), L(NL), DINV(N), U(NU)
  !
  !     CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,
  !    $    DINV, NU, IU, JU, U, NROW, NCOL )
  !
  ! *Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! NELT   :IN       Integer.
  !         Number of elements in arrays IA, JA, and A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Real A(NELT).
  !         These arrays should hold the matrix A in the SLAP Column
  !         format.  See "Description", below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! NL     :OUT      Integer.
  !         Number of non-zeros in the L array.
  ! IL     :OUT      Integer IL(NL).
  ! JL     :OUT      Integer JL(NL).
  ! L      :OUT      Real     L(NL).
  !         IL, JL, L  contain the unit lower triangular factor of  the
  !         incomplete decomposition  of some  matrix stored  in   SLAP
  !         Row format.     The   Diagonal  of ones  *IS*  stored.  See
  !         "DESCRIPTION", below for more details about the SLAP format.
  ! NU     :OUT      Integer.
  !         Number of non-zeros in the U array.
  ! IU     :OUT      Integer IU(NU).
  ! JU     :OUT      Integer JU(NU).
  ! U      :OUT      Real     U(NU).
  !         IU, JU, U contain   the unit upper triangular factor of the
  !         incomplete  decomposition    of some matrix  stored in SLAP
  !         Column  format.   The Diagonal of ones   *IS*  stored.  See
  !         "Description", below  for  more  details  about  the   SLAP
  !         format.
  ! NROW   :WORK     Integer NROW(N).
  !         NROW(I) is the number of non-zero elements in the I-th row
  !         of L.
  ! NCOL   :WORK     Integer NCOL(N).
  !         NCOL(I) is the number of non-zero elements in the I-th
  !         column of U.
  !
  ! *Description
  !       IL, JL, L should contain the unit  lower triangular factor of
  !       the incomplete decomposition of the A matrix  stored in SLAP
  !       Row format.  IU, JU, U should contain  the unit upper factor
  !       of the  incomplete decomposition of  the A matrix  stored in
  !       SLAP Column format This ILU factorization can be computed by
  !       the SSILUS routine. The diagonals (which are all one's) are
  !       stored.
  !
  !       =================== S L A P Column format ==================
  !
  !       This routine  requires that  the matrix A  be stored in  the
  !       SLAP Column format.  In this format the non-zeros are stored
  !       counting down columns (except for  the diagonal entry, which
  !       must appear first in each  "column")  and are stored  in the
  !       real array A.  In other words, for each column in the matrix
  !       put the diagonal entry in A.  Then put in the other non-zero
  !       elements going down   the  column (except  the diagonal)  in
  !       order.  The IA array holds the row  index for each non-zero.
  !       The JA array holds the offsets into the IA, A arrays for the
  !       beginning of   each    column.    That  is,    IA(JA(ICOL)),
  !       A(JA(ICOL)) points to the beginning of the ICOL-th column in
  !       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
  !       end  of   the ICOL-th  column.  Note   that  we  always have
  !       JA(N+1) = NELT+1, where  N  is the number of columns in  the
  !       matrix and  NELT   is the number of non-zeros in the matrix.
  !
  !       Here is an example of the  SLAP Column  storage format for a
  !       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
  !       column):
  !
  !           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
  !                              1  2  3    4  5    6  7    8    9 10 11
  !       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
  !       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !       ==================== S L A P Row format ====================
  !
  !       This routine requires  that the matrix A  be  stored  in the
  !       SLAP  Row format.   In this format  the non-zeros are stored
  !       counting across  rows (except for the diagonal  entry, which
  !       must appear first in each "row") and  are stored in the real
  !       array A.  In other words, for each row in the matrix put the
  !       diagonal entry in  A.   Then   put  in the   other  non-zero
  !       elements   going  across the  row (except   the diagonal) in
  !       order.   The  JA array  holds   the column   index for  each
  !       non-zero.   The IA  array holds the  offsets into  the JA, A
  !       arrays  for   the   beginning  of   each  row.   That    is,
  !       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
  !       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
  !       points to the  end of the  IROW-th row.  Note that we always
  !       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
  !       the matrix  and NELT  is the  number   of  non-zeros in  the
  !       matrix.
  !
  !       Here is an example of the SLAP Row storage format for a  5x5
  !       Matrix (in the A and JA arrays '|' denotes the end of a row):
  !
  !           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
  !                              1  2  3    4  5    6  7    8    9 10 11
  !       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
  !       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
  !       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  !***SEE ALSO  SILUR
  !***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
  !                  Johns Hopkins University Press, Baltimore, Maryland,
  !                  1983.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   920929  Corrected format of reference.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  !***END PROLOGUE  SSILUS
  !     .. Scalar Arguments ..
  INTEGER Isym, N, Nelt, Nl, Nu
  !     .. Array Arguments ..
  REAL A(Nelt), Dinv(N), L(Nl), U(Nu)
  INTEGER Ia(Nelt), Il(Nl), Iu(Nu), Ja(Nelt), Jl(Nl), Ju(Nu), Ncol(N)&
    , Nrow(N)
  !     .. Local Scalars ..
  REAL temp
  INTEGER i, ibgn, icol, iend, indx, indx1, indx2, indxc1, indxc2, &
    indxr1, indxr2, irow, itemp, j, jbgn, jend, jtemp, k, &
    kc, kr
  !***FIRST EXECUTABLE STATEMENT  SSILUS
  !
  !         Count number of elements in each row of the lower triangle.
  !
  DO i = 1, N
    Nrow(i) = 0
    Ncol(i) = 0
  ENDDO
  !VD$R NOCONCUR
  !VD$R NOVECTOR
  DO icol = 1, N
    jbgn = Ja(icol) + 1
    jend = Ja(icol+1) - 1
    IF ( jbgn<=jend ) THEN
      DO j = jbgn, jend
        IF ( Ia(j)<icol ) THEN
          Ncol(icol) = Ncol(icol) + 1
        ELSE
          Nrow(Ia(j)) = Nrow(Ia(j)) + 1
          IF ( Isym/=0 ) Ncol(Ia(j)) = Ncol(Ia(j)) + 1
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  Ju(1) = 1
  Il(1) = 1
  DO icol = 1, N
    Il(icol+1) = Il(icol) + Nrow(icol)
    Ju(icol+1) = Ju(icol) + Ncol(icol)
    Nrow(icol) = Il(icol)
    Ncol(icol) = Ju(icol)
  ENDDO
  !
  !         Copy the matrix A into the L and U structures.
  DO icol = 1, N
    Dinv(icol) = A(Ja(icol))
    jbgn = Ja(icol) + 1
    jend = Ja(icol+1) - 1
    IF ( jbgn<=jend ) THEN
      DO j = jbgn, jend
        irow = Ia(j)
        IF ( irow<icol ) THEN
          !         Part of the upper triangle.
          Iu(Ncol(icol)) = irow
          U(Ncol(icol)) = A(j)
          Ncol(icol) = Ncol(icol) + 1
        ELSE
          !         Part of the lower triangle (stored by row).
          Jl(Nrow(irow)) = icol
          L(Nrow(irow)) = A(j)
          Nrow(irow) = Nrow(irow) + 1
          IF ( Isym/=0 ) THEN
            !         Symmetric...Copy lower triangle into upper triangle as well.
            Iu(Ncol(irow)) = icol
            U(Ncol(irow)) = A(j)
            Ncol(irow) = Ncol(irow) + 1
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
  !         Sort the rows of L and the columns of U.
  DO k = 2, N
    jbgn = Ju(k)
    jend = Ju(k+1) - 1
    IF ( jbgn<jend ) THEN
      DO j = jbgn, jend - 1
        DO i = j + 1, jend
          IF ( Iu(j)>Iu(i) ) THEN
            itemp = Iu(j)
            Iu(j) = Iu(i)
            Iu(i) = itemp
            temp = U(j)
            U(j) = U(i)
            U(i) = temp
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    ibgn = Il(k)
    iend = Il(k+1) - 1
    IF ( ibgn<iend ) THEN
      DO i = ibgn, iend - 1
        DO j = i + 1, iend
          IF ( Jl(i)>Jl(j) ) THEN
            jtemp = Ju(i)
            Ju(i) = Ju(j)
            Ju(j) = jtemp
            temp = L(i)
            L(i) = L(j)
            L(j) = temp
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !
  !         Perform the incomplete LDU decomposition.
  DO i = 2, N
    !
    !           I-th row of L
    indx1 = Il(i)
    indx2 = Il(i+1) - 1
    IF ( indx1<=indx2 ) THEN
      DO indx = indx1, indx2
        IF ( indx/=indx1 ) THEN
          indxr1 = indx1
          indxr2 = indx - 1
          indxc1 = Ju(Jl(indx))
          indxc2 = Ju(Jl(indx)+1) - 1
          IF ( indxc1<=indxc2 ) THEN
            kr = Jl(indxr1)
            DO
              kc = Iu(indxc1)
              IF ( kr>kc ) THEN
                indxc1 = indxc1 + 1
                IF ( indxc1<=indxc2 ) CYCLE
              ELSEIF ( kr<kc ) THEN
                indxr1 = indxr1 + 1
                IF ( indxr1<=indxr2 ) THEN
                  kr = Jl(indxr1)
                  CYCLE
                ENDIF
              ELSEIF ( kr==kc ) THEN
                L(indx) = L(indx) - L(indxr1)*Dinv(kc)*U(indxc1)
                indxr1 = indxr1 + 1
                indxc1 = indxc1 + 1
                IF ( indxr1<=indxr2.AND.indxc1<=indxc2 ) THEN
                  kr = Jl(indxr1)
                  CYCLE
                ENDIF
              ENDIF
              EXIT
            ENDDO
          ENDIF
        ENDIF
        L(indx) = L(indx)/Dinv(Jl(indx))
      ENDDO
    ENDIF
    !
    !         I-th column of U
    indx1 = Ju(i)
    indx2 = Ju(i+1) - 1
    IF ( indx1<=indx2 ) THEN
      DO indx = indx1, indx2
        IF ( indx/=indx1 ) THEN
          indxc1 = indx1
          indxc2 = indx - 1
          indxr1 = Il(Iu(indx))
          indxr2 = Il(Iu(indx)+1) - 1
          IF ( indxr1<=indxr2 ) THEN
            kr = Jl(indxr1)
            DO
              kc = Iu(indxc1)
              IF ( kr>kc ) THEN
                indxc1 = indxc1 + 1
                IF ( indxc1<=indxc2 ) CYCLE
              ELSEIF ( kr<kc ) THEN
                indxr1 = indxr1 + 1
                IF ( indxr1<=indxr2 ) THEN
                  kr = Jl(indxr1)
                  CYCLE
                ENDIF
              ELSEIF ( kr==kc ) THEN
                U(indx) = U(indx) - L(indxr1)*Dinv(kc)*U(indxc1)
                indxr1 = indxr1 + 1
                indxc1 = indxc1 + 1
                IF ( indxr1<=indxr2.AND.indxc1<=indxc2 ) THEN
                  kr = Jl(indxr1)
                  CYCLE
                ENDIF
              ENDIF
              EXIT
            ENDDO
          ENDIF
        ENDIF
        U(indx) = U(indx)/Dinv(Iu(indx))
      ENDDO
    ENDIF
    !
    !         I-th diagonal element
    indxr1 = Il(i)
    indxr2 = Il(i+1) - 1
    IF ( indxr1<=indxr2 ) THEN
      indxc1 = Ju(i)
      indxc2 = Ju(i+1) - 1
      IF ( indxc1<=indxc2 ) THEN
        kr = Jl(indxr1)
        DO
          kc = Iu(indxc1)
          IF ( kr>kc ) THEN
            indxc1 = indxc1 + 1
            IF ( indxc1<=indxc2 ) CYCLE
          ELSEIF ( kr<kc ) THEN
            indxr1 = indxr1 + 1
            IF ( indxr1<=indxr2 ) THEN
              kr = Jl(indxr1)
              CYCLE
            ENDIF
          ELSEIF ( kr==kc ) THEN
            Dinv(i) = Dinv(i) - L(indxr1)*Dinv(kc)*U(indxc1)
            indxr1 = indxr1 + 1
            indxc1 = indxc1 + 1
            IF ( indxr1<=indxr2.AND.indxc1<=indxc2 ) THEN
              kr = Jl(indxr1)
              CYCLE
            ENDIF
          ENDIF
          EXIT
        ENDDO
      ENDIF
    ENDIF
    !
  ENDDO
  !
  !         Replace diagonal elements by their inverses.
  !VD$ VECTOR
  DO i = 1, N
    Dinv(i) = 1.0E0/Dinv(i)
  ENDDO
  !
  !------------- LAST LINE OF SSILUS FOLLOWS ----------------------------
END SUBROUTINE SSILUS
