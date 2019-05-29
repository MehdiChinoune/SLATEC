!** DSICS
SUBROUTINE DSICS(N,Nelt,Ia,Ja,A,Isym,Nel,Iel,Jel,El,D,R,Iwarn)
  !>
  !  Incompl. Cholesky Decomposition Preconditioner SLAP Set Up.
  !            Routine to generate the Incomplete Cholesky decomposition,
  !            L*D*L-trans, of a symmetric positive definite matrix, A,
  !            which is stored in SLAP Column format.  The unit lower
  !            triangular matrix L is stored by rows, and the inverse of
  !            the diagonal matrix D is stored.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2E
  !***
  ! **Type:**      DOUBLE PRECISION (SSICS-S, DSICS-D)
  !***
  ! **Keywords:**  INCOMPLETE CHOLESKY FACTORIZATION,
  !             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE
  !***
  ! **Author:**  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !
  !- Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
  !     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN
  !     DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)
  !
  !     CALL DSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
  !    $    IWARN )
  !
  !- Arguments:
  ! N      :IN       Integer.
  !         Order of the Matrix.
  ! NELT   :IN       Integer.
  !         Number of elements in arrays IA, JA, and A.
  ! IA     :INOUT    Integer IA(NELT).
  ! JA     :INOUT    Integer JA(NELT).
  ! A      :INOUT    Double Precision A(NELT).
  !         These arrays should hold the matrix A in the SLAP Column
  !         format.  See "Description", below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! NEL    :OUT      Integer.
  !         Number of non-zeros in the lower triangle of A.   Also
  !         corresponds to the length of the IEL, JEL, EL arrays.
  ! IEL    :OUT      Integer IEL(NEL).
  ! JEL    :OUT      Integer JEL(NEL).
  ! EL     :OUT      Double Precision EL(NEL).
  !         IEL, JEL, EL contain the unit lower triangular factor  of the
  !         incomplete decomposition   of the A  matrix  stored  in  SLAP
  !         Row format.   The Diagonal of   ones   *IS*   stored.     See
  !         "Description", below for more details about the SLAP Row fmt.
  ! D      :OUT      Double Precision D(N)
  !         Upon return this array holds D(I) = 1./DIAG(A).
  ! R      :WORK     Double Precision R(N).
  !         Temporary double precision workspace needed for the
  !         factorization.
  ! IWARN  :OUT      Integer.
  !         This is a warning variable and is zero if the IC factoriza-
  !         tion goes well.  It is set to the row index corresponding to
  !         the last zero pivot found.  See "Description", below.
  !
  !- Description
  !       =================== S L A P Column format ==================
  !       This routine  requires that  the matrix A  be stored in  the
  !       SLAP Column format.  In this format the non-zeros are stored
  !       counting down columns (except for  the diagonal entry, which
  !       must appear first in each  "column")  and are stored  in the
  !       double precision array A.   In other words,  for each column
  !       in the matrix put the diagonal entry in  A.  Then put in the
  !       other non-zero  elements going down  the column (except  the
  !       diagonal) in order.   The  IA array holds the  row index for
  !       each non-zero.  The JA array holds the offsets  into the IA,
  !       A arrays  for  the  beginning  of each   column.   That  is,
  !       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
  !       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
  !       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
  !       Note that we always have  JA(N+1) = NELT+1,  where N is  the
  !       number of columns in  the matrix and NELT  is the number  of
  !       non-zeros in the matrix.
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
  !       must  appear first  in each  "row")  and  are stored  in the
  !       double precision  array A.  In other words, for each row  in
  !       the matrix  put the diagonal  entry in A.   Then put in  the
  !       other  non-zero elements  going across  the row  (except the
  !       diagonal) in order.  The JA array holds the column index for
  !       each non-zero.  The IA array holds the offsets  into the JA,
  !       A  arrays  for  the   beginning  of  each  row.    That  is,
  !       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
  !       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
  !       are  the last elements  of the  IROW-th row.   Note  that we
  !       always have  IA(N+1) = NELT+1, where N is the number of rows
  !       in the matrix  and  NELT is the  number of non-zeros  in the
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
  !       With the SLAP  format some  of  the   "inner  loops" of this
  !       routine should vectorize  on  machines with hardware support
  !       for vector   gather/scatter  operations.  Your compiler  may
  !       require a compiler directive to  convince it that  there are
  !       no  implicit  vector  dependencies.  Compiler directives for
  !       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
  !       supplied with the standard SLAP distribution.
  !
  !       The IC factorization does not always exist for SPD matrices.
  !       In the event that a zero pivot is found it is set  to be 1.0
  !       and the factorization proceeds.   The integer variable IWARN
  !       is set to the last row where the Diagonal was fudged.  This
  !       eventuality hardly ever occurs in practice.
  !
  !***
  ! **See also:**  DCG, DSICCG
  !***
  ! **References:**  1. Gene Golub and Charles Van Loan, Matrix Computations,
  !                  Johns Hopkins University Press, Baltimore, Maryland,
  !                  1983.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   920929  Corrected format of reference.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  USE service, ONLY : XERMSG
  !     .. Scalar Arguments ..
  INTEGER Isym, Iwarn, N, Nel, Nelt
  !     .. Array Arguments ..
  REAL(8) :: A(Nelt), D(N), El(Nel), R(N)
  INTEGER Ia(Nelt), Iel(Nel), Ja(Nelt), Jel(Nel)
  !     .. Local Scalars ..
  REAL(8) :: eltmp
  INTEGER i, ibgn, ic, icbgn, icend, icol, iend, ir, irbgn, irend, &
    irow, irr, j, jbgn, jeltmp, jend
  CHARACTER xern1*8
  !* FIRST EXECUTABLE STATEMENT  DSICS
  !
  !         Set the lower triangle in IEL, JEL, EL
  !
  Iwarn = 0
  !
  !         All matrix elements stored in IA, JA, A.  Pick out the lower
  !         triangle (making sure that the Diagonal of EL is one) and
  !         store by rows.
  !
  Nel = 1
  Iel(1) = 1
  Jel(1) = 1
  El(1) = 1
  D(1) = A(1)
  DO irow = 2, N
    !         Put in the Diagonal.
    Nel = Nel + 1
    Iel(irow) = Nel
    Jel(Nel) = irow
    El(Nel) = 1
    D(irow) = A(Ja(irow))
    !
    !         Look in all the lower triangle columns for a matching row.
    !         Since the matrix is symmetric, we can look across the
    !         IROW-th row by looking down the IROW-th column (if it is
    !         stored ISYM=0)...
    IF ( Isym==0 ) THEN
      icbgn = Ja(irow)
      icend = Ja(irow+1) - 1
    ELSE
      icbgn = 1
      icend = irow - 1
    END IF
    DO ic = icbgn, icend
      IF ( Isym==0 ) THEN
        icol = Ia(ic)
        IF ( icol>=irow ) CYCLE
      ELSE
        icol = ic
      END IF
      jbgn = Ja(icol) + 1
      jend = Ja(icol+1) - 1
      IF ( jbgn<=jend.AND.Ia(jend)>=irow ) THEN
        DO j = jbgn, jend
          IF ( Ia(j)==irow ) THEN
            Nel = Nel + 1
            Jel(Nel) = icol
            El(Nel) = A(j)
            EXIT
          END IF
        END DO
      END IF
    END DO
  END DO
  Iel(N+1) = Nel + 1
  !
  !         Sort ROWS of lower triangle into descending order (count out
  !         along rows out from Diagonal).
  !
  DO irow = 2, N
    ibgn = Iel(irow) + 1
    iend = Iel(irow+1) - 1
    IF ( ibgn<iend ) THEN
      DO i = ibgn, iend - 1
        DO j = i + 1, iend
          IF ( Jel(i)>Jel(j) ) THEN
            jeltmp = Jel(j)
            Jel(j) = Jel(i)
            Jel(i) = jeltmp
            eltmp = El(j)
            El(j) = El(i)
            El(i) = eltmp
          END IF
        END DO
      END DO
    END IF
  END DO
  !
  !         Perform the Incomplete Cholesky decomposition by looping
  !         over the rows.
  !         Scale the first column.  Use the structure of A to pick out
  !         the rows with something in column 1.
  !
  irbgn = Ja(1) + 1
  irend = Ja(2) - 1
  DO irr = irbgn, irend
    ir = Ia(irr)
    !         Find the index into EL for EL(1,IR).
    !         Hint: it's the second entry.
    i = Iel(ir) + 1
    El(i) = El(i)/D(1)
  END DO
  !
  DO irow = 2, N
    !
    !         Update the IROW-th diagonal.
    !
    DO i = 1, irow - 1
      R(i) = 0
    END DO
    ibgn = Iel(irow) + 1
    iend = Iel(irow+1) - 1
    IF ( ibgn<=iend ) THEN
      DO i = ibgn, iend
        R(Jel(i)) = El(i)*D(Jel(i))
        D(irow) = D(irow) - El(i)*R(Jel(i))
      END DO
      !
      !         Check to see if we have a problem with the diagonal.
      !
      IF ( D(irow)<=0.0D0 ) THEN
        IF ( Iwarn==0 ) Iwarn = irow
        D(irow) = 1
      END IF
    END IF
    !
    !         Update each EL(IROW+1:N,IROW), if there are any.
    !         Use the structure of A to determine the Non-zero elements
    !         of the IROW-th column of EL.
    !
    irbgn = Ja(irow)
    irend = Ja(irow+1) - 1
    DO irr = irbgn, irend
      ir = Ia(irr)
      IF ( ir>irow ) THEN
        !         Find the index into EL for EL(IR,IROW)
        ibgn = Iel(ir) + 1
        iend = Iel(ir+1) - 1
        IF ( Jel(ibgn)<=irow ) THEN
          DO i = ibgn, iend
            IF ( Jel(i)==irow ) THEN
              icend = iend
              DO
                IF ( Jel(icend)>=irow ) THEN
                  icend = icend - 1
                  CYCLE
                END IF
                !         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
                DO ic = ibgn, icend
                  El(i) = El(i) - El(ic)*R(Jel(ic))
                END DO
                El(i) = El(i)/D(irow)
                GOTO 50
              END DO
            END IF
          END DO
          !
          !         If we get here, we have real problems...
          WRITE (xern1,'(I8)') irow
          CALL XERMSG('SLATEC','DSICS',&
            'A and EL data structure mismatch in row '//xern1,1,2)
        END IF
      END IF
      50 CONTINUE
    END DO
  END DO
  !
  !         Replace diagonals by their inverses.
  !
  DO i = 1, N
    D(i) = 1.0D0/D(i)
  END DO
  !------------- LAST LINE OF DSICS FOLLOWS ----------------------------
END SUBROUTINE DSICS
