!*==SSDSCL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SSDSCL
SUBROUTINE SSDSCL(N,Nelt,Ia,Ja,A,Isym,X,B,Dinv,Job,Itol)
  IMPLICIT NONE
  !*--SSDSCL5
  !***BEGIN PROLOGUE  SSDSCL
  !***PURPOSE  Diagonal Scaling of system Ax = b.
  !            This routine scales (and unscales) the system  Ax = b
  !            by symmetric diagonal scaling.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2E
  !***TYPE      SINGLE PRECISION (SSDSCL-S, DSDSCL-D)
  !***KEYWORDS  DIAGONAL, SLAP SPARSE
  !***AUTHOR  Greenbaum, Anne, (Courant Institute)
  !           Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  !    This routine scales (and unscales) the system Ax = b by symmetric
  !    diagonal scaling.  The new system is:
  !         -1/2  -1/2  1/2      -1/2
  !        D    AD    (D   x) = D    b
  !    when scaling is selected with the JOB parameter.  When unscaling
  !    is selected this process is reversed.  The true solution is also
  !    scaled or unscaled if ITOL is set appropriately, see below.
  !
  ! *Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL
  !     REAL    A(NELT), X(N), B(N), DINV(N)
  !
  !     CALL SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL )
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
  !         If ISYM=1, the matrix is symmetric, and only the upper
  !         or lower triangle of the matrix is stored.
  ! X      :INOUT    Real X(N).
  !         Initial guess that will be later used in the iterative
  !         solution.
  !         of the scaled system.
  ! B      :INOUT    Real B(N).
  !         Right hand side vector.
  ! DINV   :INOUT    Real DINV(N).
  !         Upon return this array holds 1./DIAG(A).
  !         This is an input if JOB = 0.
  ! JOB    :IN       Integer.
  !         Flag indicating whether to scale or not.
  !         JOB non-zero means do scaling.
  !         JOB = 0 means do unscaling.
  ! ITOL   :IN       Integer.
  !         Flag indicating what type of error estimation to do in the
  !         iterative method.  When ITOL = 11 the exact solution from
  !         common block SSLBLK will be used.  When the system is scaled
  !         then the true solution must also be scaled.  If ITOL is not
  !         11 then this vector is not referenced.
  !
  ! *Common Blocks:
  ! SOLN    :INOUT   Real SOLN(N).  COMMON BLOCK /SSLBLK/
  !         The true solution, SOLN, is scaled (or unscaled) if ITOL is
  !         set to 11, see above.
  !
  ! *Description
  !       =================== S L A P Column format ==================
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
  !       With the SLAP  format  all  of  the   "inner  loops" of this
  !       routine should vectorize  on  machines with hardware support
  !       for vector   gather/scatter  operations.  Your compiler  may
  !       require a compiler directive to  convince it that  there are
  !       no  implicit  vector  dependencies.  Compiler directives for
  !       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
  !       supplied with the standard SLAP distribution.
  !
  !
  ! *Cautions:
  !       This routine assumes that the diagonal of A is all  non-zero
  !       and that the operation DINV = 1.0/DIAG(A)  will  not  under-
  !       flow or overflow. This is done so that the loop  vectorizes.
  !       Matrices  with zero or near zero or very  large entries will
  !       have numerical difficulties  and  must  be fixed before this
  !       routine is called.
  !
  !***SEE ALSO  SSDCG
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    SSLBLK
  !***REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
  !   920407  COMMON BLOCK renamed SSLBLK.  (WRB)
  !   920511  Added complete declaration section.  (WRB)
  !   921113  Corrected C***CATEGORY line.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  !***END PROLOGUE  SSDSCL
  !     .. Scalar Arguments ..
  INTEGER Isym, Itol, Job, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt), B(N), Dinv(N), X(N)
  INTEGER Ia(Nelt), Ja(Nelt)
  !     .. Arrays in Common ..
  REAL SOLn(1)
  !     .. Local Scalars ..
  REAL di
  INTEGER icol, j, jbgn, jend
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !     .. Common blocks ..
  COMMON /SSLBLK/ SOLn
  !***FIRST EXECUTABLE STATEMENT  SSDSCL
  !
  !         SCALING...
  !
  IF ( Job/=0 ) THEN
    DO icol = 1, N
      Dinv(icol) = 1.0E0/SQRT(A(Ja(icol)))
    ENDDO
  ELSE
    !
    !         UNSCALING...
    !
    DO icol = 1, N
      Dinv(icol) = 1.0E0/Dinv(icol)
    ENDDO
  ENDIF
  !
  DO icol = 1, N
    jbgn = Ja(icol)
    jend = Ja(icol+1) - 1
    di = Dinv(icol)
    DO j = jbgn, jend
      A(j) = Dinv(Ia(j))*A(j)*di
    ENDDO
  ENDDO
  !
  DO icol = 1, N
    B(icol) = B(icol)*Dinv(icol)
    X(icol) = X(icol)/Dinv(icol)
  ENDDO
  !
  !         Check to see if we need to scale the "true solution" as well.
  !
  IF ( Itol==11 ) THEN
    DO icol = 1, N
      SOLn(icol) = SOLn(icol)/Dinv(icol)
    ENDDO
  ENDIF
  !
  !------------- LAST LINE OF SSDSCL FOLLOWS ----------------------------
END SUBROUTINE SSDSCL
