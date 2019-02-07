!*==SBHIN.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SBHIN
SUBROUTINE SBHIN(N,Nelt,Ia,Ja,A,Isym,Soln,Rhs,Iunit,Job)
  IMPLICIT NONE
  !*--SBHIN5
  !***BEGIN PROLOGUE  SBHIN
  !***PURPOSE  Read a Sparse Linear System in the Boeing/Harwell Format.
  !            The matrix is read in and if the right hand side is also
  !            present in the input file then it too is read in.  The
  !            matrix is then modified to be in the SLAP Column format.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  N1
  !***TYPE      SINGLE PRECISION (SBHIN-S, DBHIN-D)
  !***KEYWORDS  LINEAR SYSTEM, MATRIX READ, SLAP SPARSE
  !***AUTHOR  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  ! *Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
  !     REAL    A(NELT), SOLN(N), RHS(N)
  !
  !     CALL SBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
  !
  ! *Arguments:
  ! N      :OUT      Integer
  !         Order of the Matrix.
  ! NELT   :INOUT    Integer.
  !         On input NELT is the maximum number of non-zeros that
  !         can be stored in the IA, JA, A arrays.
  !         On output NELT is the number of non-zeros stored in A.
  ! IA     :OUT      Integer IA(NELT).
  ! JA     :OUT      Integer JA(NELT).
  ! A      :OUT      Real A(NELT).
  !         On output these arrays hold the matrix A in the SLAP
  !         Triad format.  See "Description", below.
  ! ISYM   :OUT      Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! SOLN   :OUT      Real SOLN(N).
  !         The solution to the linear system, if present.  This array
  !         is accessed if and only if JOB is set to read it in, see
  !         below.  If the user requests that SOLN be read in, but it is
  !         not in the file, then it is simply zeroed out.
  ! RHS    :OUT      Real RHS(N).
  !         The right hand side vector.  This array is accessed if and
  !         only if JOB is set to read it in, see below.
  !         If the user requests that RHS be read in, but it is not in
  !         the file, then it is simply zeroed out.
  ! IUNIT  :IN       Integer.
  !         Fortran logical I/O device unit number to read the matrix
  !         from.  This unit must be connected in a system dependent
  !         fashion to a file, or you will get a nasty message
  !         from the Fortran I/O libraries.
  ! JOB    :INOUT    Integer.
  !         Flag indicating what I/O operations to perform.
  !         On input JOB indicates what Input operations to try to
  !         perform.
  !         JOB = 0 => Read only the matrix.
  !         JOB = 1 => Read matrix and RHS (if present).
  !         JOB = 2 => Read matrix and SOLN (if present).
  !         JOB = 3 => Read matrix, RHS and SOLN (if present).
  !         On output JOB indicates what operations were actually
  !         performed.
  !         JOB = -3 => Unable to parse matrix "CODE" from input file
  !                     to determine if only the lower triangle of matrix
  !                     is stored.
  !         JOB = -2 => Number of non-zeros (NELT) too large.
  !         JOB = -1 => System size (N) too large.
  !         JOB =  0 => Read in only the matrix.
  !         JOB =  1 => Read in the matrix and RHS.
  !         JOB =  2 => Read in the matrix and SOLN.
  !         JOB =  3 => Read in the matrix, RHS and SOLN.
  !         JOB = 10 => Read in only the matrix *STRUCTURE*, but no
  !                     non-zero entries.  Hence, A(*) is not referenced
  !                     and has the return values the same as the input.
  !         JOB = 11 => Read in the matrix *STRUCTURE* and RHS.
  !         JOB = 12 => Read in the matrix *STRUCTURE* and SOLN.
  !         JOB = 13 => Read in the matrix *STRUCTURE*, RHS and SOLN.
  !
  ! *Description:
  !       The format for the input is as follows.  The first line contains
  !       a title to identify the data file.  On the second line (5I4) are
  !       counters: NLINE, NPLS, NRILS, NNVLS, NRHSLS.
  !        NLINE  Number of data lines (after the header) in the file.
  !        NPLS   Number of lines for the Column Pointer data in the file.
  !        NRILS  Number of lines for the Row indices in the file.
  !        NNVLS  Number of lines for the Matrix elements in the file.
  !        NRHSLS Number of lines for the RHS in the file.
  !       The third line (A3,11X,4I4) contains a symmetry code and some
  !       additional counters: CODE, NROW, NCOL, NIND, NELE.
  !       On the fourth line (2A16,2A20) are formats to be used to read
  !       the following data: PNTFNT, RINFMT, NVLFMT, RHSFMT.
  !       Following that are the blocks of data in the order indicated.
  !
  !       =================== S L A P Triad format ===================
  !       This routine requires that the  matrix A be   stored in  the
  !       SLAP  Triad format.  In  this format only the non-zeros  are
  !       stored.  They may appear in  *ANY* order.  The user supplies
  !       three arrays of  length NELT, where  NELT is  the number  of
  !       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
  !       each non-zero the user puts the row and column index of that
  !       matrix element  in the IA and  JA arrays.  The  value of the
  !       non-zero  matrix  element is  placed   in  the corresponding
  !       location of the A array.   This is  an  extremely  easy data
  !       structure to generate.  On  the  other hand it   is  not too
  !       efficient on vector computers for  the iterative solution of
  !       linear systems.  Hence,   SLAP changes   this  input    data
  !       structure to the SLAP Column format  for  the iteration (but
  !       does not change it back).
  !
  !       Here is an example of the  SLAP Triad   storage format for a
  !       5x5 Matrix.  Recall that the entries may appear in any order.
  !
  !           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
  !                              1  2  3  4  5  6  7  8  9 10 11
  !       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
  !       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
  !       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
  !       | 0  0  0 44  0|
  !       |51  0 53  0 55|
  !
  ! *Portability:
  !         You must make sure that IUNIT is a valid Fortran logical
  !         I/O device unit number and that the unit number has been
  !         associated with a file or the console.  This is a system
  !         dependent function.
  !
  ! *Implementation note:
  !         SOLN is not read by this version.  It will simply be
  !         zeroed out if JOB = 2 or 3 and the returned value of
  !         JOB will indicate SOLN has not been read.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   881107  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   911122  Added loop to zero out RHS if user wants to read RHS, but
  !           it's not in the input file. (MKS)
  !   911125  Minor improvements to prologue.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   921007  Corrected description of input format.  (FNF)
  !   921208  Added Implementation Note and code to zero out SOLN.  (FNF)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  !***END PROLOGUE  SBHIN
  !     .. Scalar Arguments ..
  INTEGER Isym, Iunit, Job, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt), Rhs(N), Soln(N)
  INTEGER Ia(Nelt), Ja(Nelt)
  !     .. Local Scalars ..
  REAL temp
  INTEGER i, ibgn, icol, iend, itemp, j, jobret, ncol, nele, nind, &
    nline, nnvls, npls, nrhsls, nrils, nrow
  CHARACTER code*3, pntfmt*16, rinfmt*16, nvlfmt*20, rhsfmt*20, &
    title*80
  !     .. Intrinsic Functions ..
  INTRINSIC MOD
  !***FIRST EXECUTABLE STATEMENT  SBHIN
  !
  !         Read Matrices In BOEING-HARWELL format.
  !
  ! TITLE  Header line to identify data file.
  ! NLINE  Number of data lines (after the header) in the file.
  ! NPLS   Number of lines for the Column Pointer data in the file.
  ! NRILS  Number of lines for the Row indices in the data file.
  ! NNVLS  Number of lines for the Matrix elements in the data file.
  ! NRHSLS Number of lines for the RHS in the data file.
  ! ---- Only those variables needed by SLAP are referenced. ----
  !
  READ (Iunit,99001) title
  99001 FORMAT (A80)
  READ (Iunit,99002) nline, npls, nrils, nnvls, nrhsls
  99002 FORMAT (5I14)
  READ (Iunit,99003) code, nrow, ncol, nind, nele
  99003 FORMAT (A3,11X,4I14)
  READ (Iunit,99004) pntfmt, rinfmt, nvlfmt, rhsfmt
  99004 FORMAT (2A16,2A20)
  !
  IF ( nrow>N ) THEN
    N = nrow
    jobret = -1
    GOTO 100
  ENDIF
  IF ( nind>Nelt ) THEN
    Nelt = nind
    jobret = -2
    GOTO 100
  ENDIF
  !
  !         Set the parameters.
  !
  N = nrow
  Nelt = nind
  IF ( code=='RUA' ) THEN
    Isym = 0
  ELSEIF ( code=='RSA' ) THEN
    Isym = 1
  ELSE
    jobret = -3
    GOTO 100
  ENDIF
  READ (Iunit,pntfmt) (Ja(i),i=1,N+1)
  READ (Iunit,rinfmt) (Ia(i),i=1,Nelt)
  jobret = 10
  IF ( nnvls>0 ) THEN
    READ (Iunit,nvlfmt) (A(i),i=1,Nelt)
    jobret = 0
  ENDIF
  IF ( MOD(Job,2)==1 ) THEN
    !
    !         User requests that the RHS be read in.  If it is in the input
    !         file, read it in; otherwise just zero it out.
    !
    IF ( nrhsls>0 ) THEN
      READ (5,rhsfmt) (Rhs(i),i=1,N)
      jobret = jobret + 1
    ELSE
      DO i = 1, N
        Rhs(i) = 0
      ENDDO
    ENDIF
  ENDIF
  IF ( (Job==2).OR.(Job==3) ) THEN
    !
    !         User requests that the SOLN be read in.
    !         Just zero out the array.
    !
    DO i = 1, N
      Soln(i) = 0
    ENDDO
  ENDIF
  !
  !         Now loop through the IA array making sure that the diagonal
  !         matrix element appears first in the column.  Then sort the
  !         rest of the column in ascending order.
  !
  !VD$R NOCONCUR
  !VD$R NOVECTOR
  DO icol = 1, N
    ibgn = Ja(icol)
    iend = Ja(icol+1) - 1
    DO i = ibgn, iend
      IF ( Ia(i)==icol ) THEN
        !
        !              Swap the diagonal element with the first element in the
        !              column.
        !
        itemp = Ia(i)
        Ia(i) = Ia(ibgn)
        Ia(ibgn) = itemp
        temp = A(i)
        A(i) = A(ibgn)
        A(ibgn) = temp
        EXIT
      ENDIF
    ENDDO
    ibgn = ibgn + 1
    IF ( ibgn<iend ) THEN
      DO i = ibgn, iend
        DO j = i + 1, iend
          IF ( Ia(i)>Ia(j) ) THEN
            itemp = Ia(i)
            Ia(i) = Ia(j)
            Ia(j) = itemp
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !
  !         Set return flag.
  100  Job = jobret
  RETURN
  !------------- LAST LINE OF SBHIN FOLLOWS ------------------------------
END SUBROUTINE SBHIN
