!** STOUT
SUBROUTINE STOUT(N,Nelt,Ia,Ja,A,Isym,Soln,Rhs,Iunit,Job)
  IMPLICIT NONE
  !>
  !***
  !  Write out SLAP Triad Format Linear System.
  !            Routine to write out a SLAP Triad format matrix and right
  !            hand side and solution to the system, if known.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  N1
  !***
  ! **Type:**      SINGLE PRECISION (STOUT-S, DTOUT-D)
  !***
  ! **Keywords:**  DIAGNOSTICS, LINEAR SYSTEM, SLAP SPARSE
  !***
  ! **Author:**  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***
  ! **Description:**
  !
  !- Usage:
  !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB
  !     REAL    A(NELT), SOLN(N), RHS(N)
  !
  !     CALL STOUT( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB )
  !
  !- Arguments:
  ! N      :IN       Integer
  !         Order of the Matrix.
  ! NELT   :IN       Integer.
  !         Number of non-zeros stored in A.
  ! IA     :IN       Integer IA(NELT).
  ! JA     :IN       Integer JA(NELT).
  ! A      :IN       Real A(NELT).
  !         These arrays should hold the matrix A in the SLAP
  !         Triad format.  See "Description", below.
  ! ISYM   :IN       Integer.
  !         Flag to indicate symmetric storage format.
  !         If ISYM=0, all non-zero entries of the matrix are stored.
  !         If ISYM=1, the matrix is symmetric, and only the lower
  !         triangle of the matrix is stored.
  ! SOLN   :IN       Real SOLN(N).
  !         The solution to the linear system, if known.  This array
  !         is accessed if and only if JOB is set to print it out,
  !         see below.
  ! RHS    :IN       Real RHS(N).
  !         The right hand side vector.  This array is accessed if and
  !         only if JOB is set to print it out, see below.
  ! IUNIT  :IN       Integer.
  !         Fortran logical I/O device unit number to write the matrix
  !         to.  This unit must be connected in a system dependent fashion
  !         to a file or the console or you will get a nasty message
  !         from the Fortran I/O libraries.
  ! JOB    :IN       Integer.
  !         Flag indicating what I/O operations to perform.
  !         JOB = 0 => Print only the matrix.
  !             = 1 => Print matrix and RHS.
  !             = 2 => Print matrix and SOLN.
  !             = 3 => Print matrix, RHS and SOLN.
  !
  !- Description:
  !       The format for the output is as follows.  On  the first line
  !       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT
  !       and ISYM are described above.  IRHS is  a flag indicating if
  !       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a
  !       flag indicating if the SOLN was written out  (1 is yes, 0 is
  !       no).  The format for the fist line is: 5i10.  Then comes the
  !       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format
  !       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes
  !       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1,
  !       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7.
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
  !- Cautions:
  !     This routine will attempt to write to the Fortran logical output
  !     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
  !     this logical unit is attached to a file or terminal before calling
  !     this routine with a non-zero value for IUNIT.  This routine does
  !     not check for the validity of a non-zero IUNIT unit number.
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   871119  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   920511  Added complete declaration section.  (WRB)
  !   930701  Updated CATEGORY section.  (FNF, WRB)
  
  !     .. Scalar Arguments ..
  INTEGER Isym, Iunit, Job, N, Nelt
  !     .. Array Arguments ..
  REAL A(Nelt), Rhs(N), Soln(N)
  INTEGER Ia(Nelt), Ja(Nelt)
  !     .. Local Scalars ..
  INTEGER i, irhs, isoln
  !* FIRST EXECUTABLE STATEMENT  STOUT
  !
  !         If RHS and SOLN are to be printed also.
  !         Write out the information heading.
  !
  irhs = 0
  isoln = 0
  IF ( Job==1.OR.Job==3 ) irhs = 1
  IF ( Job>1 ) isoln = 1
  WRITE (Iunit,99001) N, Nelt, Isym, irhs, isoln
  99001 FORMAT (5I10)
  !
  !         Write out the matrix non-zeros in Triad format.
  DO i = 1, Nelt
    WRITE (Iunit,99002) Ia(i), Ja(i), A(i)
    99002 FORMAT (1X,I5,1X,I5,1X,E16.7)
  END DO
  !
  !         If requested, write out the rhs.
  IF ( irhs==1 ) WRITE (Iunit,99003) (Rhs(i),i=1,N)
  !
  !         If requested, write out the solution.
  IF ( isoln==1 ) WRITE (Iunit,99003) (Soln(i),i=1,N)
  RETURN
  99003 FORMAT (1X,E16.7)
  !------------- LAST LINE OF STOUT FOLLOWS ----------------------------
END SUBROUTINE STOUT
