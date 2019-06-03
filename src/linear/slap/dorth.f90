!** DORTH
SUBROUTINE DORTH(Vnew,V,Hes,N,Ll,Ldhes,Kmp,Snormw)
  !>
  !  Internal routine for DGMRES.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (SORTH-S, DORTH-D)
  !***
  ! **Keywords:**  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
  !             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
  !***
  ! **Author:**  Brown, Peter, (LLNL), pnbrown@llnl.gov
  !           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
  !           Seager, Mark K., (LLNL), seager@llnl.gov
  !             Lawrence Livermore National Laboratory
  !             PO Box 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !***
  ! **Description:**
  !        This routine  orthogonalizes  the  vector  VNEW  against the
  !        previous KMP  vectors in the   V array.  It uses  a modified
  !        Gram-Schmidt   orthogonalization procedure with  conditional
  !        reorthogonalization.
  !
  !- Usage:
  !      INTEGER N, LL, LDHES, KMP
  !      DOUBLE PRECISION VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
  !
  !      CALL DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
  !
  !- Arguments:
  ! VNEW   :INOUT    Double Precision VNEW(N)
  !         On input, the vector of length N containing a scaled
  !         product of the Jacobian and the vector V(*,LL).
  !         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
  !         where i0 = max(1, LL-KMP+1).
  ! V      :IN       Double Precision V(N,LL)
  !         The N x LL array containing the previous LL
  !         orthogonal vectors V(*,1) to V(*,LL).
  ! HES    :INOUT    Double Precision HES(LDHES,LL)
  !         On input, an LL x LL upper Hessenberg matrix containing,
  !         in HES(I,K), K.lt.LL, the scaled inner products of
  !         A*V(*,K) and V(*,i).
  !         On return, column LL of HES is filled in with
  !         the scaled inner products of A*V(*,LL) and V(*,i).
  ! N      :IN       Integer
  !         The order of the matrix A, and the length of VNEW.
  ! LL     :IN       Integer
  !         The current order of the matrix HES.
  ! LDHES  :IN       Integer
  !         The leading dimension of the HES array.
  ! KMP    :IN       Integer
  !         The number of previous vectors the new vector VNEW
  !         must be made orthogonal to (KMP .le. MAXL).
  ! SNORMW :OUT      DOUBLE PRECISION
  !         Scalar containing the l-2 norm of VNEW.
  !
  !***
  ! **See also:**  DGMRES
  !***
  ! **Routines called:**  DAXPY, DDOT, DNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910506  Made subsidiary to DGMRES.  (FNF)
  !   920511  Added complete declaration section.  (WRB)

  !     .. Scalar Arguments ..
  REAL(DP) :: Snormw
  INTEGER Kmp, Ldhes, Ll, N
  !     .. Array Arguments ..
  REAL(DP) :: Hes(Ldhes,*), V(N,*), Vnew(*)
  !     .. Local Scalars ..
  REAL(DP) :: arg, sumdsq, tem, vnrm
  INTEGER i, i0
  !     .. Intrinsic Functions ..
  INTRINSIC MAX, SQRT
  !* FIRST EXECUTABLE STATEMENT  DORTH
  !
  !         Get norm of unaltered VNEW for later use.
  !
  vnrm = DNRM2(N,Vnew,1)
  !   -------------------------------------------------------------------
  !         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
  !         Scaled inner products give new column of HES.
  !         Projections of earlier vectors are subtracted from VNEW.
  !   -------------------------------------------------------------------
  i0 = MAX(1,Ll-Kmp+1)
  DO i = i0, Ll
    Hes(i,Ll) = DDOT(N,V(1,i),1,Vnew,1)
    tem = -Hes(i,Ll)
    CALL DAXPY(N,tem,V(1,i),1,Vnew,1)
  END DO
  !   -------------------------------------------------------------------
  !         Compute SNORMW = norm of VNEW.  If VNEW is small compared
  !         to its input value (in norm), then reorthogonalize VNEW to
  !         V(*,1) through V(*,LL).  Correct if relative correction
  !         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
  !         the dot products involved.
  !   -------------------------------------------------------------------
  Snormw = DNRM2(N,Vnew,1)
  IF ( vnrm+0.001D0*Snormw/=vnrm ) RETURN
  sumdsq = 0
  DO i = i0, Ll
    tem = -DDOT(N,V(1,i),1,Vnew,1)
    IF ( Hes(i,Ll)+0.001D0*tem/=Hes(i,Ll) ) THEN
      Hes(i,Ll) = Hes(i,Ll) - tem
      CALL DAXPY(N,tem,V(1,i),1,Vnew,1)
      sumdsq = sumdsq + tem**2
    END IF
  END DO
  IF ( sumdsq==0.0D0 ) RETURN
  arg = MAX(0.0D0,Snormw**2-sumdsq)
  Snormw = SQRT(arg)
  !
  !------------- LAST LINE OF DORTH FOLLOWS ----------------------------
END SUBROUTINE DORTH
