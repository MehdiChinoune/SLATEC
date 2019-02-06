!*==SORTH.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SORTH
      SUBROUTINE SORTH(Vnew,V,Hes,N,Ll,Ldhes,Kmp,Snormw)
      IMPLICIT NONE
!*--SORTH5
!***BEGIN PROLOGUE  SORTH
!***SUBSIDIARY
!***PURPOSE  Internal routine for SGMRES.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      SINGLE PRECISION (SORTH-S, DORTH-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This routine  orthogonalizes  the  vector  VNEW  against the
!        previous KMP  vectors in the   V array.  It uses  a modified
!        Gram-Schmidt   orthogonalization procedure with  conditional
!        reorthogonalization.
!
! *Usage:
!      INTEGER N, LL, LDHES, KMP
!      REAL VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
!
!      CALL SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!
! *Arguments:
! VNEW   :INOUT    Real VNEW(N)
!         On input, the vector of length N containing a scaled
!         product of the Jacobian and the vector V(*,LL).
!         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
!         where i0 = max(1, LL-KMP+1).
! V      :IN       Real V(N,LL)
!         The N x LL array containing the previous LL
!         orthogonal vectors V(*,1) to V(*,LL).
! HES    :INOUT    Real HES(LDHES,LL)
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
! SNORMW :OUT      REAL
!         Scalar containing the l-2 norm of VNEW.
!
!***SEE ALSO  SGMRES
!***ROUTINES CALLED  SAXPY, SDOT, SNRM2
!***REVISION HISTORY  (YYMMDD)
!   871001  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to SGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  SORTH
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
      REAL Snormw
      INTEGER Kmp , Ldhes , Ll , N
!     .. Array Arguments ..
      REAL Hes(Ldhes,*) , V(N,*) , Vnew(*)
!     .. Local Scalars ..
      REAL arg , sumdsq , tem , vnrm
      INTEGER i , i0
!     .. External Functions ..
      REAL SDOT , SNRM2
      EXTERNAL SDOT , SNRM2
!     .. External Subroutines ..
      EXTERNAL SAXPY
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!***FIRST EXECUTABLE STATEMENT  SORTH
!
!         Get norm of unaltered VNEW for later use.
!
      vnrm = SNRM2(N,Vnew,1)
!   -------------------------------------------------------------------
!         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
!         Scaled inner products give new column of HES.
!         Projections of earlier vectors are subtracted from VNEW.
!   -------------------------------------------------------------------
      i0 = MAX(1,Ll-Kmp+1)
      DO i = i0 , Ll
        Hes(i,Ll) = SDOT(N,V(1,i),1,Vnew,1)
        tem = -Hes(i,Ll)
        CALL SAXPY(N,tem,V(1,i),1,Vnew,1)
      ENDDO
!   -------------------------------------------------------------------
!         Compute SNORMW = norm of VNEW.  If VNEW is small compared
!         to its input value (in norm), then reorthogonalize VNEW to
!         V(*,1) through V(*,LL).  Correct if relative correction
!         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
!         the dot products involved.
!   -------------------------------------------------------------------
      Snormw = SNRM2(N,Vnew,1)
      IF ( vnrm+0.001E0*Snormw/=vnrm ) RETURN
      sumdsq = 0
      DO i = i0 , Ll
        tem = -SDOT(N,V(1,i),1,Vnew,1)
        IF ( Hes(i,Ll)+0.001E0*tem/=Hes(i,Ll) ) THEN
          Hes(i,Ll) = Hes(i,Ll) - tem
          CALL SAXPY(N,tem,V(1,i),1,Vnew,1)
          sumdsq = sumdsq + tem**2
        ENDIF
      ENDDO
      IF ( sumdsq==0.0E0 ) RETURN
      arg = MAX(0.0E0,Snormw**2-sumdsq)
      Snormw = SQRT(arg)
!
!------------- LAST LINE OF SORTH FOLLOWS ----------------------------
      END SUBROUTINE SORTH
