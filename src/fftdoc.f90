!*==FFTDOC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK FFTDOC
      SUBROUTINE FFTDOC
      IMPLICIT NONE
!*--FFTDOC5
!***BEGIN PROLOGUE  FFTDOC
!***PURPOSE  Documentation for FFTPACK, a collection of Fast Fourier
!            Transform routines.
!***LIBRARY   SLATEC
!***CATEGORY  J1, Z
!***TYPE      ALL (FFTDOC-A)
!***KEYWORDS  DOCUMENTATION, FAST FOURIER TRANSFORM, FFT
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                       Version 3  June 1979
!
!          A Package of Fortran Subprograms for The Fast Fourier
!           Transform of Periodic and Other Symmetric Sequences
!                              By
!                       Paul N Swarztrauber
!
!    National Center For Atmospheric Research, Boulder, Colorado 80307
!        which is sponsored by the National Science Foundation
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     This package consists of programs which perform Fast Fourier
!     Transforms for both complex and real periodic sequences and
!     certain other symmetric sequences that are listed below.
!
!     1.   RFFTI     Initialize RFFTF and RFFTB
!     2.   RFFTF     Forward transform of a real periodic sequence
!     3.   RFFTB     Backward transform of a real coefficient array
!
!     4.   EZFFTI    Initialize EZFFTF and EZFFTB
!     5.   EZFFTF    A simplified real periodic forward transform
!     6.   EZFFTB    A simplified real periodic backward transform
!
!     7.   SINTI     Initialize SINT
!     8.   SINT      Sine transform of a real odd sequence
!
!     9.   COSTI     Initialize COST
!     10.  COST      Cosine transform of a real even sequence
!
!     11.  SINQI     Initialize SINQF and SINQB
!     12.  SINQF     Forward sine transform with odd wave numbers
!     13.  SINQB     Unnormalized inverse of SINQF
!
!     14.  COSQI     Initialize COSQF and COSQB
!     15.  COSQF     Forward cosine transform with odd wave numbers
!     16.  COSQB     Unnormalized inverse of COSQF
!
!     17.  CFFTI     Initialize CFFTF and CFFTB
!     18.  CFFTF     Forward transform of a complex periodic sequence
!     19.  CFFTB     Unnormalized inverse of CFFTF
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780201  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900723  PURPOSE section revised.  (WRB)
!***END PROLOGUE  FFTDOC
!***FIRST EXECUTABLE STATEMENT  FFTDOC
      END SUBROUTINE FFTDOC
