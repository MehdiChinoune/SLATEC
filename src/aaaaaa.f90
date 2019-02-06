!*==AAAAAA.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK AAAAAA
      SUBROUTINE AAAAAA(Ver)
      IMPLICIT NONE
!*--AAAAAA5
!***BEGIN PROLOGUE  AAAAAA
!***PURPOSE  SLATEC Common Mathematical Library disclaimer and version.
!***LIBRARY   SLATEC
!***CATEGORY  Z
!***TYPE      ALL (AAAAAA-A)
!***KEYWORDS  DISCLAIMER, DOCUMENTATION, VERSION
!***AUTHOR  SLATEC Common Mathematical Library Committee
!***DESCRIPTION
!
!   The SLATEC Common Mathematical Library is issued by the following
!
!           Air Force Weapons Laboratory, Albuquerque
!           Lawrence Livermore National Laboratory, Livermore
!           Los Alamos National Laboratory, Los Alamos
!           National Institute of Standards and Technology, Washington
!           National Energy Research Supercomputer Center, Livermore
!           Oak Ridge National Laboratory, Oak Ridge
!           Sandia National Laboratories, Albuquerque
!           Sandia National Laboratories, Livermore
!
!   All questions concerning the distribution of the library should be
!   directed to the NATIONAL ENERGY SOFTWARE CENTER, 9700 Cass Ave.,
!   Argonne, Illinois  60439, and not to the authors of the subprograms.
!
!                    * * * * * Notice * * * * *
!
!   This material was prepared as an account of work sponsored by the
!   United States Government.  Neither the United States, nor the
!   Department of Energy, nor the Department of Defense, nor any of
!   their employees, nor any of their contractors, subcontractors, or
!   their employees, makes any warranty, expressed or implied, or
!   assumes any legal liability or responsibility for the accuracy,
!   completeness, or usefulness of any information, apparatus, product,
!   or process disclosed, or represents that its use would not infringe
!   upon privately owned rights.
!
! *Usage:
!
!        CHARACTER * 16 VER
!
!        CALL AAAAAA (VER)
!
! *Arguments:
!
!     VER:OUT   will contain the version number of the SLATEC CML.
!
! *Description:
!
!   This routine contains the SLATEC Common Mathematical Library
!   disclaimer and can be used to return the library version number.
!
!***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
!                 and Lee Walton, Guide to the SLATEC Common Mathema-
!                 tical Library, April 10, 1990.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800424  DATE WRITTEN
!   890414  REVISION DATE from Version 3.2
!   890713  Routine modified to return version number.  (WRB)
!   900330  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   921215  Updated for Version 4.0.  (WRB)
!   930701  Updated for Version 4.1.  (WRB)
!***END PROLOGUE  AAAAAA
      CHARACTER*(*) Ver
!***FIRST EXECUTABLE STATEMENT  AAAAAA
      Ver = ' 4.1'
      END SUBROUTINE AAAAAA
