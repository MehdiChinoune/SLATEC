!** AAAAAA
SUBROUTINE AAAAAA(Ver)
  IMPLICIT NONE
  !>
  !***
  !  SLATEC Common Mathematical Library disclaimer and version.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  Z
  !***
  ! **Type:**      ALL (AAAAAA-A)
  !***
  ! **Keywords:**  DISCLAIMER, DOCUMENTATION, VERSION
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
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
  !- Usage:
  !
  !        CHARACTER * 16 VER
  !
  !        CALL AAAAAA (VER)
  !
  !- Arguments:
  !
  !     VER:OUT   will contain the version number of the SLATEC CML.
  !
  !- Description:
  !
  !   This routine contains the SLATEC Common Mathematical Library
  !   disclaimer and can be used to return the library version number.
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800424  DATE WRITTEN
  !   890414  REVISION DATE from Version 3.2
  !   890713  Routine modified to return version number.  (WRB)
  !   900330  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   921215  Updated for Version 4.0.  (WRB)
  !   930701  Updated for Version 4.1.  (WRB)
  
  CHARACTER*(*) Ver
  !* FIRST EXECUTABLE STATEMENT  AAAAAA
  Ver = ' 4.2'
END SUBROUTINE AAAAAA
