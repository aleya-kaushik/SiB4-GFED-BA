        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_PHENT__genmod
          INTERFACE 
            SUBROUTINE SETUP_PHENT(GNUM,LNUM,PNUM,PREF,PHENCONT,POOLCONT&
     &,DAYLEN,DAYLENDT,DAYLENMAX,LAI,TM,VMAX,POOLDT,POOLLT,PHENT,       &
     &PHYSCONT)
              USE MODULE_SIB, ONLY :                                    &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE,                                             &
     &          PHEN_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHEN_PARAM,                                             &
     &          POOL_PARAM,                                             &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: GNUM
              INTEGER(KIND=4), INTENT(IN) :: LNUM
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              INTEGER(KIND=4), INTENT(IN) :: PREF
              TYPE (PHEN_PARAM), INTENT(IN) :: PHENCONT
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: DAYLEN
              REAL(KIND=8), INTENT(IN) :: DAYLENDT
              REAL(KIND=4), INTENT(IN) :: DAYLENMAX
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: TM
              REAL(KIND=8), INTENT(INOUT) :: VMAX
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (PHEN_TYPE), INTENT(INOUT) :: PHENT
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
            END SUBROUTINE SETUP_PHENT
          END INTERFACE 
        END MODULE SETUP_PHENT__genmod
