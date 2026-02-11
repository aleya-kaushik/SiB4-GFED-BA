        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:28 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_GRAZE__genmod
          INTERFACE 
            SUBROUTINE POOL_GRAZE(POOLCONT,CLIM_LAI,LAI,POOLLT,POOLDT,  &
     &FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: CLIM_LAI
              REAL(KIND=8), INTENT(IN) :: LAI
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_GRAZE
          END INTERFACE 
        END MODULE POOL_GRAZE__genmod
