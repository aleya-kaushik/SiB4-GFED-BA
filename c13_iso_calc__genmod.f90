        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:16 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE C13_ISO_CALC__genmod
          INTERFACE 
            SUBROUTINE C13_ISO_CALC(POOLCONT,POOLLT,POOLDT,FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          POOLL_TYPE,                                             &
     &          POOLD_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(INOUT) :: POOLCONT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (FRACT_TYPE), INTENT(INOUT) :: FRACT
            END SUBROUTINE C13_ISO_CALC
          END INTERFACE 
        END MODULE C13_ISO_CALC__genmod
