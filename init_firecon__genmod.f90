        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:07 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_FIRECON__genmod
          INTERFACE 
            SUBROUTINE INIT_FIRECON(NPFT,NPOOLPFT,NPOOLLU,FIRECON)
              USE MODULE_PARAM, ONLY :                                  &
     &          FIRE_PARAM
              INTEGER(KIND=4), INTENT(IN) :: NPFT
              INTEGER(KIND=4), INTENT(IN) :: NPOOLPFT
              INTEGER(KIND=4), INTENT(IN) :: NPOOLLU
              TYPE (FIRE_PARAM), INTENT(INOUT) :: FIRECON(NPFT)
            END SUBROUTINE INIT_FIRECON
          END INTERFACE 
        END MODULE INIT_FIRECON__genmod
