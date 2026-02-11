        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:10 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_POOLCON__genmod
          INTERFACE 
            SUBROUTINE INIT_POOLCON(NPFT,NPOOLPFT,NPOOLLU,POOLCON)
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              INTEGER(KIND=4), INTENT(IN) :: NPFT
              INTEGER(KIND=4), INTENT(IN) :: NPOOLPFT
              INTEGER(KIND=4), INTENT(IN) :: NPOOLLU
              TYPE (POOL_PARAM), INTENT(INOUT) :: POOLCON(NPFT)
            END SUBROUTINE INIT_POOLCON
          END INTERFACE 
        END MODULE INIT_POOLCON__genmod
