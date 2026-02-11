        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:08 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_PHENCON__genmod
          INTERFACE 
            SUBROUTINE INIT_PHENCON(NPFT,NPOOLPFT,NPSTGMAX,PHENCON)
              USE MODULE_PARAM, ONLY :                                  &
     &          PHEN_PARAM
              INTEGER(KIND=4) :: NPFT
              INTEGER(KIND=4) :: NPOOLPFT
              INTEGER(KIND=4) :: NPSTGMAX
              TYPE (PHEN_PARAM), INTENT(INOUT) :: PHENCON(NPFT)
            END SUBROUTINE INIT_PHENCON
          END INTERFACE 
        END MODULE INIT_PHENCON__genmod
