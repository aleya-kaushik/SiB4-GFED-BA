        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_FRACT__genmod
          INTERFACE 
            SUBROUTINE INIT_FRACT(FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          FRACT_TYPE
              TYPE (FRACT_TYPE), INTENT(INOUT) :: FRACT
            END SUBROUTINE INIT_FRACT
          END INTERFACE 
        END MODULE INIT_FRACT__genmod
