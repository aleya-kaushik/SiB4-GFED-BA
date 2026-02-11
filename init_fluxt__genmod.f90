        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT_FLUXT__genmod
          INTERFACE 
            SUBROUTINE INIT_FLUXT(FLUXT)
              USE MODULE_SIB, ONLY :                                    &
     &          FLUX_TYPE
              TYPE (FLUX_TYPE), INTENT(INOUT) :: FLUXT
            END SUBROUTINE INIT_FLUXT
          END INTERFACE 
        END MODULE INIT_FLUXT__genmod
