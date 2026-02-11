        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADNET__genmod
          INTERFACE 
            SUBROUTINE RADNET(NSL,TSFC,TC,RADT)
              USE MODULE_SIB, ONLY :                                    &
     &          RAD_TYPE
              INTEGER(KIND=1), INTENT(IN) :: NSL
              REAL(KIND=8), INTENT(IN) :: TSFC
              REAL(KIND=8), INTENT(IN) :: TC
              TYPE (RAD_TYPE), INTENT(INOUT) :: RADT
            END SUBROUTINE RADNET
          END INTERFACE 
        END MODULE RADNET__genmod
