        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SNOW_COMPACT__genmod
          INTERFACE 
            SUBROUTINE SNOW_COMPACT(NSL,TD,WWW_LIQ,WWW_ICE,DZ)
              INTEGER(KIND=1), INTENT(IN) :: NSL
              REAL(KIND=8), INTENT(IN) :: TD(NSL+1:0)
              REAL(KIND=8), INTENT(IN) :: WWW_LIQ(NSL+1:0)
              REAL(KIND=8), INTENT(IN) :: WWW_ICE(NSL+1:0)
              REAL(KIND=8), INTENT(INOUT) :: DZ(NSL+1:0)
            END SUBROUTINE SNOW_COMPACT
          END INTERFACE 
        END MODULE SNOW_COMPACT__genmod
