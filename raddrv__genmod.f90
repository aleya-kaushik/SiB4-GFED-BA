        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:48 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADDRV__genmod
          INTERFACE 
            SUBROUTINE RADDRV(SWDOWN,SUNANG,RADVBC,RADVDC,RADNBC,RADNDC)
              REAL(KIND=8), INTENT(IN) :: SWDOWN
              REAL(KIND=8), INTENT(IN) :: SUNANG
              REAL(KIND=8), INTENT(INOUT) :: RADVBC
              REAL(KIND=8), INTENT(INOUT) :: RADVDC
              REAL(KIND=8), INTENT(INOUT) :: RADNBC
              REAL(KIND=8), INTENT(INOUT) :: RADNDC
            END SUBROUTINE RADDRV
          END INTERFACE 
        END MODULE RADDRV__genmod
