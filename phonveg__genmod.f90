        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHONVEG__genmod
          INTERFACE 
            SUBROUTINE PHONVEG(RADVBC,RADVDC,GPROGT,CO2T,PS)
              USE MODULE_SIB, ONLY :                                    &
     &          CO2_TYPE,                                               &
     &          GPROG_TYPE
              REAL(KIND=8), INTENT(IN) :: RADVBC
              REAL(KIND=8), INTENT(IN) :: RADVDC
              TYPE (GPROG_TYPE), INTENT(IN) :: GPROGT
              TYPE (CO2_TYPE), INTENT(INOUT) :: CO2T
              REAL(KIND=8), INTENT(IN) :: PS
            END SUBROUTINE PHONVEG
          END INTERFACE 
        END MODULE PHONVEG__genmod
