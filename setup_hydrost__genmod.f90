        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_HYDROST__genmod
          INTERFACE 
            SUBROUTINE SETUP_HYDROST(Z1,Z2,SSCOLT,HYDROST)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE,                                              &
     &          SSCOL_TYPE,                                             &
     &          HYDROS_TYPE
              REAL(KIND=4), INTENT(IN) :: Z1
              REAL(KIND=4), INTENT(IN) :: Z2
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
            END SUBROUTINE SETUP_HYDROST
          END INTERFACE 
        END MODULE SETUP_HYDROST__genmod
