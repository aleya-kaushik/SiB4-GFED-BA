        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_SOIL__genmod
          INTERFACE 
            SUBROUTINE HYDRO_SOIL(ECT,EGS,EGSMAX,ROOTF,WP_EFF,SOILT,    &
     &HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE,                                              &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              REAL(KIND=8), INTENT(IN) :: ECT
              REAL(KIND=8), INTENT(IN) :: EGS
              REAL(KIND=8), INTENT(IN) :: EGSMAX
              REAL(KIND=8), INTENT(IN) :: ROOTF(10)
              REAL(KIND=8), INTENT(IN) :: WP_EFF(10)
              TYPE (SOIL_TYPE), INTENT(IN) :: SOILT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE HYDRO_SOIL
          END INTERFACE 
        END MODULE HYDRO_SOIL__genmod
