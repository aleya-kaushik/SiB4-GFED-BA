        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_NVEG__genmod
          INTERFACE 
            SUBROUTINE HYDRO_NVEG(TM,TCAS,CUPRT,LSPRT,HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              REAL(KIND=8), INTENT(IN) :: TM
              REAL(KIND=8), INTENT(IN) :: TCAS
              REAL(KIND=8), INTENT(INOUT) :: CUPRT
              REAL(KIND=8), INTENT(INOUT) :: LSPRT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE HYDRO_NVEG
          END INTERFACE 
        END MODULE HYDRO_NVEG__genmod
