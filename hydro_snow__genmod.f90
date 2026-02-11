        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_SNOW__genmod
          INTERFACE 
            SUBROUTINE HYDRO_SNOW(POROS,HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              REAL(KIND=8), INTENT(IN) :: POROS
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE HYDRO_SNOW
          END INTERFACE 
        END MODULE HYDRO_SNOW__genmod
