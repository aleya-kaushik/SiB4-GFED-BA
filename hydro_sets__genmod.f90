        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_SETS__genmod
          INTERFACE 
            SUBROUTINE HYDRO_SETS(SOILT,HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE,                                              &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              TYPE (SOIL_TYPE), INTENT(IN) :: SOILT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE HYDRO_SETS
          END INTERFACE 
        END MODULE HYDRO_SETS__genmod
