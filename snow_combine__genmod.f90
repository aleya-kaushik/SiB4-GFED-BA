        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:32 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SNOW_COMBINE__genmod
          INTERFACE 
            SUBROUTINE SNOW_COMBINE(GREF,GLON,GLAT,PREF,LAI,POROS,      &
     &HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: GLON
              REAL(KIND=4), INTENT(IN) :: GLAT
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: POROS
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE SNOW_COMBINE
          END INTERFACE 
        END MODULE SNOW_COMBINE__genmod
