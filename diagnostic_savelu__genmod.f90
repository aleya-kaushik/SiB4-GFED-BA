        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 13 14:42:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_SAVELU__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_SAVELU(PNUM,NVARS,NSAVE,REFSAVE,      &
     &REFVARS,SOILT,CAST,CO2T,COST,FRACT,FLUXT,HYDROST,HYDROVT,PHENT,   &
     &POOLDT,POOLLT,RADT,SIFT,SSCOLT,VEGT,OUTSAVE)
              USE MODULE_SIB
              INTEGER(KIND=4), INTENT(IN) :: NSAVE
              INTEGER(KIND=4), INTENT(IN) :: NVARS
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              INTEGER(KIND=4), INTENT(IN) :: REFSAVE
              INTEGER(KIND=4), INTENT(IN) :: REFVARS(NVARS)
              TYPE (SOIL_TYPE), INTENT(IN) :: SOILT
              TYPE (CAS_TYPE), INTENT(IN) :: CAST
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (COS_TYPE), INTENT(IN) :: COST
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
              TYPE (FLUX_TYPE), INTENT(IN) :: FLUXT
              TYPE (HYDROS_TYPE), INTENT(IN) :: HYDROST
              TYPE (HYDROV_TYPE), INTENT(IN) :: HYDROVT
              TYPE (PHEN_TYPE), INTENT(IN) :: PHENT
              TYPE (POOLD_TYPE), INTENT(IN) :: POOLDT
              TYPE (POOLL_TYPE), INTENT(IN) :: POOLLT
              TYPE (RAD_TYPE), INTENT(IN) :: RADT
              TYPE (SIF_TYPE), INTENT(IN) :: SIFT
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
              TYPE (VEG_TYPE), INTENT(IN) :: VEGT
              REAL(KIND=8), INTENT(INOUT) :: OUTSAVE(NVARS,NSAVE)
            END SUBROUTINE DIAGNOSTIC_SAVELU
          END INTERFACE 
        END MODULE DIAGNOSTIC_SAVELU__genmod
