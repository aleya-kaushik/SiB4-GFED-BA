        !COMPILER-GENERATED INTERFACE MODULE: Tue Jan 13 14:42:01 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIAGNOSTIC_SAVE_LALL__genmod
          INTERFACE 
            SUBROUTINE DIAGNOSTIC_SAVE_LALL(GNUM,LNUM,PNUM,SOILT,CAST,  &
     &CO2T,COST,FRACT,FLUXT,HYDROST,HYDROVT,PHENT,POOLDT,POOLLT,RADT,   &
     &SIFT,SSCOLT,VEGT)
              USE MODULE_SIB
              INTEGER(KIND=4), INTENT(IN) :: GNUM
              INTEGER(KIND=4), INTENT(IN) :: LNUM
              INTEGER(KIND=4), INTENT(IN) :: PNUM
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
            END SUBROUTINE DIAGNOSTIC_SAVE_LALL
          END INTERFACE 
        END MODULE DIAGNOSTIC_SAVE_LALL__genmod
