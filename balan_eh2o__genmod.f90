        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BALAN_EH2O__genmod
          INTERFACE 
            SUBROUTINE BALAN_EH2O(GREF,GLON,GLAT,PREF,LSPR,CUPR,LAI,    &
     &FLUXT,SSCOLT,HYDROST,EBALNUM,WBALNUM)
              USE MODULE_SIB, ONLY :                                    &
     &          FLUX_TYPE,                                              &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: GLON
              REAL(KIND=4), INTENT(IN) :: GLAT
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=8), INTENT(IN) :: LSPR
              REAL(KIND=8), INTENT(IN) :: CUPR
              REAL(KIND=8), INTENT(IN) :: LAI
              TYPE (FLUX_TYPE), INTENT(IN) :: FLUXT
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              INTEGER(KIND=4), INTENT(INOUT) :: EBALNUM
              INTEGER(KIND=4), INTENT(INOUT) :: WBALNUM
            END SUBROUTINE BALAN_EH2O
          END INTERFACE 
        END MODULE BALAN_EH2O__genmod
