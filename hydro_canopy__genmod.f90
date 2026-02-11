        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:21 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_CANOPY__genmod
          INTERFACE 
            SUBROUTINE HYDRO_CANOPY(GREF,GLON,GLAT,PREF,CHIL,TM,TCAS,LAI&
     &,VCOVER,CUPRT,LSPRT,TC,HYDROST,SSCOLT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: GLON
              REAL(KIND=4), INTENT(IN) :: GLAT
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=4), INTENT(IN) :: CHIL
              REAL(KIND=8), INTENT(IN) :: TM
              REAL(KIND=8), INTENT(IN) :: TCAS
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: VCOVER
              REAL(KIND=8), INTENT(INOUT) :: CUPRT
              REAL(KIND=8), INTENT(INOUT) :: LSPRT
              REAL(KIND=8), INTENT(INOUT) :: TC
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
            END SUBROUTINE HYDRO_CANOPY
          END INTERFACE 
        END MODULE HYDRO_CANOPY__genmod
