        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HYDRO_SETV__genmod
          INTERFACE 
            SUBROUTINE HYDRO_SETV(FC_EFF,WP_EFF,ROOTF,SSCOLT,HYDROVT)
              USE MODULE_SIB, ONLY :                                    &
     &          SSCOL_TYPE,                                             &
     &          HYDROV_TYPE
              REAL(KIND=8), INTENT(IN) :: FC_EFF(10)
              REAL(KIND=8), INTENT(IN) :: WP_EFF(10)
              REAL(KIND=8), INTENT(IN) :: ROOTF(10)
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
              TYPE (HYDROV_TYPE), INTENT(INOUT) :: HYDROVT
            END SUBROUTINE HYDRO_SETV
          END INTERFACE 
        END MODULE HYDRO_SETV__genmod
