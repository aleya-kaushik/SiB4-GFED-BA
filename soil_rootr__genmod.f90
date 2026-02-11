        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOIL_ROOTR__genmod
          INTERFACE 
            SUBROUTINE SOIL_ROOTR(FIELDCAP,ECT,EGS,EGSMAX,ROOTF,WP_EFF, &
     &TD,VOL_LIQ,WWW_LIQ,ROOTR)
              REAL(KIND=8), INTENT(IN) :: FIELDCAP
              REAL(KIND=8), INTENT(IN) :: ECT
              REAL(KIND=8), INTENT(IN) :: EGS
              REAL(KIND=8), INTENT(IN) :: EGSMAX
              REAL(KIND=8), INTENT(IN) :: ROOTF(10)
              REAL(KIND=8), INTENT(IN) :: WP_EFF(10)
              REAL(KIND=8), INTENT(IN) :: TD(10)
              REAL(KIND=8), INTENT(IN) :: VOL_LIQ(10)
              REAL(KIND=8), INTENT(IN) :: WWW_LIQ(10)
              REAL(KIND=8), INTENT(INOUT) :: ROOTR(10)
            END SUBROUTINE SOIL_ROOTR
          END INTERFACE 
        END MODULE SOIL_ROOTR__genmod
