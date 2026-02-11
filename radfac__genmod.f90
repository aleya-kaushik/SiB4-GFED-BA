        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADFAC__genmod
          INTERFACE 
            SUBROUTINE RADFAC(CHIL,REF,TRAN,Z1,Z2,COSZ,SOILT,SSCOLT,VEGT&
     &,HYDROST,RADT)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE,                                              &
     &          SSCOL_TYPE,                                             &
     &          VEG_TYPE,                                               &
     &          HYDROS_TYPE,                                            &
     &          RAD_TYPE
              REAL(KIND=4), INTENT(IN) :: CHIL
              REAL(KIND=4), INTENT(IN) :: REF(2,2)
              REAL(KIND=4), INTENT(IN) :: TRAN(2,2)
              REAL(KIND=4), INTENT(IN) :: Z1
              REAL(KIND=4), INTENT(IN) :: Z2
              REAL(KIND=8), INTENT(IN) :: COSZ
              TYPE (SOIL_TYPE), INTENT(IN) :: SOILT
              TYPE (SSCOL_TYPE), INTENT(IN) :: SSCOLT
              TYPE (VEG_TYPE), INTENT(IN) :: VEGT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (RAD_TYPE), INTENT(INOUT) :: RADT
            END SUBROUTINE RADFAC
          END INTERFACE 
        END MODULE RADFAC__genmod
