        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:46:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SETUP_SOILT__genmod
          INTERFACE 
            SUBROUTINE SETUP_SOILT(CLAYFRAC,SANDFRAC,SOREF_VIS,SOREF_NIR&
     &,FC_MIN,WP_MIN,SOILT)
              USE MODULE_SIB, ONLY :                                    &
     &          SOIL_TYPE
              REAL(KIND=8), INTENT(IN) :: CLAYFRAC
              REAL(KIND=8), INTENT(IN) :: SANDFRAC
              REAL(KIND=8), INTENT(IN) :: SOREF_VIS
              REAL(KIND=8), INTENT(IN) :: SOREF_NIR
              REAL(KIND=4), INTENT(IN) :: FC_MIN(5)
              REAL(KIND=4), INTENT(IN) :: WP_MIN(5)
              TYPE (SOIL_TYPE), INTENT(INOUT) :: SOILT
            END SUBROUTINE SETUP_SOILT
          END INTERFACE 
        END MODULE SETUP_SOILT__genmod
