        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:16 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CAS_UPDATE__genmod
          INTERFACE 
            SUBROUTINE CAS_UPDATE(BPS1,PSY,ROS,CASD,CAPACC_LIQ,         &
     &CAPACC_SNOW,LAI,CAST)
              USE MODULE_SIB, ONLY :                                    &
     &          CAS_TYPE
              REAL(KIND=8), INTENT(IN) :: BPS1
              REAL(KIND=8), INTENT(IN) :: PSY
              REAL(KIND=8), INTENT(IN) :: ROS
              REAL(KIND=8), INTENT(IN) :: CASD
              REAL(KIND=8), INTENT(IN) :: CAPACC_LIQ
              REAL(KIND=8), INTENT(IN) :: CAPACC_SNOW
              REAL(KIND=8), INTENT(IN) :: LAI
              TYPE (CAS_TYPE), INTENT(INOUT) :: CAST
            END SUBROUTINE CAS_UPDATE
          END INTERFACE 
        END MODULE CAS_UPDATE__genmod
