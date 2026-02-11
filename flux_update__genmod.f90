        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLUX_UPDATE__genmod
          INTERFACE 
            SUBROUTINE FLUX_UPDATE(PSY,ROS,RADTC,RADTG,RADTS,RADC3C,    &
     &RADC3G,POROS,PAW_LAY,CAST,FLUXT,HYDROST,SSCOLT,PRESS,SPDM,DLWBOT)
              USE MODULE_SIB, ONLY :                                    &
     &          CAS_TYPE,                                               &
     &          FLUX_TYPE,                                              &
     &          HYDROS_TYPE,                                            &
     &          SSCOL_TYPE
              REAL(KIND=8), INTENT(IN) :: PSY
              REAL(KIND=8), INTENT(IN) :: ROS
              REAL(KIND=8), INTENT(IN) :: RADTC
              REAL(KIND=8), INTENT(IN) :: RADTG
              REAL(KIND=8), INTENT(IN) :: RADTS
              REAL(KIND=8), INTENT(IN) :: RADC3C
              REAL(KIND=8), INTENT(IN) :: RADC3G
              REAL(KIND=8), INTENT(IN) :: POROS
              REAL(KIND=8), INTENT(IN) :: PAW_LAY(10)
              TYPE (CAS_TYPE), INTENT(INOUT) :: CAST
              TYPE (FLUX_TYPE), INTENT(INOUT) :: FLUXT
              TYPE (HYDROS_TYPE), INTENT(INOUT) :: HYDROST
              TYPE (SSCOL_TYPE), INTENT(INOUT) :: SSCOLT
              REAL(KIND=8), INTENT(IN) :: PRESS
              REAL(KIND=8), INTENT(IN) :: SPDM
              REAL(KIND=8), INTENT(IN) :: DLWBOT
            END SUBROUTINE FLUX_UPDATE
          END INTERFACE 
        END MODULE FLUX_UPDATE__genmod
