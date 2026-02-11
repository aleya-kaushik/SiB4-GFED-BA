        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:02 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BALAN_CARBON__genmod
          INTERFACE 
            SUBROUTINE BALAN_CARBON(SIBPT,SIBLON,SIBLAT,PREF,ASSIMIN,   &
     &C13ASSIMIN,LAIIN,FPARIN,POOLDT,POOLLT,GRZ_TRANSFER,CO2T,FRACT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLCAN,                                               &
     &          NPOOLLU,                                                &
     &          NPOOLPFT,                                               &
     &          CARBONB_PRINT,                                          &
     &          CARBONB_STOP,                                           &
     &          CARBONB_THRESH,                                         &
     &          CARBONB_THRESHC13,                                      &
     &          NPOOLCANC13
              USE MODULE_SIB, ONLY :                                    &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE,                                             &
     &          CO2_TYPE,                                               &
     &          FRACT_TYPE
              INTEGER(KIND=4), INTENT(IN) :: SIBPT
              REAL(KIND=4), INTENT(IN) :: SIBLON
              REAL(KIND=4), INTENT(IN) :: SIBLAT
              INTEGER(KIND=4), INTENT(IN) :: PREF
              REAL(KIND=8), INTENT(IN) :: ASSIMIN
              REAL(KIND=8), INTENT(IN) :: C13ASSIMIN
              REAL(KIND=8), INTENT(IN) :: LAIIN
              REAL(KIND=8), INTENT(IN) :: FPARIN
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              REAL(KIND=4), INTENT(IN) :: GRZ_TRANSFER(NPOOLLU+2)
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE BALAN_CARBON
          END INTERFACE 
        END MODULE BALAN_CARBON__genmod
