        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:29 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_UPDATE__genmod
          INTERFACE 
            SUBROUTINE POOL_UPDATE(SIBPT,SIBLON,SIBLAT,PREF,ASSIMIN,    &
     &C13ASSIMIN,LAIIN,FPARIN,EQUIBDT,POOLDT,EQUIBLT,POOLLT,GRZ_TRANSFER&
     &,CO2T,FRACT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLLU,                                                &
     &          NPOOLPFT
              USE MODULE_SIB, ONLY :                                    &
     &          EQUIBD_TYPE,                                            &
     &          POOLD_TYPE,                                             &
     &          EQUIBL_TYPE,                                            &
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
              TYPE (EQUIBD_TYPE), INTENT(INOUT) :: EQUIBDT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (EQUIBL_TYPE), INTENT(INOUT) :: EQUIBLT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              REAL(KIND=4), INTENT(IN) :: GRZ_TRANSFER(NPOOLLU+2)
              TYPE (CO2_TYPE), INTENT(IN) :: CO2T
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_UPDATE
          END INTERFACE 
        END MODULE POOL_UPDATE__genmod
