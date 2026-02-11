        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:28 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POOL_HET_RETR__genmod
          INTERFACE 
            SUBROUTINE POOL_HET_RETR(POOLCONT,ZM,WOPTZM,WSAT,SEAS_PRECIP&
     &,CLIM_PRECIP,CLIM_ASSIM,ASSIMD,ROOTF_LAY,PAWFRAC_LAY,TD_LAY,      &
     &SATFRAC_LAY,POOLDT,FRACT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NSOIL,                                                  &
     &          NPOOLPFT,                                               &
     &          NPOOLLU,                                                &
     &          NPOOLSFC,                                               &
     &          NPOOLSOIL,                                              &
     &          NPOOLSFCC13,                                            &
     &          NPOOLSOILC13
              USE MODULE_SIB, ONLY :                                    &
     &          POOLD_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          POOL_PARAM
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: ZM
              REAL(KIND=8), INTENT(IN) :: WOPTZM
              REAL(KIND=8), INTENT(IN) :: WSAT
              REAL(KIND=8), INTENT(IN) :: SEAS_PRECIP
              REAL(KIND=8), INTENT(IN) :: CLIM_PRECIP
              REAL(KIND=8), INTENT(IN) :: CLIM_ASSIM
              REAL(KIND=8), INTENT(IN) :: ASSIMD
              REAL(KIND=8), INTENT(IN) :: ROOTF_LAY(10)
              REAL(KIND=8), INTENT(IN) :: PAWFRAC_LAY(10)
              REAL(KIND=8), INTENT(IN) :: TD_LAY(10)
              REAL(KIND=8), INTENT(IN) :: SATFRAC_LAY(10)
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE POOL_HET_RETR
          END INTERFACE 
        END MODULE POOL_HET_RETR__genmod
