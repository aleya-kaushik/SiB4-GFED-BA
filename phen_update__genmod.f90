        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:24 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHEN_UPDATE__genmod
          INTERFACE 
            SUBROUTINE PHEN_UPDATE(GNUM,LNUM,PNUM,IPFT,PHENCONT,POOLCONT&
     &,DAYNEW,DAYLMAX,DAYLEN,DAYLENDT,CLIM_CUPR,CLIM_PRECIP,CUPR,LSPR,  &
     &TMDF,TM,ASSIM,RSTFAC2,RSTFAC4,LAI,HYDROVT,PHENT,POOLDT,POOLLT,VMAX&
     &,C13ASSIM,PHYSCONT,FRACT)
              USE MODULE_SIB, ONLY :                                    &
     &          HYDROV_TYPE,                                            &
     &          PHEN_TYPE,                                              &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE,                                             &
     &          FRACT_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHEN_PARAM,                                             &
     &          POOL_PARAM,                                             &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: GNUM
              INTEGER(KIND=4), INTENT(IN) :: LNUM
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              INTEGER(KIND=4), INTENT(INOUT) :: IPFT
              TYPE (PHEN_PARAM), INTENT(IN) :: PHENCONT
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              LOGICAL(KIND=4), INTENT(IN) :: DAYNEW
              REAL(KIND=4), INTENT(IN) :: DAYLMAX
              REAL(KIND=8), INTENT(IN) :: DAYLEN
              REAL(KIND=8), INTENT(IN) :: DAYLENDT
              REAL(KIND=8), INTENT(IN) :: CLIM_CUPR
              REAL(KIND=8), INTENT(IN) :: CLIM_PRECIP
              REAL(KIND=8), INTENT(IN) :: CUPR
              REAL(KIND=8), INTENT(IN) :: LSPR
              REAL(KIND=8), INTENT(IN) :: TMDF
              REAL(KIND=8), INTENT(IN) :: TM
              REAL(KIND=8), INTENT(IN) :: ASSIM
              REAL(KIND=8), INTENT(IN) :: RSTFAC2
              REAL(KIND=8), INTENT(IN) :: RSTFAC4
              REAL(KIND=8), INTENT(IN) :: LAI
              TYPE (HYDROV_TYPE), INTENT(IN) :: HYDROVT
              TYPE (PHEN_TYPE), INTENT(INOUT) :: PHENT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              REAL(KIND=8), INTENT(INOUT) :: VMAX
              REAL(KIND=8), INTENT(IN) :: C13ASSIM
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              TYPE (FRACT_TYPE), INTENT(IN) :: FRACT
            END SUBROUTINE PHEN_UPDATE
          END INTERFACE 
        END MODULE PHEN_UPDATE__genmod
