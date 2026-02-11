        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:32 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VEG_UPDATE__genmod
          INTERFACE 
            SUBROUTINE VEG_UPDATE(DOY,GREF,GLON,GLAT,PNUM,PREF,ISCROP,  &
     &ISGRASS,PHYSCONT,POOLLU,SNOW_CVFC,NODE_Z,POOLLT,VEGT)
              USE MODULE_SIBCONST, ONLY :                               &
     &          NPOOLPFT,                                               &
     &          NPOOLLU,                                                &
     &          NSOIL,                                                  &
     &          GREEN_SWITCH,                                           &
     &          PRINT_STOP,                                             &
     &          PRINT_VEG
              USE MODULE_SIB, ONLY :                                    &
     &          POOLL_TYPE,                                             &
     &          VEG_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: DOY
              INTEGER(KIND=4), INTENT(IN) :: GREF
              REAL(KIND=4), INTENT(IN) :: GLON
              REAL(KIND=4), INTENT(IN) :: GLAT
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              INTEGER(KIND=4), INTENT(IN) :: PREF
              LOGICAL(KIND=4), INTENT(IN) :: ISCROP
              LOGICAL(KIND=4), INTENT(IN) :: ISGRASS
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
              REAL(KIND=8), INTENT(IN) :: POOLLU(NPOOLLU)
              REAL(KIND=8), INTENT(IN) :: SNOW_CVFC
              REAL(KIND=8), INTENT(IN) :: NODE_Z(10)
              TYPE (POOLL_TYPE), INTENT(IN) :: POOLLT
              TYPE (VEG_TYPE), INTENT(INOUT) :: VEGT
            END SUBROUTINE VEG_UPDATE
          END INTERFACE 
        END MODULE VEG_UPDATE__genmod
