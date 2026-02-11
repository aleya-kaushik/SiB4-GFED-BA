        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:23 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHEN_DEFINED__genmod
          INTERFACE 
            SUBROUTINE PHEN_DEFINED(SUBPT,SUBL,IPFT,PNUM,PHENCONT,      &
     &POOLCONT,LAI,TMDF,PHENT,POOLDT,POOLLT,PHYSCONT)
              USE MODULE_SIB, ONLY :                                    &
     &          PHEN_TYPE,                                              &
     &          POOLD_TYPE,                                             &
     &          POOLL_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHEN_PARAM,                                             &
     &          POOL_PARAM,                                             &
     &          PHYS_PARAM
              INTEGER(KIND=4), INTENT(IN) :: SUBPT
              INTEGER(KIND=4), INTENT(IN) :: SUBL
              INTEGER(KIND=4), INTENT(INOUT) :: IPFT
              INTEGER(KIND=4), INTENT(IN) :: PNUM
              TYPE (PHEN_PARAM), INTENT(IN) :: PHENCONT
              TYPE (POOL_PARAM), INTENT(IN) :: POOLCONT
              REAL(KIND=8), INTENT(IN) :: LAI
              REAL(KIND=8), INTENT(IN) :: TMDF
              TYPE (PHEN_TYPE), INTENT(INOUT) :: PHENT
              TYPE (POOLD_TYPE), INTENT(INOUT) :: POOLDT
              TYPE (POOLL_TYPE), INTENT(INOUT) :: POOLLT
              TYPE (PHYS_PARAM), INTENT(IN) :: PHYSCONT
            END SUBROUTINE PHEN_DEFINED
          END INTERFACE 
        END MODULE PHEN_DEFINED__genmod
