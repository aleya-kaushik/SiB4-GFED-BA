        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:24 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHEN_GSS__genmod
          INTERFACE 
            SUBROUTINE PHEN_GSS(PHENCONT,DAYLENMAX,DAYLEN,DAYLENDT,PHENT&
     &)
              USE MODULE_SIB, ONLY :                                    &
     &          PHEN_TYPE
              USE MODULE_PARAM, ONLY :                                  &
     &          PHEN_PARAM
              TYPE (PHEN_PARAM), INTENT(IN) :: PHENCONT
              REAL(KIND=4), INTENT(IN) :: DAYLENMAX
              REAL(KIND=8), INTENT(IN) :: DAYLEN
              REAL(KIND=8), INTENT(IN) :: DAYLENDT
              TYPE (PHEN_TYPE), INTENT(INOUT) :: PHENT
            END SUBROUTINE PHEN_GSS
          END INTERFACE 
        END MODULE PHEN_GSS__genmod
