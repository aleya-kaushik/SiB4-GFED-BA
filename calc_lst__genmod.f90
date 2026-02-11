        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALC_LST__genmod
          INTERFACE 
            SUBROUTINE CALC_LST(DAYFRAC,SUBCOUNT,LONSIB,DAYFRAC_LST)
              INTEGER(KIND=4), INTENT(IN) :: SUBCOUNT
              REAL(KIND=8), INTENT(IN) :: DAYFRAC
              REAL(KIND=4), INTENT(IN) :: LONSIB(SUBCOUNT)
              REAL(KIND=8), INTENT(INOUT) :: DAYFRAC_LST(SUBCOUNT)
            END SUBROUTINE CALC_LST
          END INTERFACE 
        END MODULE CALC_LST__genmod
