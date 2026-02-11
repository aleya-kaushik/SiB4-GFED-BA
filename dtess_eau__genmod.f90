        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DTESS_EAU__genmod
          INTERFACE 
            SUBROUTINE DTESS_EAU(LEN,PL,TL,ESS,DTESS)
              INTEGER(KIND=4), INTENT(IN) :: LEN
              REAL(KIND=8), INTENT(IN) :: PL(LEN)
              REAL(KIND=8), INTENT(IN) :: TL(LEN)
              REAL(KIND=8), INTENT(OUT) :: ESS(LEN)
              REAL(KIND=8), INTENT(OUT) :: DTESS(LEN)
            END SUBROUTINE DTESS_EAU
          END INTERFACE 
        END MODULE DTESS_EAU__genmod
