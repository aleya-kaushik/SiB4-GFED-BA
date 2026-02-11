        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 26 16:45:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CLM_TRIDIA__genmod
          INTERFACE 
            SUBROUTINE CLM_TRIDIA(N,A,B,C,R,U)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A(1:N)
              REAL(KIND=8), INTENT(IN) :: B(1:N)
              REAL(KIND=8), INTENT(IN) :: C(1:N)
              REAL(KIND=8), INTENT(IN) :: R(1:N)
              REAL(KIND=8), INTENT(OUT) :: U(1:N)
            END SUBROUTINE CLM_TRIDIA
          END INTERFACE 
        END MODULE CLM_TRIDIA__genmod
