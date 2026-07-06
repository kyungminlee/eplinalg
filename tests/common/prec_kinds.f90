module prec_kinds
    implicit none
    private
    public :: ep, dp

    ! Universal reference precision for the differential precision
    ! tests — every comparison happens in REAL(KIND=ep) regardless of
    ! which target is being tested.
    integer, parameter :: ep = 16

    ! Double precision — the kind for interfaces that stay REAL(KIND=8)
    ! regardless of the working precision: multifloats wrappers bridging
    ! REAL(KIND=8) and TYPE(real64x2), and fixed-ABI values such as
    ! MPI_WTIME results in the MUMPS drivers.
    integer, parameter :: dp = 8
end module prec_kinds
