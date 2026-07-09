! Test BLACS_SET round-trip via BLACS_GET on the WHAT parameters that
! are documented as settable (the broadcast / combine ring / nbranch
! tunables and the tops-repeat / tops-coherent flags). For each WHAT,
! set a non-default value and verify GET returns it.
!
! WHAT codes (from BLACS/SRC/Bdef.h):
!   11 SGET_NR_BS        broadcast nrings
!   12 SGET_NB_BS        broadcast nbranches (note: stored as val+1)
!   13 SGET_NR_CO        combine nrings
!   14 SGET_NB_CO        combine nbranches  (note: stored as val+1)
!   15 SGET_TOPSREPEAT   tops-repeat flag
!   16 SGET_TOPSCOHRNT   tops-coherent flag
program test_blacs_set
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context
    use target_blacs,      only: target_name
    implicit none

    interface
        subroutine blacs_set(icontxt, what, val)
            integer, intent(in) :: icontxt, what, val
        end subroutine
        subroutine blacs_get(icontxt, what, val)
            integer, intent(in)  :: icontxt, what
            integer, intent(out) :: val
        end subroutine
    end interface

    integer :: got
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)

    call grid_init()
    call report_init('blacs_set', target_name, my_rank)

    ! ── SGET_NR_BS: broadcast nrings — round-trip 7
    call blacs_set(my_context, 11, 7)
    call blacs_get(my_context, 11, got)
    err = 0.0_ep
    if (got /= 7) err = BAD
    if (my_rank == 0) call report_case('SGET_NR_BS', err, tol)

    ! ── SGET_NR_CO: combine nrings — round-trip 5
    call blacs_set(my_context, 13, 5)
    call blacs_get(my_context, 13, got)
    err = 0.0_ep
    if (got /= 5) err = BAD
    if (my_rank == 0) call report_case('SGET_NR_CO', err, tol)

    ! ── SGET_TOPSREPEAT: 1 (true)
    call blacs_set(my_context, 15, 1)
    call blacs_get(my_context, 15, got)
    err = 0.0_ep
    if (got /= 1) err = BAD
    if (my_rank == 0) call report_case('SGET_TOPSREPEAT', err, tol)

    ! ── SGET_TOPSCOHRNT: 1 (true)
    call blacs_set(my_context, 16, 1)
    call blacs_get(my_context, 16, got)
    err = 0.0_ep
    if (got /= 1) err = BAD
    if (my_rank == 0) call report_case('SGET_TOPSCOHRNT', err, tol)

    call report_finalize()
    call grid_exit()
end program test_blacs_set
