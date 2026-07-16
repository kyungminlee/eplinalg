! dlamrg: merge two sorted runs into one index permutation.
program test_dlamrg
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use target_lapack,   only: target_name, target_dlamrg
    use ref_quad_lapack, only: dlamrg
    implicit none

    integer, parameter :: ns(*) = [4, 8, 16]
    integer :: i, n1, n2, j
    real(ep), allocatable :: A(:)
    integer, allocatable :: idx_r(:), idx_g(:)
    integer :: nfail
    character(len=48) :: label

    call report_init('dlamrg', target_name)
    do i = 1, size(ns)
        n1 = ns(i); n2 = ns(i)
        allocate(A(n1+n2), idx_r(n1+n2), idx_g(n1+n2))
        ! Two sorted ascending runs
        do j = 1, n1; A(j)    = real(j, ep);     end do
        do j = 1, n2; A(n1+j) = real(j, ep) + 0.5_ep; end do
        call dlamrg(n1, n2, A, 1, 1, idx_r)
        call target_dlamrg(n1, n2, A, 1, 1, idx_g)
        nfail = count(idx_r /= idx_g)
        write(label, '(a,i0)') 'n=', n1+n2
        ! integer comparison: 0 differences = pass
        if (nfail == 0) then
            call report_case(trim(label), 0.0_ep, 1.0_ep)
        else
            call report_case(trim(label), 1.0_ep, 0.0_ep)
        end if
        deallocate(A, idx_r, idx_g)
    end do
    call report_finalize()
end program test_dlamrg
