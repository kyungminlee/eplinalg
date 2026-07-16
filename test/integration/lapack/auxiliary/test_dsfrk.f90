! dsfrk: rank-k update of a symmetric matrix in RFP form.
!   C := alpha*A*A^T + beta*C  (TRANS='N')
!   C := alpha*A^T*A + beta*C  (TRANS='T')
program test_dsfrk
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad, gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsfrk
    use ref_quad_lapack, only: dtrttf, dsfrk
    implicit none

    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: ks(*) = [4, 12]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, t, u, j, n, k, info, nt
    real(ep), allocatable :: A(:,:), Cmat(:,:), C0(:), C_ref(:), C_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dsfrk', target_name)
    do i = 1, size(ns)
        n = ns(i); k = ks(i); nt = n*(n+1)/2
        call gen_spd_matrix_quad(n, Cmat, seed = 21501 + 67 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(C0(nt), C_ref(nt), C_got(nt))
                call dtrttf(transrs(t), uplos(u), n, Cmat, n, C0, info)
                do j = 1, size(transes)
                    if (transes(j) == 'N') then
                        call gen_matrix_quad(n, k, A, seed = 21511 + 47 * i + j)
                    else
                        call gen_matrix_quad(k, n, A, seed = 21511 + 47 * i + j)
                    end if
                    C_ref = C0; C_got = C0
                    call dsfrk(transrs(t), uplos(u), transes(j), n, k, &
                               0.5_ep, A, size(A, 1), 0.7_ep, C_ref)
                    call target_dsfrk(transrs(t), uplos(u), transes(j), n, k, &
                                      0.5_ep, A, size(A, 1), 0.7_ep, C_got)
                    err = max_rel_err_vec(C_got, C_ref)
                    tol = 16.0_ep * real(n + k, ep) * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0)') 'tR=', transrs(t), &
                        ',u=', uplos(u), ',tA=', transes(j), ',n=', n
                    call report_case(trim(label), err, tol)
                    deallocate(A)
                end do
                deallocate(C0, C_ref, C_got)
            end do
        end do
        deallocate(Cmat)
    end do
    call report_finalize()
end program test_dsfrk
