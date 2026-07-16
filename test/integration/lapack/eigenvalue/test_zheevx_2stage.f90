! 2-stage selective Hermitian eigensolver. JOBZ='N', RANGE='A'.
program test_zheevx_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zheevx_2stage
    use ref_quad_lapack, only: zheevx_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, lwork, mref, mgot
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), z_ref(:,:), z_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    integer,     allocatable :: iwork(:), ifail(:)
    complex(ep) :: wopt(1)
    real(ep) :: abstol, vl, vu, err, tol
    character(len=48) :: label

    call report_init('zheevx_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 87001 + 47 * i)
        abstol = 0.0_ep; vl = 0.0_ep; vu = 0.0_ep
        allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
        allocate(z_ref(n,n), z_got(n,n), rwork(7*n), iwork(5*n), ifail(n))
        A_ref = A0; A_got = A0
        call zheevx_2stage('N', 'A', uplos(i), n, A_ref, n, vl, vu, 1, n, &
                           abstol, mref, w_ref, z_ref, n, &
                           wopt, -1, rwork, iwork, ifail, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zheevx_2stage('N', 'A', uplos(i), n, A_ref, n, vl, vu, 1, n, &
                           abstol, mref, w_ref, z_ref, n, &
                           work, lwork, rwork, iwork, ifail, info)
        deallocate(work)
        call target_zheevx_2stage('N', 'A', uplos(i), n, A_got, n, vl, vu, &
                                  1, n, abstol, mgot, w_got, z_got, n, &
                                  ifail, info)
        err = max_rel_err_vec(w_got(1:mgot), w_ref(1:mref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,a,a,i0)') 'jobz=N,uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, w_ref, w_got, z_ref, z_got, rwork, iwork, ifail)
    end do
    call report_finalize()
end program test_zheevx_2stage
