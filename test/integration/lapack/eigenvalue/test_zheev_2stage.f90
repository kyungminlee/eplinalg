! 2-stage Hermitian eigensolver. Only JOBZ='N' is supported.
program test_zheev_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zheev_2stage
    use ref_quad_lapack, only: zheev_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, lwork, lrwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work_ref(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork_ref(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zheev_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 81001 + 73 * i)
        lrwork = max(1, 3*n - 2)
        allocate(A_ref(n, n), A_got(n, n), w_ref(n), w_got(n), rwork_ref(lrwork))
        A_ref = A0;  A_got = A0
        call zheev_2stage('N', uplos(i), n, A_ref, n, w_ref, wopt, -1, &
                          rwork_ref, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work_ref(lwork))
        call zheev_2stage('N', uplos(i), n, A_ref, n, w_ref, work_ref, lwork, &
                          rwork_ref, info)
        call target_zheev_2stage('N', uplos(i), n, A_got, n, w_got, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,a,a,i0)') 'jobz=N,uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, w_ref, w_got, work_ref, rwork_ref)
    end do
    call report_finalize()
end program test_zheev_2stage
