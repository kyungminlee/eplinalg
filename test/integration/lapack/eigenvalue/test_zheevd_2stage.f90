! 2-stage divide-and-conquer Hermitian eigensolver. JOBZ='N' only.
program test_zheevd_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zheevd_2stage
    use ref_quad_lapack, only: zheevd_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, lwork, lrwork, liwork, iwopt(1)
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    integer, allocatable :: iwork(:)
    character(len=48) :: label

    call report_init('zheevd_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 83001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
        A_ref = A0; A_got = A0
        call zheevd_2stage('N', uplos(i), n, A_ref, n, w_ref, &
                           wopt, -1, rwopt, -1, iwopt, -1, info)
        lwork  = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        call zheevd_2stage('N', uplos(i), n, A_ref, n, w_ref, work, lwork, &
                           rwork, lrwork, iwork, liwork, info)
        deallocate(work, rwork, iwork)
        call target_zheevd_2stage('N', uplos(i), n, A_got, n, w_got, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,a,a,i0)') 'jobz=N,uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, w_ref, w_got)
    end do
    call report_finalize()
end program test_zheevd_2stage
