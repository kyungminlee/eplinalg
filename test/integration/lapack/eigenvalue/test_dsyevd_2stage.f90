! 2-stage divide-and-conquer symmetric eigensolver. JOBZ='N' only.
program test_dsyevd_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsyevd_2stage
    use ref_quad_lapack, only: dsyevd_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, lwork, liwork, iwopt(1)
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: w_ref(:), w_got(:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsyevd_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 82001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
        A_ref = A0; A_got = A0
        call dsyevd_2stage('N', uplos(i), n, A_ref, n, w_ref, &
                           wopt, -1, iwopt, -1, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dsyevd_2stage('N', uplos(i), n, A_ref, n, w_ref, &
                           work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        call target_dsyevd_2stage('N', uplos(i), n, A_got, n, w_got, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,a,a,i0)') 'jobz=N,uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, w_ref, w_got)
    end do
    call report_finalize()
end program test_dsyevd_2stage
