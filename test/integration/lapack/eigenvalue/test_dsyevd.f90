program test_dsyevd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsyevd
    use ref_quad_lapack, only: dsyevd
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, jz, ju, lwork, liwork, iwopt(1)
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), w_ref(:), w_got(:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsyevd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 19001 + 47 * i)
        do jz = 1, size(jobzs)
            do ju = 1, size(uplos)
                allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
                A_ref = A0; A_got = A0
                call dsyevd(jobzs(jz), uplos(ju), n, A_ref, n, w_ref, &
                            wopt, -1, iwopt, -1, info)
                lwork  = max(1, int(wopt(1)))
                liwork = max(1, iwopt(1))
                allocate(work(lwork), iwork(liwork))
                call dsyevd(jobzs(jz), uplos(ju), n, A_ref, n, w_ref, &
                            work, lwork, iwork, liwork, info)
                deallocate(work, iwork)
                call target_dsyevd(jobzs(jz), uplos(ju), n, A_got, n, w_got, info)
                err = max_rel_err_vec(w_got, w_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,a,a,i0)') 'jobz=', jobzs(jz), ',uplo=', &
                    uplos(ju), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(A_ref, A_got, w_ref, w_got)
            end do
        end do
    end do
    call report_finalize()
end program test_dsyevd
