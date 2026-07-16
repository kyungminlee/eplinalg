! dtrtri: triangular inverse. Build a well-conditioned triangular A by
! starting from a random matrix with diagonally-dominant magnitude.
program test_dtrtri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrtri
    use ref_quad_lapack, only: dtrtri
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    character(len=1), parameter :: diags(2) = ['N', 'U']
    integer :: i, n, info, ju, jd, j, k
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtrtri', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 24001 + 47 * i)
        do j = 1, n
            A0(j, j) = A0(j, j) + real(2*n, ep)  ! diagonal dominance
        end do
        do ju = 1, size(uplos)
            do jd = 1, size(diags)
                allocate(A_ref(n,n), A_got(n,n))
                A_ref = A0; A_got = A0
                call dtrtri(uplos(ju), diags(jd), n, A_ref, n, info)
                call target_dtrtri(uplos(ju), diags(jd), n, A_got, n, info)
                if (uplos(ju) == 'U') then
                    do j = 1, n
                        do k = j+1, n
                            A_ref(k, j) = 0.0_ep; A_got(k, j) = 0.0_ep
                        end do
                    end do
                else
                    do j = 1, n
                        do k = 1, j-1
                            A_ref(k, j) = 0.0_ep; A_got(k, j) = 0.0_ep
                        end do
                    end do
                end if
                err = max_rel_err_mat(A_got, A_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,a,a,i0)') 'uplo=', uplos(ju), &
                    ',diag=', diags(jd), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(A_ref, A_got)
            end do
        end do
    end do
    call report_finalize()
end program test_dtrtri
