program test_zgemv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zgemv
    use ref_quad_blas, only: zgemv
    implicit none

    integer, parameter :: ms(*) = [5, 50, 100]
    integer, parameter :: ns(*) = [7, 60, 120]
    ! Cycle TRANS — 'C' (conjugate-transpose) is its own code path
    ! and the most bug-prone in complex BLAS. For 'T'/'C' the y output
    ! has length n (not m), so dimensions swap.
    character(len=1), parameter :: transes(*) = ['N', 'C', 'T']
    integer :: i, m, n, len_x, len_y
    complex(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgemv', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        if (transes(i) == 'N') then
            len_x = n; len_y = m
        else
            len_x = m; len_y = n
        end if
        call gen_matrix_complex(m, n, A,  seed = 711 + 19 * i)
        call gen_vector_complex(len_x, x, seed = 721 + 19 * i)
        call gen_vector_complex(len_y, y0, seed = 731 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.3_ep, -0.1_ep, ep)
        allocate(y_ref(len_y), y_got(len_y))
        y_ref = y0
        y_got = y0
        call zgemv(transes(i), m, n, alpha, A, m, x, 1, beta, y_ref, 1)
        call target_zgemv(transes(i), m, n, alpha, A, m, x, 1, beta, &
                          y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * 4.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0,a,i0)') &
            'trans=', transes(i), ',m=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x, y0, y_ref, y_got)
    end do
    call report_finalize()
end program test_zgemv
