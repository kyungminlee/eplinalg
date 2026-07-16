! dorm22: apply 2x2 block-triangular orthogonal matrix Q to C.
program test_dorm22
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorm22
    use ref_quad_lapack, only: dorm22
    implicit none

    integer, parameter :: ns(*) = [10, 14]
    character(len=1), parameter :: sides(2) = ['L', 'R']
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, n, n1, n2, info, lwork, js, jt
    real(ep), allocatable :: Q(:,:), C0(:,:), Cr(:,:), Cg(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=64) :: label

    call report_init('dorm22', target_name)
    do i = 1, size(ns)
        n = ns(i); n1 = n/2; n2 = n - n1
        call gen_matrix_quad(n, n, Q, seed = 27001 + 79 * i)
        ! Zero off-diagonal blocks to enforce 2x2 block-triangular shape.
        Q(1:n1, n1+1:n) = 0.0_ep
        Q(n1+1:n, 1:n1) = 0.0_ep
        do js = 1, size(sides)
            do jt = 1, size(transes)
                call gen_matrix_quad(n, n, C0, seed = 27011 + 79*i + js + 3*jt)
                allocate(Cr(n, n), Cg(n, n))
                Cr = C0; Cg = C0
                call dorm22(sides(js), transes(jt), n, n, n1, n2, Q, n, Cr, n, &
                            wopt, -1, info)
                lwork = max(1, int(wopt(1))); allocate(work(lwork))
                call dorm22(sides(js), transes(jt), n, n, n1, n2, Q, n, Cr, n, &
                            work, lwork, info)
                deallocate(work)
                call target_dorm22(sides(js), transes(jt), n, n, n1, n2, Q, n, Cg, n, info)
                err = max_rel_err_mat(Cg, Cr)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,a,a,i0)') 'side=', sides(js), &
                    ',trans=', transes(jt), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(Cr, Cg, C0)
            end do
        end do
        deallocate(Q)
    end do
    call report_finalize()
end program test_dorm22
