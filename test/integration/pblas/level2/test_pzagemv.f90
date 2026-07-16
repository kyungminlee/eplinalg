! pzagemv: row/column 1-norm of |alpha|*|op(A)|*|X|, plus |beta*Y|.
!   sub(Y)_i := |alpha| * sum_j cabs1(op(A)_{i,j}) * cabs1(X_j) + |beta * sub(Y)_i|
! A and X are complex, alpha/beta/Y are real. The routine uses the
! Manhattan / L1 norm `cabs1(z) = |Re(z)| + |Im(z)|` (NOT the Euclidean
! magnitude — see ZAGEMV in PBLAS PTZBLAS), and conjugate transpose
! ('C') agrees with plain transpose ('T') under cabs1 since
! cabs1(conjg(z)) = cabs1(z).
program test_pzagemv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gen_distrib_vector, gather_vector
    use target_pblas,  only: target_name, target_eps, target_pzagemv
    implicit none

    integer, parameter :: ms(*) = [32, 80]
    integer, parameter :: ns(*) = [40, 60]
    character(len=1), parameter :: transes(*) = ['N', 'T', 'C']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, it, m, n, info, ii, jj
    integer :: lenx, leny, locm_a, locn_a, locn_x, locn_y, lld_a, lld_x, lld_y
    integer :: desca(9), descx(9), descy(9)
    character(len=1) :: trans
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), A_glob(:,:), x_glob(:)
    real(ep),    allocatable :: y_loc(:), y_glob(:), y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol, acc
    character(len=48) :: label

    ! cabs1: Manhattan norm of complex scalar — matches the routine's
    ! `CABS1` statement function in PTZBLAS/zagemv.f.
    !
    !   cabs1(z) = |Re(z)| + |Im(z)|

    call grid_init()
    call report_init('pzagemv', target_name, my_rank)

    alpha = 0.6_ep; beta = 0.3_ep
    do it = 1, size(transes)
        trans = transes(it)
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (trans == 'N') then
                lenx = n; leny = m
            else
                lenx = m; leny = n
            end if

            call gen_distrib_matrix_z(m, n, mb, nb, A_loc, A_glob, &
                                      seed = 7101 + 11 * i + 113 * it)
            call gen_distrib_vector_z(lenx, nb, x_loc, x_glob, &
                                      seed = 7121 + 11 * i + 113 * it)
            call gen_distrib_vector(leny, mb, y_loc, y_glob, &
                                    seed = 7141 + 11 * i + 113 * it)

            locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(n, nb, my_col, 0, my_npcol)
            lld_a  = max(1, locm_a)
            locn_x = numroc_local(lenx, nb, my_row, 0, my_nprow); lld_x = max(1, locn_x)
            locn_y = numroc_local(leny, mb, my_row, 0, my_nprow); lld_y = max(1, locn_y)

            call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)
            call descinit_local(descx, lenx, 1, nb, 1, 0, 0, my_context, lld_x, info)
            call descinit_local(descy, leny, 1, mb, 1, 0, 0, my_context, lld_y, info)

            call target_pzagemv(trans, m, n, alpha, A_loc, 1, 1, desca, &
                                x_loc, 1, 1, descx, 1, beta, &
                                y_loc, 1, 1, descy, 1)
            call gather_vector(leny, mb, y_loc, y_got)

            if (my_rank == 0) then
                allocate(y_ref(leny))
                if (trans == 'N') then
                    do ii = 1, m
                        acc = 0.0_ep
                        do jj = 1, n
                            acc = acc + (abs(real(A_glob(ii, jj), ep)) + abs(aimag(A_glob(ii, jj)))) * &
                                        (abs(real(x_glob(jj), ep)) + abs(aimag(x_glob(jj))))
                        end do
                        y_ref(ii) = abs(alpha) * acc + abs(beta * y_glob(ii))
                    end do
                else
                    ! TRANS='T' or 'C': op(A)_{i,j} = A_{j,i} (conjg
                    ! drops out under cabs1).
                    do jj = 1, n
                        acc = 0.0_ep
                        do ii = 1, m
                            acc = acc + (abs(real(A_glob(ii, jj), ep)) + abs(aimag(A_glob(ii, jj)))) * &
                                        (abs(real(x_glob(ii), ep)) + abs(aimag(x_glob(ii))))
                        end do
                        y_ref(jj) = abs(alpha) * acc + abs(beta * y_glob(jj))
                    end do
                end if
                err = max_rel_err_vec(y_got, y_ref)
                tol = 32.0_ep * 2.0_ep * real(max(m, n), ep) * target_eps
                write(label, '(a,a,a,i0,a,i0)') 'trans=', trans, ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(y_ref, y_got)
            end if
            deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzagemv
