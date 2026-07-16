! pzahemv: row 1-norm of |alpha|*|A_herm|*|X|, plus |beta*Y|.
! A is Hermitian (only UPLO triangle stored); the mirror element is
! conjg(A_{j,i}) for the off-diagonal, but cabs1(conjg(z)) = cabs1(z)
! so the reference is symmetric. Routine uses CABS1 (= |Re|+|Im|),
! NOT Euclidean — see PHEMV / ZAHEMV in PBLAS PTZBLAS.
program test_pzahemv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gen_distrib_vector, gather_vector
    use target_pblas,  only: target_name, target_eps, target_pzahemv
    implicit none

    integer, parameter :: ns(*) = [40, 80]
    character(len=1), parameter :: uplos(*) = ['U', 'L']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, iu, n, info, ii, jj
    integer :: locm_a, locn_a, locn_x, locn_y, lld_a, lld_x, lld_y
    integer :: desca(9), descx(9), descy(9)
    character(len=1) :: uplo
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), A_glob(:,:), x_glob(:)
    real(ep),    allocatable :: y_loc(:), y_glob(:), y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol, acc, aij_abs
    character(len=48) :: label

    call grid_init()
    call report_init('pzahemv', target_name, my_rank)

    alpha = 0.6_ep; beta = 0.3_ep
    do iu = 1, size(uplos)
        uplo = uplos(iu)
        do i = 1, size(ns)
            n = ns(i)

            call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, &
                                      seed = 7201 + 11 * i + 113 * iu)
            call gen_distrib_vector_z(n, nb, x_loc, x_glob, &
                                      seed = 7221 + 11 * i + 113 * iu)
            call gen_distrib_vector(n, mb, y_loc, y_glob, &
                                    seed = 7241 + 11 * i + 113 * iu)

            locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(n, nb, my_col, 0, my_npcol)
            lld_a  = max(1, locm_a)
            locn_x = numroc_local(n, nb, my_row, 0, my_nprow); lld_x = max(1, locn_x)
            locn_y = numroc_local(n, mb, my_row, 0, my_nprow); lld_y = max(1, locn_y)

            call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)
            call descinit_local(descx, n, 1, nb, 1, 0, 0, my_context, lld_x, info)
            call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld_y, info)

            call target_pzahemv(uplo, n, alpha, A_loc, 1, 1, desca, &
                                x_loc, 1, 1, descx, 1, beta, &
                                y_loc, 1, 1, descy, 1)
            call gather_vector(n, mb, y_loc, y_got)

            if (my_rank == 0) then
                allocate(y_ref(n))
                ! Off-diagonal: cabs1(A_{i,j}) (for (i,j) in UPLO triangle)
                ! or cabs1(conjg(A_{j,i})) = cabs1(A_{j,i}) for the mirror.
                ! Diagonal: routine uses ABS(REAL(A(j,j))) — Hermitian
                ! invariant says diagonal is real, so the imaginary part
                ! is ignored even if the test feeds non-Hermitian noise.
                ! cabs1(z) = |Re(z)| + |Im(z)| (matches CABS1 in zahemv.f).
                do ii = 1, n
                    acc = 0.0_ep
                    do jj = 1, n
                        if (ii == jj) then
                            aij_abs = abs(real(A_glob(jj, jj), ep))   ! diagonal: real part only
                        else if (uplo == 'U') then
                            if (ii < jj) then
                                aij_abs = abs(real(A_glob(ii, jj), ep)) + abs(aimag(A_glob(ii, jj)))
                            else
                                aij_abs = abs(real(A_glob(jj, ii), ep)) + abs(aimag(A_glob(jj, ii)))
                            end if
                        else
                            if (ii > jj) then
                                aij_abs = abs(real(A_glob(ii, jj), ep)) + abs(aimag(A_glob(ii, jj)))
                            else
                                aij_abs = abs(real(A_glob(jj, ii), ep)) + abs(aimag(A_glob(jj, ii)))
                            end if
                        end if
                        acc = acc + aij_abs * (abs(real(x_glob(jj), ep)) + abs(aimag(x_glob(jj))))
                    end do
                    y_ref(ii) = abs(alpha) * acc + abs(beta * y_glob(ii))
                end do
                err = max_rel_err_vec(y_got, y_ref)
                tol = 32.0_ep * 2.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,i0)') 'uplo=', uplo, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(y_ref, y_got)
            end if
            deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzahemv
