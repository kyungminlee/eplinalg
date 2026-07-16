program test_pzhemv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zhemv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzhemv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    character(len=1), parameter :: uplos(*) = ['U', 'L']
    integer, parameter :: mb = 8
    integer :: i, iu, n, info
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9), descy(9)
    character(len=1) :: uplo
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), &
                                y_got(:), y_ref(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzhemv', target_name, my_rank)

    alpha = cmplx(0.6_ep, 0.2_ep, ep); beta = cmplx(0.3_ep, -0.1_ep, ep)
    do iu = 1, size(uplos)
        uplo = uplos(iu)
        do i = 1, size(ns)
            n = ns(i)
            call gen_distrib_matrix_z(n, n, mb, mb, A_loc, A_glob, &
                                      seed = 5401 + 19 * i + 113 * iu)
            call gen_distrib_vector_z(n, mb, x_loc, x_glob, &
                                      seed = 5411 + 19 * i + 113 * iu)
            call gen_distrib_vector_z(n, mb, y_loc, y_glob, &
                                      seed = 5421 + 19 * i + 113 * iu)

            locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(n, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locn_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)

            call descinit_local(desca, n, n, mb, mb, 0, 0, my_context, lld_a, info)
            call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld_x, info)
            call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld_x, info)

            call target_pzhemv(uplo, n, alpha, A_loc, 1, 1, desca, &
                               x_loc, 1, 1, descx, 1, beta, &
                               y_loc, 1, 1, descy, 1)
            call gather_vector_z(n, mb, y_loc, y_got)

            if (my_rank == 0) then
                allocate(y_ref(n))
                y_ref = y_glob
                call zhemv(uplo, n, alpha, A_glob, n, x_glob, 1, beta, y_ref, 1)
                err = max_rel_err_vec_z(y_got, y_ref)
                tol = 32.0_ep * 8.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,i0)') 'uplo=', uplo, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(y_ref, y_got)
            end if
            deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzhemv
