program test_pzgeequ
    ! Z-mirror of pdgeequ.
    use prec_kinds,       only: ep
    use compare,          only: rel_err_scalar, max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgeequ
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_vector, &
                                gather_vector_row
    use target_scalapack, only: target_name, target_eps, target_pzgeequ
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    real(ep), allocatable :: R_loc(:), C_loc(:), R_got(:), C_got(:)
    real(ep), allocatable :: R_ref(:), C_ref(:)
    real(ep) :: rowcnd, colcnd, amax
    real(ep) :: rowcnd_ref, colcnd_ref, amax_ref
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzgeequ', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 22101 + 31*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(R_loc(max(1, locm_a)), C_loc(max(1, locn_a)))
        call target_pzgeequ(n, n, A_loc, 1, 1, desca, R_loc, C_loc, &
                            rowcnd, colcnd, amax, info)

        call gather_vector(n, mb, R_loc, R_got)
        call gather_vector_row(n, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(R_ref(n), C_ref(n))
            call zgeequ(n, n, A_glob, n, R_ref, C_ref, &
                        rowcnd_ref, colcnd_ref, amax_ref, info_ref)
            err = max(max_rel_err_vec(R_got, R_ref), &
                      max_rel_err_vec(C_got, C_ref), &
                      rel_err_scalar(rowcnd, rowcnd_ref), &
                      rel_err_scalar(colcnd, colcnd_ref), &
                      rel_err_scalar(amax, amax_ref))
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(R_ref, C_ref, R_got, C_got)
        end if
        deallocate(A_loc, A_glob, R_loc, C_loc)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgeequ
