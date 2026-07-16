! Test PBZMATADD — complex local block matrix add.
! Same per-rank-local pattern as test_pbdmatadd, but with the 'C'
! conjugate-transpose mode exercised for complex data (where it
! actually does something — conjg + transpose).
program test_pbzmatadd
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_mat_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_matrix_complex
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context
    use target_pbblas,      only: target_name, target_eps, target_pbzmatadd
    use mpi
    implicit none

    integer :: m, n, i, j, ierr, fail_local, fail_global
    complex(ep), allocatable :: A(:,:), B0(:,:), B_got(:,:), B_ref(:,:)
    complex(ep) :: alpha, beta
    real(ep)    :: err, tol
    character(len=4), parameter :: modes(*) = ['U   ', 'L   ', 'V   ', 'T   ', 'C   ']
    character(len=64) :: label
    integer :: k

    call grid_init()
    call report_init('pbzmatadd', target_name, my_rank)

    alpha = (1.25_ep, 0.5_ep)
    beta  = (-0.75_ep, 0.25_ep)

    do k = 1, size(modes)
        m = 24; n = 32

        if (modes(k)(1:1) == 'T' .or. modes(k)(1:1) == 't' .or. &
            modes(k)(1:1) == 'C' .or. modes(k)(1:1) == 'c') then
            call gen_matrix_complex(n, m, A, seed = 1310 + 7 * k + 100 * my_rank)
        else
            call gen_matrix_complex(m, n, A, seed = 1310 + 7 * k + 100 * my_rank)
        end if
        call gen_matrix_complex(m, n, B0,    seed = 2410 + 7 * k + 100 * my_rank)
        allocate(B_got(m, n), B_ref(m, n))
        B_got = B0; B_ref = B0

        call target_pbzmatadd(my_context, modes(k)(1:1), m, n, alpha, &
                              A, size(A, 1), beta, B_got, m)

        select case (modes(k)(1:1))
        case ('U', 'u')
            do j = 1, n
                do i = 1, min(j, m)
                    B_ref(i, j) = alpha * A(i, j) + beta * B_ref(i, j)
                end do
            end do
        case ('L', 'l')
            do j = 1, n
                do i = j, m
                    B_ref(i, j) = alpha * A(i, j) + beta * B_ref(i, j)
                end do
            end do
        case ('T', 't')
            do j = 1, n
                do i = 1, m
                    B_ref(i, j) = alpha * A(j, i) + beta * B_ref(i, j)
                end do
            end do
        case ('C', 'c')
            do j = 1, n
                do i = 1, m
                    B_ref(i, j) = alpha * conjg(A(j, i)) + beta * B_ref(i, j)
                end do
            end do
        case default
            do j = 1, n
                do i = 1, m
                    B_ref(i, j) = alpha * A(i, j) + beta * B_ref(i, j)
                end do
            end do
        end select

        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 64.0_ep * real(m * n, ep) * target_eps

        fail_local = 0
        if (err > tol) fail_local = 1
        call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                           mpi_max, mpi_comm_world, ierr)

        if (my_rank == 0) then
            write(label, '(a,a,a,i0,a,i0,a,i0)') 'mode=', trim(modes(k)), &
                ',m=', m, ',n=', n, ',anyfail=', fail_global
            if (fail_global /= 0 .and. err <= tol) err = 2.0_ep * tol
            call report_case(trim(label), err, tol)
        end if
        deallocate(A, B0, B_got, B_ref)
    end do

    call report_finalize()
    call grid_exit()
end program test_pbzmatadd
