! Test PBZVECADD — complex local vector add Y := alpha*op(X) + beta*Y.
! 'C' mode does Y := alpha * conjg(X) + beta*Y for complex data.
program test_pbzvecadd
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_vec_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_vector_complex
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context
    use target_pbblas,      only: target_name, target_eps, target_pbzvecadd
    use mpi
    implicit none

    integer, parameter :: cases(*) = [16, 100, 1000]
    character(len=2), parameter :: modes(*) = ['V ', 'C ']
    ! Stride pairs (incx, incy); see test_pbdvecadd for rationale.
    integer, parameter :: incx_list(*) = [1, 2, 3]
    integer, parameter :: incy_list(*) = [1, 3, 2]
    integer :: n, i, k, s, incx, incy, nx, ny, ix, iy, j, ierr
    integer :: fail_local, fail_global
    complex(ep), allocatable :: X(:), Y0(:), Y_got(:), Y_ref(:)
    complex(ep) :: alpha, beta
    real(ep)    :: err, tol
    character(len=80) :: label

    call grid_init()
    call report_init('pbzvecadd', target_name, my_rank)

    alpha = (1.7_ep, 0.3_ep)
    beta  = (-0.5_ep, 0.2_ep)

    do s = 1, size(incx_list)
        incx = incx_list(s)
        incy = incy_list(s)
        do k = 1, size(modes)
            do i = 1, size(cases)
                n = cases(i)
                nx = 1 + (n - 1) * abs(incx)
                ny = 1 + (n - 1) * abs(incy)
                call gen_vector_complex(nx, X,  seed = 911 + 13 * i + 100 * my_rank + 1000 * k + 7 * s)
                call gen_vector_complex(ny, Y0, seed = 1011 + 13 * i + 100 * my_rank + 1000 * k + 7 * s)
                allocate(Y_got(ny), Y_ref(ny))
                Y_got = Y0; Y_ref = Y0

                call target_pbzvecadd(my_context, modes(k)(1:1), n, alpha, X, incx, &
                                      beta, Y_got, incy)

                ix = 1; iy = 1
                if (modes(k)(1:1) == 'C' .or. modes(k)(1:1) == 'c') then
                    do j = 1, n
                        Y_ref(iy) = alpha * conjg(X(ix)) + beta * Y_ref(iy)
                        ix = ix + incx
                        iy = iy + incy
                    end do
                else
                    do j = 1, n
                        Y_ref(iy) = alpha * X(ix) + beta * Y_ref(iy)
                        ix = ix + incx
                        iy = iy + incy
                    end do
                end if

                err = max_rel_err_vec_z(Y_got, Y_ref)
                tol = 64.0_ep * real(n, ep) * target_eps

                fail_local = 0
                if (err > tol) fail_local = 1
                call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                                   mpi_max, mpi_comm_world, ierr)

                if (my_rank == 0) then
                    write(label, '(a,a,a,i0,a,i0,a,i0,a,i0)') 'mode=', trim(modes(k)), &
                        ',n=', n, ',incx=', incx, ',incy=', incy, ',anyfail=', fail_global
                    if (fail_global /= 0 .and. err <= tol) err = 2.0_ep * tol
                    call report_case(trim(label), err, tol)
                end if
                deallocate(X, Y0, Y_got, Y_ref)
            end do
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pbzvecadd
