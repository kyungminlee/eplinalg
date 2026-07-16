! dtfsm: triangular solve where A is in RFP storage.
program test_dtfsm
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad, gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtfsm
    use ref_quad_lapack, only: dtrttf, dtfsm
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [4,  8]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    character(len=1), parameter :: sides(2)   = ['L', 'R']
    character(len=1), parameter :: transes(2) = ['N', 'T']
    character(len=1), parameter :: diags(2)   = ['N', 'U']
    integer :: i, m, n, info, nq, nt, t, u, s, tt, dd
    real(ep), allocatable :: A(:,:), Atri(:,:), B0(:,:), ARF(:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtfsm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        do s = 1, size(sides)
            nq = merge(m, n, sides(s) == 'L')
            nt = nq*(nq+1)/2
            call gen_spd_matrix_quad(nq, A, seed = 21601 + 61 * i + s)
            allocate(Atri(nq, nq))
            Atri = A
            call gen_matrix_quad(m, n, B0, seed = 21611 + 61 * i + s)
            do t = 1, size(transrs)
                do u = 1, size(uplos)
                    allocate(ARF(nt))
                    call dtrttf(transrs(t), uplos(u), nq, Atri, nq, ARF, info)
                    do tt = 1, size(transes)
                        do dd = 1, size(diags)
                            allocate(B_ref(m, n), B_got(m, n))
                            B_ref = B0; B_got = B0
                            call dtfsm(transrs(t), sides(s), uplos(u), transes(tt), &
                                       diags(dd), m, n, 1.0_ep, ARF, B_ref, m)
                            call target_dtfsm(transrs(t), sides(s), uplos(u), transes(tt), &
                                              diags(dd), m, n, 1.0_ep, ARF, B_got, m)
                            err = max_rel_err_mat(B_got, B_ref)
                            tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
                            write(label, '(a,a,a,a,a,a,a,a,a,a,a,i0)') &
                                'tR=', transrs(t), ',s=', sides(s), ',u=', uplos(u), &
                                ',t=', transes(tt), ',d=', diags(dd), ',m=', m
                            call report_case(trim(label), err, tol)
                            deallocate(B_ref, B_got)
                        end do
                    end do
                    deallocate(ARF)
                end do
            end do
            deallocate(Atri, A, B0)
        end do
    end do
    call report_finalize()
end program test_dtfsm
