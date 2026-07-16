! Explicit interfaces to refblas_quad — vendored Netlib BLAS compiled
! at REAL(KIND=16). PTZBLAS tests use the BLAS analogues for vasum /
! vvdot / dzvasum which match dasum / ddot / dzasum byte-for-byte.
!
! Module name is ptzblas_ref_quad_blas (not ref_quad_blas) to avoid
! .mod collisions with tests/blas and tests/pblas equivalents — all
! targets write into the same ${PROJECT_BINARY_DIR}/fmod directory.
!
! refblas_quad's exports are renamed to *_quad_ (see
! tests/blas/refblas/CMakeLists.txt) so they do not collide with the
! migrated libblas's native-precision symbols. We declare the *_quad
! externals here and add thin module-procedure wrappers under the
! original names so test sources can keep `use ptzblas_ref_quad_blas,
! only: dasum` unchanged.
module ptzblas_ref_quad_blas
    use prec_kinds, only: ep
    implicit none

    interface
        function ddot_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,  intent(in) :: n, incx, incy
            real(ep), intent(in) :: x(*), y(*)
            real(ep) :: r
        end function ddot_quad

        function dasum_quad(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dasum_quad

        function zdotc_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotc_quad

        function zdotu_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotu_quad

        function dzasum_quad(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dzasum_quad
    end interface

contains

    function ddot(n, x, incx, y, incy) result(r)
        integer,  intent(in) :: n, incx, incy
        real(ep), intent(in) :: x(*), y(*)
        real(ep) :: r
        r = ddot_quad(n, x, incx, y, incy)
    end function ddot

    function dasum(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        r = dasum_quad(n, x, incx)
    end function dasum

    function zdotc(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = zdotc_quad(n, x, incx, y, incy)
    end function zdotc

    function zdotu(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = zdotu_quad(n, x, incx, y, incy)
    end function zdotu

    function dzasum(n, x, incx) result(r)
        integer,     intent(in) :: n, incx
        complex(ep), intent(in) :: x(*)
        real(ep) :: r
        r = dzasum_quad(n, x, incx)
    end function dzasum
end module ptzblas_ref_quad_blas
