! Compatibility re-export: the PBLAS test drivers were written against
! a suite-renamed copy of the refblas_quad interface module from the
! era when every test tree wrote its .mod files into one shared fmod
! directory. The interfaces now live in ref_quad_blas.f90 (module
! ref_quad_blas); this shim keeps the drivers'
! `use pblas_ref_quad_blas` lines working unchanged.
module pblas_ref_quad_blas
    use ref_quad_blas
    implicit none
    public
end module pblas_ref_quad_blas
