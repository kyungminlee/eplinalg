! Compatibility re-export: the PBLAS and ScaLAPACK test drivers were
! written against a suite-renamed reporter module from the era when
! every test tree wrote its .mod files into one shared fmod directory.
! The implementation now lives in prec_report.F90 (module prec_report,
! compiled per suite with PREC_REPORT_MPI); this shim keeps the
! drivers' `use pblas_prec_report` lines working unchanged.
module pblas_prec_report
    use prec_report
    implicit none
    public
end module pblas_prec_report
