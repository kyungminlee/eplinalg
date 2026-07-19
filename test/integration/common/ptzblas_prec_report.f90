! Compatibility re-export: the PTZBLAS test drivers were written
! against a suite-renamed reporter module from the era when every test
! tree wrote its .mod files into one shared fmod directory. The
! implementation now lives in prec_report.F90 (module prec_report,
! compiled per suite with no variant macros — PTZBLAS is serial and
! its drivers pass rank 0 explicitly); this shim keeps the drivers'
! `use ptzblas_prec_report` lines working unchanged.
module ptzblas_prec_report
    use prec_report
    implicit none
    public
end module ptzblas_prec_report
