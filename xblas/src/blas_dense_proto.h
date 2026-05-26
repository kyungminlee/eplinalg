#ifndef BLAS_DENSE_PROTO_H
#define BLAS_DENSE_PROTO_H

  /* Chapter 2 prototypes */

  /* Reduction Operations */

void BLAS_sdot( enum blas_conj_type conj, int n, float alpha,
                const float *x, int incx, float beta, const float *y,
                int incy, float *r );
void BLAS_edot( enum blas_conj_type conj, int n, EREAL alpha,
                const EREAL *x, int incx, EREAL beta, const EREAL *y,
                int incy, EREAL *r );
void BLAS_ddot( enum blas_conj_type conj, int n, double alpha,
                const double *x, int incx, double beta, const double *y,
                int incy, double *r );
void BLAS_edot( enum blas_conj_type conj, int n, EREAL alpha,
                const EREAL *x, int incx, EREAL beta, const EREAL *y,
                int incy, EREAL *r );
void BLAS_cdot( enum blas_conj_type conj, int n, const void *alpha,
                const void *x, int incx, const void *beta, const void *y,
                int incy, void *r );
void BLAS_ydot( enum blas_conj_type conj, int n, const void *alpha,
                const void *x, int incx, const void *beta, const void *y,
                int incy, void *r );
void BLAS_zdot( enum blas_conj_type conj, int n, const void *alpha,
                const void *x, int incx, const void *beta, const void *y,
                int incy, void *r );
void BLAS_ydot( enum blas_conj_type conj, int n, const void *alpha,
                const void *x, int incx, const void *beta, const void *y,
                int incy, void *r );

void BLAS_snorm( enum blas_norm_type norm, int n, const float *x, int incx,
                 float *r );
void BLAS_enorm( enum blas_norm_type norm, int n, const EREAL *x, int incx,
                 EREAL *r );
void BLAS_dnorm( enum blas_norm_type norm, int n, const double *x, int incx,
                 double *r );
void BLAS_enorm( enum blas_norm_type norm, int n, const EREAL *x, int incx,
                 EREAL *r );
void BLAS_cnorm( enum blas_norm_type norm, int n, const void *x, int incx,
                 float *r );
void BLAS_ynorm( enum blas_norm_type norm, int n, const void *x, int incx,
                 EREAL *r );
void BLAS_znorm( enum blas_norm_type norm, int n, const void *x, int incx,
                 double *r );
void BLAS_ynorm( enum blas_norm_type norm, int n, const void *x, int incx,
                 EREAL *r );

void BLAS_ssum( int n, const float *x, int incx, float *sum );
void BLAS_esum( int n, const EREAL *x, int incx, EREAL *sum );
void BLAS_dsum( int n, const double *x, int incx, double *sum );
void BLAS_esum( int n, const EREAL *x, int incx, EREAL *sum );
void BLAS_csum( int n, const void *x, int incx, void *sum );
void BLAS_ysum( int n, const void *x, int incx, void *sum );
void BLAS_zsum( int n, const void *x, int incx, void *sum );
void BLAS_ysum( int n, const void *x, int incx, void *sum );

void BLAS_smin_val( int n, const float *x, int incx, int k, float *r );
void BLAS_emin_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_dmin_val( int n, const double *x, int incx, int k, double *r );
void BLAS_emin_val( int n, const EREAL *x, int incx, int k, EREAL *r );

void BLAS_samin_val( int n, const float *x, int incx, int k, float *r );
void BLAS_eamin_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_damin_val( int n, const double *x, int incx, int k, double *r );
void BLAS_eamin_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_camin_val( int n, const void *x, int incx, int k, float *r );
void BLAS_yamin_val( int n, const void *x, int incx, int k, EREAL *r );
void BLAS_zamin_val( int n, const void *x, int incx, int k, double *r );
void BLAS_yamin_val( int n, const void *x, int incx, int k, EREAL *r );

void BLAS_smax_val( int n, const float *x, int incx, int k, float *r );
void BLAS_emax_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_dmax_val( int n, const double *x, int incx, int k, double *r );
void BLAS_emax_val( int n, const EREAL *x, int incx, int k, EREAL *r );

void BLAS_samax_val( int n, const float *x, int incx, int k, float *r );
void BLAS_eamax_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_damax_val( int n, const double *x, int incx, int k, double *r );
void BLAS_eamax_val( int n, const EREAL *x, int incx, int k, EREAL *r );
void BLAS_camax_val( int n, const void *x, int incx, int k, float *r );
void BLAS_yamax_val( int n, const void *x, int incx, int k, EREAL *r );
void BLAS_zamax_val( int n, const void *x, int incx, int k, double *r );
void BLAS_yamax_val( int n, const void *x, int incx, int k, EREAL *r );

void BLAS_ssumsq( int n, const float *x, int incx, float *ssq, float *scl );
void BLAS_esumsq( int n, const EREAL *x, int incx, EREAL *ssq, EREAL *scl );
void BLAS_dsumsq( int n, const double *x, int incx, double *ssq,
                  double *scl );
void BLAS_esumsq( int n, const EREAL *x, int incx, EREAL *ssq,
                  EREAL *scl );
void BLAS_csumsq( int n, const void *x, int incx, float *ssq, float *scl );
void BLAS_ysumsq( int n, const void *x, int incx, EREAL *ssq, EREAL *scl );
void BLAS_zsumsq( int n, const void *x, int incx, double *ssq, double *scl );
void BLAS_ysumsq( int n, const void *x, int incx, EREAL *ssq, EREAL *scl );

  /* Generate Transformations */

void BLAS_sgen_grot( float a, float b, float *c, float *s, float *r );
void BLAS_egen_grot( EREAL a, EREAL b, EREAL *c, EREAL *s, EREAL *r );
void BLAS_dgen_grot( double a, double b, double *c, double *s, double *r );
void BLAS_egen_grot( EREAL a, EREAL b, EREAL *c, EREAL *s, EREAL *r );
void BLAS_cgen_grot( const void *a, const void *b, float *c, void *s,
                     void *r );
void BLAS_ygen_grot( const void *a, const void *b, EREAL *c, void *s,
                     void *r );
void BLAS_zgen_grot( const void *a, const void *b, double *c, void *s,
                     void *r );
void BLAS_ygen_grot( const void *a, const void *b, EREAL *c, void *s,
                     void *r );

void BLAS_sgen_jrot( enum blas_jrot_type jrot, float *x, float y,
                     float *z, float *c, float *s );
void BLAS_egen_jrot( enum blas_jrot_type jrot, EREAL *x, EREAL y,
                     EREAL *z, EREAL *c, EREAL *s );
void BLAS_dgen_jrot( enum blas_jrot_type jrot, double *x, double y,
                     double *z, double *c, double *s );
void BLAS_egen_jrot( enum blas_jrot_type jrot, EREAL *x, EREAL y,
                     EREAL *z, EREAL *c, EREAL *s );
void BLAS_cgen_jrot( enum blas_jrot_type jrot, float *x, const void *y,
                     float *z, float *c, void *s );
void BLAS_ygen_jrot( enum blas_jrot_type jrot, EREAL *x, const void *y,
                     EREAL *z, EREAL *c, void *s );
void BLAS_zgen_jrot( enum blas_jrot_type jrot, double *x, const void *y,
                     double *z, double *c, void *s );
void BLAS_ygen_jrot( enum blas_jrot_type jrot, EREAL *x, const void *y,
                     EREAL *z, EREAL *c, void *s );

void BLAS_sgen_house( int n, float *alpha, float *x, int incx, float *tau );
void BLAS_egen_house( int n, EREAL *alpha, EREAL *x, int incx, EREAL *tau );
void BLAS_dgen_house( int n, double *alpha, double *x, int incx,
                      double *tau );
void BLAS_egen_house( int n, EREAL *alpha, EREAL *x, int incx,
                      EREAL *tau );
void BLAS_cgen_house( int n, void *alpha, void *x, int incx, void *tau );
void BLAS_ygen_house( int n, void *alpha, void *x, int incx, void *tau );
void BLAS_zgen_house( int n, void *alpha, void *x, int incx, void *tau );
void BLAS_ygen_house( int n, void *alpha, void *x, int incx, void *tau );

  /* Vector Operations */

void BLAS_srscale( int n, float alpha, float *x, int incx );
void BLAS_erscale( int n, EREAL alpha, EREAL *x, int incx );
void BLAS_drscale( int n, double alpha, double *x, int incx );
void BLAS_erscale( int n, EREAL alpha, EREAL *x, int incx );
void BLAS_crscale( int n, float alpha, void *x, int incx );
void BLAS_yrscale( int n, EREAL alpha, void *x, int incx );
void BLAS_zrscale( int n, double alpha, void *x, int incx );
void BLAS_yrscale( int n, EREAL alpha, void *x, int incx );

void BLAS_saxpby( int n, float alpha, const float *x, int incx,
                  float beta, float *y, int incy );
void BLAS_eaxpby( int n, EREAL alpha, const EREAL *x, int incx,
                  EREAL beta, EREAL *y, int incy );
void BLAS_daxpby( int n, double alpha, const double *x, int incx,
                  double beta, double *y, int incy );
void BLAS_eaxpby( int n, EREAL alpha, const EREAL *x, int incx,
                  EREAL beta, EREAL *y, int incy );
void BLAS_caxpby( int n, const void *alpha, const void *x, int incx,
                  const void *beta, void *y, int incy );
void BLAS_yaxpby( int n, const void *alpha, const void *x, int incx,
                  const void *beta, void *y, int incy );
void BLAS_zaxpby( int n, const void *alpha, const void *x, int incx,
                  const void *beta, void *y, int incy );
void BLAS_yaxpby( int n, const void *alpha, const void *x, int incx,
                  const void *beta, void *y, int incy );

void BLAS_swaxpby( int n, float alpha, const float *x, int incx,
                   float beta, const float *y, int incy, float *w,
                   int incw );
void BLAS_ewaxpby( int n, EREAL alpha, const EREAL *x, int incx,
                   EREAL beta, const EREAL *y, int incy, EREAL *w,
                   int incw );
void BLAS_dwaxpby( int n, double alpha, const double *x, int incx,
                   double beta, const double *y, int incy, double *w,
                   int incw );
void BLAS_ewaxpby( int n, EREAL alpha, const EREAL *x, int incx,
                   EREAL beta, const EREAL *y, int incy, EREAL *w,
                   int incw );
void BLAS_cwaxpby( int n, const void *alpha, const void *x, int incx,
                   const void *beta, const void *y, int incy, void *w,
                   int incw );
void BLAS_ywaxpby( int n, const void *alpha, const void *x, int incx,
                   const void *beta, const void *y, int incy, void *w,
                   int incw );
void BLAS_zwaxpby( int n, const void *alpha, const void *x, int incx,
                   const void *beta, const void *y, int incy, void *w,
                   int incw );
void BLAS_ywaxpby( int n, const void *alpha, const void *x, int incx,
                   const void *beta, const void *y, int incy, void *w,
                   int incw );

void BLAS_saxpy_dot( int n, float alpha, float *w, int incw, const float *v,
                     int incv, const float *u, int incu, float *r );
void BLAS_eaxpy_dot( int n, EREAL alpha, EREAL *w, int incw, const EREAL *v,
                     int incv, const EREAL *u, int incu, EREAL *r );
void BLAS_daxpy_dot( int n, double alpha, double *w, int incw,
                     const double *v, int incv, const double *u,
                     int incu, double *r );
void BLAS_eaxpy_dot( int n, EREAL alpha, EREAL *w, int incw,
                     const EREAL *v, int incv, const EREAL *u,
                     int incu, EREAL *r );
void BLAS_caxpy_dot( int n, const void *alpha, void *w, int incw,
                     const void *v, int incv, const void *u, int incu,
                     void *r );
void BLAS_yaxpy_dot( int n, const void *alpha, void *w, int incw,
                     const void *v, int incv, const void *u, int incu,
                     void *r );
void BLAS_zaxpy_dot( int n, const void *alpha, void *w, int incw,
                     const void *v, int incv, const void *u, int incu, 
                     void *r );
void BLAS_yaxpy_dot( int n, const void *alpha, void *w, int incw,
                     const void *v, int incv, const void *u, int incu, 
                     void *r );

void BLAS_sapply_grot( int n, float c, float s, float *x, int incx,
                       float *y, int incy );
void BLAS_eapply_grot( int n, EREAL c, EREAL s, EREAL *x, int incx,
                       EREAL *y, int incy );
void BLAS_dapply_grot( int n, double c, double s, double *x, int incx,
                       double *y, int incy );
void BLAS_eapply_grot( int n, EREAL c, EREAL s, EREAL *x, int incx,
                       EREAL *y, int incy );
void BLAS_capply_grot( int n, float c, const void *s, void *x, int incx,
                       void *y, int incy );
void BLAS_yapply_grot( int n, EREAL c, const void *s, void *x, int incx,
                       void *y, int incy );
void BLAS_zapply_grot( int n, double c, const void *s, void *x, int incx,
                       void *y, int incy );
void BLAS_yapply_grot( int n, EREAL c, const void *s, void *x, int incx,
                       void *y, int incy );

  /* Data Movement with Vectors */

void BLAS_scopy( int n, const float *x, int incx, float *y, int incy );
void BLAS_ecopy( int n, const EREAL *x, int incx, EREAL *y, int incy );
void BLAS_dcopy( int n, const double *x, int incx, double *y, int incy );
void BLAS_ecopy( int n, const EREAL *x, int incx, EREAL *y, int incy );
void BLAS_ccopy( int n, const void *x, int incx, void *y, int incy );
void BLAS_ycopy( int n, const void *x, int incx, void *y, int incy );
void BLAS_zcopy( int n, const void *x, int incx, void *y, int incy );
void BLAS_ycopy( int n, const void *x, int incx, void *y, int incy );

void BLAS_sswap( int n, float *x, int incx, float *y, int incy );
void BLAS_eswap( int n, EREAL *x, int incx, EREAL *y, int incy );
void BLAS_dswap( int n, double *x, int incx, double *y, int incy );
void BLAS_eswap( int n, EREAL *x, int incx, EREAL *y, int incy );
void BLAS_cswap( int n, void *x, int incx, void *y, int incy );
void BLAS_yswap( int n, void *x, int incx, void *y, int incy );
void BLAS_zswap( int n, void *x, int incx, void *y, int incy );
void BLAS_yswap( int n, void *x, int incx, void *y, int incy );

void BLAS_ssort( enum blas_sort_type sort, int n, float *x, int incx );
void BLAS_esort( enum blas_sort_type sort, int n, EREAL *x, int incx );
void BLAS_dsort( enum blas_sort_type sort, int n, double *x, int incx );
void BLAS_esort( enum blas_sort_type sort, int n, EREAL *x, int incx );

void BLAS_ssortv( enum blas_sort_type sort, int n, float *x, int incx,
                  int *p, int incp );
void BLAS_esortv( enum blas_sort_type sort, int n, EREAL *x, int incx,
                  int *p, int incp );
void BLAS_dsortv( enum blas_sort_type sort, int n, double *x, int incx,
                  int *p, int incp );
void BLAS_esortv( enum blas_sort_type sort, int n, EREAL *x, int incx,
                  int *p, int incp );

void BLAS_spermute( int n, const int *p, int incp, float *x, int incx );
void BLAS_epermute( int n, const int *p, int incp, EREAL *x, int incx );
void BLAS_dpermute( int n, const int *p, int incp, double *x, int incx );
void BLAS_epermute( int n, const int *p, int incp, EREAL *x, int incx );
void BLAS_cpermute( int n, const int *p, int incp, void *x, int incx );
void BLAS_ypermute( int n, const int *p, int incp, void *x, int incx );
void BLAS_zpermute( int n, const int *p, int incp, void *x, int incx );
void BLAS_ypermute( int n, const int *p, int incp, void *x, int incx );

  /* Matrix-Vector Operations */

void BLAS_sgemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, float alpha, const float *a, int lda,
                 const float *x, int incx, float beta, float *y, int incy );
void BLAS_egemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y, int incy );
void BLAS_dgemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, double alpha, const double *a, int lda,
                 const double *x, int incx, double beta, double *y,
                 int incy );
void BLAS_egemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y,
                 int incy );
void BLAS_cgemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ygemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_zgemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ygemv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );

void BLAS_sgbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, float alpha, const float *a,
                 int lda, const float *x, int incx, float beta,
                 float *y, int incy );
void BLAS_egbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, EREAL alpha, const EREAL *a,
                 int lda, const EREAL *x, int incx, EREAL beta,
                 EREAL *y, int incy );
void BLAS_dgbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, double alpha, const double *a,
                 int lda, const double *x, int incx, double beta,
                 double *y, int incy );
void BLAS_egbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, EREAL alpha, const EREAL *a,
                 int lda, const EREAL *x, int incx, EREAL beta,
                 EREAL *y, int incy );
void BLAS_cgbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, const void *alpha,
                 const void *a, int lda, const void *x, int incx,
                 const void *beta, void *y, int incy );
void BLAS_ygbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, const void *alpha,
                 const void *a, int lda, const void *x, int incx,
                 const void *beta, void *y, int incy );
void BLAS_zgbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, const void *alpha,
                 const void *a, int lda, const void *x, int incx,
                 const void *beta, void *y, int incy );
void BLAS_ygbmv( enum blas_order_type order, enum blas_trans_type trans,
                 int m, int n, int kl, int ku, const void *alpha,
                 const void *a, int lda, const void *x, int incx,
                 const void *beta, void *y, int incy );

void BLAS_ssymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, float alpha, const float *a, int lda,
                 const float *x, int incx, float beta, float *y, int incy );
void BLAS_esymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y, int incy );
void BLAS_dsymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, double alpha, const double *a, int lda,
                 const double *x, int incx, double beta, double *y,
                 int incy );
void BLAS_esymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y,
                 int incy );
void BLAS_csymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ysymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_zsymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ysymv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );

void BLAS_ssbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, float alpha, const float *a, int lda,
                 const float *x, int incx, float beta, float *y, int incy );
void BLAS_esbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y, int incy );
void BLAS_dsbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, double alpha, const double *a, int lda,
                 const double *x, int incx, double beta, double *y,
                 int incy );
void BLAS_esbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, EREAL alpha, const EREAL *a, int lda,
                 const EREAL *x, int incx, EREAL beta, EREAL *y,
                 int incy );
void BLAS_csbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ysbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_zsbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_ysbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );

void BLAS_sspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, float alpha, const float *ap, const float *x,
                 int incx, float beta, float *y, int incy );
void BLAS_espmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, EREAL alpha, const EREAL *ap, const EREAL *x,
                 int incx, EREAL beta, EREAL *y, int incy );
void BLAS_dspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, double alpha, const double *ap, const double *x,
                 int incx, double beta, double *y, int incy );
void BLAS_espmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, EREAL alpha, const EREAL *ap, const EREAL *x,
                 int incx, EREAL beta, EREAL *y, int incy );
void BLAS_cspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_yspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_zspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_yspmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );

void BLAS_chemv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_yhemv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_zhemv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_yhemv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );

void BLAS_chbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_yhbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_zhbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );
void BLAS_yhbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, int k, const void *alpha, const void *a, int lda,
                 const void *x, int incx, const void *beta, void *y,
                 int incy );

void BLAS_chpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_yhpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_zhpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );
void BLAS_yhpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 int n, const void *alpha, const void *ap, const void *x,
                 int incx, const void *beta, void *y, int incy );

void BLAS_strmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, float alpha, const float *t, int ldt, float *x,
                 int incx );
void BLAS_etrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *t, int ldt, EREAL *x,
                 int incx );
void BLAS_dtrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, double alpha, const double *t, int ldt, double *x,
                 int incx );
void BLAS_etrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *t, int ldt, EREAL *x,
                 int incx );
void BLAS_ctrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt, void *x,
                 int incx );
void BLAS_ytrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt, void *x,
                 int incx );
void BLAS_ztrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt, void *x,
                 int incx );
void BLAS_ytrmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt, void *x,
                 int incx );

void BLAS_stbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, float alpha, const float *t, int ldt,
                 float *x, int incx );
void BLAS_etbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_dtbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, double alpha, const double *t, int ldt,
                 double *x, int incx );
void BLAS_etbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_ctbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ztbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytbmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );

void BLAS_stpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, float alpha, const float *tp, float *x, int incx );
void BLAS_etpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *tp, EREAL *x, int incx );
void BLAS_dtpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, double alpha, const double *tp, double *x,
                 int incx );
void BLAS_etpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *tp, EREAL *x,
                 int incx );
void BLAS_ctpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x, 
                 int incx );
void BLAS_ytpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x, 
                 int incx );
void BLAS_ztpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );
void BLAS_ytpmv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );

void BLAS_sge_sum_mv( enum blas_order_type order, int m, int n,
                      float alpha, const float *a, int lda, const float *x,
                      int incx, float beta, const float *B, int ldb,
                      float *y, int incy );
void BLAS_ege_sum_mv( enum blas_order_type order, int m, int n,
                      EREAL alpha, const EREAL *a, int lda, const EREAL *x,
                      int incx, EREAL beta, const EREAL *B, int ldb,
                      EREAL *y, int incy );
void BLAS_dge_sum_mv( enum blas_order_type order, int m, int n,
                      double alpha, const double *a, int lda,
                      const double *x, int incx, double beta,
                      const double *B, int ldb, double *y, int incy );
void BLAS_ege_sum_mv( enum blas_order_type order, int m, int n,
                      EREAL alpha, const EREAL *a, int lda,
                      const EREAL *x, int incx, EREAL beta,
                      const EREAL *B, int ldb, EREAL *y, int incy );
void BLAS_cge_sum_mv( enum blas_order_type order, int m, int n,
                      const void *alpha, const void *a, int lda,
                      const void *x, int incx, const void *beta,
                      const void *B, int ldb, void *y, int incy );
void BLAS_yge_sum_mv( enum blas_order_type order, int m, int n,
                      const void *alpha, const void *a, int lda,
                      const void *x, int incx, const void *beta,
                      const void *B, int ldb, void *y, int incy );
void BLAS_zge_sum_mv( enum blas_order_type order, int m, int n,
                      const void *alpha, const void *a, int lda,
                      const void *x, int incx, const void *beta,
                      const void *B, int ldb, void *y, int incy );
void BLAS_yge_sum_mv( enum blas_order_type order, int m, int n,
                      const void *alpha, const void *a, int lda,
                      const void *x, int incx, const void *beta,
                      const void *B, int ldb, void *y, int incy );

void BLAS_sgemvt( enum blas_order_type order, int m, int n, float alpha,
                  const float *a, int lda, float *x, int incx,
                  const float *y, int incy, float beta, float *w, int incw,
                  const float *z, int incz );
void BLAS_egemvt( enum blas_order_type order, int m, int n, EREAL alpha,
                  const EREAL *a, int lda, EREAL *x, int incx,
                  const EREAL *y, int incy, EREAL beta, EREAL *w, int incw,
                  const EREAL *z, int incz );
void BLAS_dgemvt( enum blas_order_type order, int m, int n, double alpha,
                  const double *a, int lda, double *x, int incx,
                  const double *y, int incy, double beta, double *w,
                  int incw, const double *z, int incz );
void BLAS_egemvt( enum blas_order_type order, int m, int n, EREAL alpha,
                  const EREAL *a, int lda, EREAL *x, int incx,
                  const EREAL *y, int incy, EREAL beta, EREAL *w,
                  int incw, const EREAL *z, int incz );
void BLAS_cgemvt( enum blas_order_type order, int m, int n,
                  const void *alpha, const void *a, int lda, void *x,
                  int incx, const void *y, int incy, const void *beta,
                  void *w, int incw, const void *z, int incz );
void BLAS_ygemvt( enum blas_order_type order, int m, int n,
                  const void *alpha, const void *a, int lda, void *x,
                  int incx, const void *y, int incy, const void *beta,
                  void *w, int incw, const void *z, int incz );
void BLAS_zgemvt( enum blas_order_type order, int m, int n,
                  const void *alpha, const void *a, int lda, void *x,
                  int incx, const void *y, int incy, const void *beta,
                  void *w, int incw, const void *z, int incz );
void BLAS_ygemvt( enum blas_order_type order, int m, int n,
                  const void *alpha, const void *a, int lda, void *x,
                  int incx, const void *y, int incy, const void *beta,
                  void *w, int incw, const void *z, int incz );

void BLAS_strmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const float *t, int ldt, float *x, int incx,
                  const float *y, int incy, float *w, int incw,
                  const float *z, int incz );
void BLAS_etrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const EREAL *t, int ldt, EREAL *x, int incx,
                  const EREAL *y, int incy, EREAL *w, int incw,
                  const EREAL *z, int incz );
void BLAS_dtrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const float *t, int ldt, float *x, int incx,
                  const float *y, int incy, float *w, int incw,
                  const float *z, int incz );
void BLAS_etrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const EREAL *t, int ldt, EREAL *x, int incx,
                  const EREAL *y, int incy, EREAL *w, int incw,
                  const EREAL *z, int incz );
void BLAS_ctrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const void *t, int ldt, void *x, int incx,
                  const void *y, int incy, void *w, int incw,
                  const void *z, int incz );
void BLAS_ytrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const void *t, int ldt, void *x, int incx,
                  const void *y, int incy, void *w, int incw,
                  const void *z, int incz );
void BLAS_ztrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const void *t, int ldt, void *x, int incx,
                  const void *y, int incy, void *w, int incw,
                  const void *z, int incz );
void BLAS_ytrmvt( enum blas_order_type order, enum blas_uplo_type uplo,
                  int n, const void *t, int ldt, void *x, int incx,
                  const void *y, int incy, void *w, int incw,
                  const void *z, int incz );

void BLAS_sgemver( enum blas_order_type order, int m, int n, float *a,
                   int lda, const float *u1, const float *v1,
                   const float *u2, const float *v2, float alpha, float *x,
                   int incx, const float *y, int incy, float *w, int incw,
                   float beta, const float *z, int incz );
void BLAS_egemver( enum blas_order_type order, int m, int n, EREAL *a,
                   int lda, const EREAL *u1, const EREAL *v1,
                   const EREAL *u2, const EREAL *v2, EREAL alpha, EREAL *x,
                   int incx, const EREAL *y, int incy, EREAL *w, int incw,
                   EREAL beta, const EREAL *z, int incz );
void BLAS_dgemver( enum blas_order_type order, int m, int n, double *a,
                   int lda, const double *u1, const double *v1,
                   const double *u2, const double *v2, double alpha,
                   double *x, int incx, const double *y, int incy,
                   double *w, int incw, double beta, const double *z,
                   int incz );
void BLAS_egemver( enum blas_order_type order, int m, int n, EREAL *a,
                   int lda, const EREAL *u1, const EREAL *v1,
                   const EREAL *u2, const EREAL *v2, EREAL alpha,
                   EREAL *x, int incx, const EREAL *y, int incy,
                   EREAL *w, int incw, EREAL beta, const EREAL *z,
                   int incz );
void BLAS_cgemver( enum blas_order_type order, int m, int n, void *a,
                   int lda, const void *u1, const void *v1,
                   const void *u2, const void *v2, const void *alpha,
                   void *x, int incx, const void *y, int incy, void *w,
                   int incw, const void *beta, const void *z, int incz );
void BLAS_ygemver( enum blas_order_type order, int m, int n, void *a,
                   int lda, const void *u1, const void *v1,
                   const void *u2, const void *v2, const void *alpha,
                   void *x, int incx, const void *y, int incy, void *w,
                   int incw, const void *beta, const void *z, int incz );
void BLAS_zgemver( enum blas_order_type order, int m, int n, void *a,
                   int lda, const void *u1, const void *v1,
                   const void *u2, const void *v2, const void *alpha,
                   void *x, int incx, const void *y, int incy, void *w,
                   int incw, const void *beta, const void *z, int incz );
void BLAS_ygemver( enum blas_order_type order, int m, int n, void *a,
                   int lda, const void *u1, const void *v1,
                   const void *u2, const void *v2, const void *alpha,
                   void *x, int incx, const void *y, int incy, void *w,
                   int incw, const void *beta, const void *z, int incz );

void BLAS_strsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, float alpha, const float *t, int ldt,
                 float *x, int incx );
void BLAS_etrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_dtrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, double alpha, const double *t, int ldt,
                 double *x, int incx );
void BLAS_etrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_ctrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ztrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytrsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *t, int ldt,
                 void *x, int incx );

void BLAS_stbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, float alpha, const float *t, int ldt,
                 float *x, int incx );
void BLAS_etbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_dtbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, double alpha, const double *t, int ldt,
                 double *x, int incx );
void BLAS_etbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, EREAL alpha, const EREAL *t, int ldt,
                 EREAL *x, int incx );
void BLAS_ctbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ztbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );
void BLAS_ytbsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, int k, const void *alpha, const void *t, int ldt,
                 void *x, int incx );

void BLAS_stpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, float alpha, const float *tp, float *x, int incx );
void BLAS_etpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *tp, EREAL *x, int incx );
void BLAS_dtpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, double alpha, const double *tp, double *x,
                 int incx );
void BLAS_etpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, EREAL alpha, const EREAL *tp, EREAL *x,
                 int incx );
void BLAS_ctpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );
void BLAS_ytpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );
void BLAS_ztpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );
void BLAS_ytpsv( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, enum blas_diag_type diag,
                 int n, const void *alpha, const void *tp, void *x,
                 int incx );

void BLAS_sger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, float alpha, const float *x, int incx,
                const float *y, int incy, float beta, float *a, int lda );
void BLAS_eger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, EREAL alpha, const EREAL *x, int incx,
                const EREAL *y, int incy, EREAL beta, EREAL *a, int lda );
void BLAS_dger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, double alpha, const double *x, int incx,
                const double *y, int incy, double beta, double *a, int lda );
void BLAS_eger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, EREAL alpha, const EREAL *x, int incx,
                const EREAL *y, int incy, EREAL beta, EREAL *a, int lda );
void BLAS_cger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, const void *alpha, const void *x, int incx,
                const void *y, int incy, const void *beta, void *a,
                int lda );
void BLAS_yger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, const void *alpha, const void *x, int incx,
                const void *y, int incy, const void *beta, void *a,
                int lda );
void BLAS_zger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, const void *alpha, const void *x, int incx,
                const void *y, int incy, const void *beta, void *a,
                int lda );
void BLAS_yger( enum blas_order_type order, enum blas_conj_type conj,
                int m, int n, const void *alpha, const void *x, int incx,
                const void *y, int incy, const void *beta, void *a,
                int lda );

void BLAS_ssyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                float alpha, const float *x, int incx, float beta,
                float *a, int lda );
void BLAS_esyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const EREAL *x, int incx, EREAL beta,
                EREAL *a, int lda );
void BLAS_dsyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                double alpha, const double *x, int incx, double beta,
                double *a, int lda );
void BLAS_esyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const EREAL *x, int incx, EREAL beta,
                EREAL *a, int lda );
void BLAS_csyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *a, int lda );
void BLAS_ysyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *a, int lda );
void BLAS_zsyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *a, int lda );
void BLAS_ysyr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *a, int lda );

void BLAS_sspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                float alpha, const float *x, int incx, float beta,
                float *ap );
void BLAS_espr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const EREAL *x, int incx, EREAL beta,
                EREAL *ap );
void BLAS_dspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                double alpha, const double *x, int incx, double beta,
                double *ap );
void BLAS_espr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const EREAL *x, int incx, EREAL beta,
                EREAL *ap );
void BLAS_cspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *ap );
void BLAS_yspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *ap );
void BLAS_zspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *ap );
void BLAS_yspr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                const void *alpha, const void *x, int incx, const void *beta,
                void *ap );

void BLAS_cher( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                float alpha, const void *x, int incx, float beta,
                void *a, int lda );
void BLAS_yher( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const void *x, int incx, EREAL beta,
                void *a, int lda );
void BLAS_zher( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                double alpha, const void *x, int incx, double beta,
                void *a, int lda );
void BLAS_yher( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const void *x, int incx, EREAL beta,
                void *a, int lda );

void BLAS_chpr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                float alpha, const void *x, int incx, float beta,
                void *ap );
void BLAS_yhpr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const void *x, int incx, EREAL beta,
                void *ap );
void BLAS_zhpr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                double alpha, const void *x, int incx, double beta,
                void *ap );
void BLAS_yhpr( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                EREAL alpha, const void *x, int incx, EREAL beta,
                void *ap );

void BLAS_ssyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 float alpha, const float *x, int incx, const float *y,
                 int incy, float beta, float *a, int lda );
void BLAS_esyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 EREAL alpha, const EREAL *x, int incx, const EREAL *y,
                 int incy, EREAL beta, EREAL *a, int lda );
void BLAS_dsyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 double alpha, const double *x, int incx, const double *y,
                 int incy, double beta, double *a, int lda );
void BLAS_esyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 EREAL alpha, const EREAL *x, int incx, const EREAL *y,
                 int incy, EREAL beta, EREAL *a, int lda );
void BLAS_csyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *a, int lda );
void BLAS_ysyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *a, int lda );
void BLAS_zsyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *a, int lda );
void BLAS_ysyr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *a, int lda );

void BLAS_sspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 float alpha, const float *x, int incx, const float *y,
                 int incy, float beta, float *ap );
void BLAS_espr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 EREAL alpha, const EREAL *x, int incx, const EREAL *y,
                 int incy, EREAL beta, EREAL *ap );
void BLAS_dspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 double alpha, const double *x, int incx, const double *y,
                 int incy, double beta, double *ap );
void BLAS_espr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 EREAL alpha, const EREAL *x, int incx, const EREAL *y,
                 int incy, EREAL beta, EREAL *ap );
void BLAS_cspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *ap );
void BLAS_yspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *ap );
void BLAS_zspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *ap );
void BLAS_yspr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, const void *beta, void *ap );

void BLAS_cher2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, float beta, void *a, int lda );
void BLAS_yher2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, EREAL beta, void *a, int lda );
void BLAS_zher2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, double beta, void *a, int lda );
void BLAS_yher2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, EREAL beta, void *a, int lda );

void BLAS_chpr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, float beta, void *ap );
void BLAS_yhpr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, EREAL beta, void *ap );
void BLAS_zhpr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, double beta, void *ap );
void BLAS_yhpr2( enum blas_order_type order, enum blas_uplo_type uplo, int n,
                 const void *alpha, const void *x, int incx, const void *y,
                 int incy, EREAL beta, void *ap );

  /* Matrix Operations */

void BLAS_sge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const float *a, int lda, float *r );
void BLAS_ege_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const EREAL *a, int lda, EREAL *r );
void BLAS_dge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const double *a, int lda, double *r );
void BLAS_ege_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const EREAL *a, int lda, EREAL *r );
void BLAS_cge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const void *a, int lda, float *r );
void BLAS_yge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const void *a, int lda, EREAL *r );
void BLAS_zge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const void *a, int lda, double *r );
void BLAS_yge_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, const void *a, int lda, EREAL *r );

void BLAS_sgb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const float *a, int lda,
                    float *r );
void BLAS_egb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const EREAL *a, int lda,
                    EREAL *r );
void BLAS_dgb_norm( enum blas_order_type order, enum blas_norm_type norm, 
                    int m, int n, int kl, int ku, const double *a, int lda,
                    double *r );
void BLAS_egb_norm( enum blas_order_type order, enum blas_norm_type norm, 
                    int m, int n, int kl, int ku, const EREAL *a, int lda,
                    EREAL *r );
void BLAS_cgb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    float *r );
void BLAS_ygb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    EREAL *r );
void BLAS_zgb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    double *r );
void BLAS_ygb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    EREAL *r );

void BLAS_ssy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const float *a,
                    int lda, float *r );
void BLAS_esy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const EREAL *a,
                    int lda, EREAL *r );
void BLAS_dsy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const double *a,
                    int lda, double *r );
void BLAS_esy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const EREAL *a,
                    int lda, EREAL *r );
void BLAS_csy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a,
                    int lda, float *r );
void BLAS_ysy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a,
                    int lda, EREAL *r );
void BLAS_zsy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    double *r );
void BLAS_ysy_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    EREAL *r );

void BLAS_che_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    float *r );
void BLAS_yhe_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    EREAL *r );
void BLAS_zhe_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    double *r );
void BLAS_yhe_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *a, int lda,
                    EREAL *r );

void BLAS_ssb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const float *a,
                    int lda, float *r );
void BLAS_esb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const EREAL *a,
                    int lda, EREAL *r );
void BLAS_dsb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const double *a,
                    int lda, double *r );
void BLAS_esb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const EREAL *a,
                    int lda, EREAL *r );
void BLAS_csb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, float *r );
void BLAS_ysb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, EREAL *r );
void BLAS_zsb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, double *r );
void BLAS_ysb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, EREAL *r );

void BLAS_chb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, float *r );
void BLAS_yhb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, EREAL *r );
void BLAS_zhb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, double *r );
void BLAS_yhb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, int k, const void *a,
                    int lda, EREAL *r );

void BLAS_ssp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const float *ap,
                    float *r );
void BLAS_esp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const EREAL *ap,
                    EREAL *r );
void BLAS_dsp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const double *ap,
                    double *r );
void BLAS_esp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const EREAL *ap,
                    EREAL *r );
void BLAS_csp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    float *r );
void BLAS_ysp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    EREAL *r );
void BLAS_zsp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    double *r );
void BLAS_ysp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    EREAL *r );

void BLAS_chp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    float *r );
void BLAS_yhp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    EREAL *r );
void BLAS_zhp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    double *r );
void BLAS_yhp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, int n, const void *ap,
                    EREAL *r );

void BLAS_str_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const float *a, int lda, float *r );
void BLAS_etr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const EREAL *a, int lda, EREAL *r );
void BLAS_dtr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const double *a, int lda, double *r );
void BLAS_etr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const EREAL *a, int lda, EREAL *r );
void BLAS_ctr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *a, int lda, float *r );
void BLAS_ytr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *a, int lda, EREAL *r );
void BLAS_ztr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *a, int lda, double *r );
void BLAS_ytr_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *a, int lda, EREAL *r );

void BLAS_stb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const float *a, int lda, float *r );
void BLAS_etb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const EREAL *a, int lda, EREAL *r );
void BLAS_dtb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const double *a, int lda, double *r );
void BLAS_etb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const EREAL *a, int lda, EREAL *r );
void BLAS_ctb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, float *r );
void BLAS_ytb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, EREAL *r );
void BLAS_ztb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, double *r );
void BLAS_ytb_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, EREAL *r );

void BLAS_stp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const float *ap, float *r );
void BLAS_etp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const EREAL *ap, EREAL *r );
void BLAS_dtp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const double *ap, double *r );
void BLAS_etp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const EREAL *ap, EREAL *r );
void BLAS_ctp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *ap, float *r );
void BLAS_ytp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *ap, EREAL *r );
void BLAS_ztp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *ap, double *r );
void BLAS_ytp_norm( enum blas_order_type order, enum blas_norm_type norm,
                    enum blas_uplo_type uplo, enum blas_diag_type diag,
                    int n, const void *ap, EREAL *r );

void BLAS_sge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const float *d, int incd, float *a, int lda );
void BLAS_ege_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const EREAL *d, int incd, EREAL *a, int lda );
void BLAS_dge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const double *d, int incd, double *a, int lda );
void BLAS_ege_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const EREAL *d, int incd, EREAL *a, int lda );
void BLAS_cge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const void *d, int incd, void *a, int lda );
void BLAS_yge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          const void *d, int incd, void *a, int lda );
void BLAS_zge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n, 
                          const void *d, int incd, void *a, int lda );
void BLAS_yge_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n, 
                          const void *d, int incd, void *a, int lda );

void BLAS_sgb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const float *d, int incd,
                          float *a, int lda );
void BLAS_egb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const EREAL *d, int incd,
                          EREAL *a, int lda );
void BLAS_dgb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const double *d, int incd,
                          double *a, int lda );
void BLAS_egb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const EREAL *d, int incd,
                          EREAL *a, int lda );
void BLAS_cgb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const void *d, int incd,
                          void *a, int lda );
void BLAS_ygb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const void *d, int incd,
                          void *a, int lda );
void BLAS_zgb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const void *d, int incd,
                          void *a, int lda );
void BLAS_ygb_diag_scale( enum blas_order_type order,
                          enum blas_side_type side, int m, int n,
                          int kl, int ku, const void *d, int incd,
                          void *a, int lda );

void BLAS_sge_lrscale( enum blas_order_type order, int m, int n,
                       const float *dl, int incdl, const float *dr,
                       int incdr, float *a, int lda );
void BLAS_ege_lrscale( enum blas_order_type order, int m, int n,
                       const EREAL *dl, int incdl, const EREAL *dr,
                       int incdr, EREAL *a, int lda );
void BLAS_dge_lrscale( enum blas_order_type order, int m, int n,
                       const double *dl, int incdl, const double *dr,
                       int incdr, double *a, int lda );
void BLAS_ege_lrscale( enum blas_order_type order, int m, int n,
                       const EREAL *dl, int incdl, const EREAL *dr,
                       int incdr, EREAL *a, int lda );
void BLAS_cge_lrscale( enum blas_order_type order, int m, int n,
                       const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_yge_lrscale( enum blas_order_type order, int m, int n,
                       const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_zge_lrscale( enum blas_order_type order, int m, int n,
                       const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_yge_lrscale( enum blas_order_type order, int m, int n,
                       const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );

void BLAS_sgb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const float *dl, int incdl, const float *dr,
                       int incdr, float *a, int lda );
void BLAS_egb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const EREAL *dl, int incdl, const EREAL *dr,
                       int incdr, EREAL *a, int lda );
void BLAS_dgb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const double *dl, int incdl,
                       const double *dr, int incdr, double *a, int lda );
void BLAS_egb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const EREAL *dl, int incdl,
                       const EREAL *dr, int incdr, EREAL *a, int lda );
void BLAS_cgb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_ygb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_zgb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );
void BLAS_ygb_lrscale( enum blas_order_type order, int m, int n, int kl,
                       int ku, const void *dl, int incdl, const void *dr,
                       int incdr, void *a, int lda );

void BLAS_ssy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const float *d, int incd, float *a, int lda );
void BLAS_esy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const EREAL *d, int incd, EREAL *a, int lda );
void BLAS_dsy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const double *d, int incd, double *a,
                       int lda );
void BLAS_esy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const EREAL *d, int incd, EREAL *a,
                       int lda );
void BLAS_csy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_ysy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_zsy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_ysy_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );

void BLAS_ssb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const float *d, int incd, float *a,
                       int lda );
void BLAS_esb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const EREAL *d, int incd, EREAL *a,
                       int lda );
void BLAS_dsb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const double *d, int incd, double *a,
                       int lda );
void BLAS_esb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const EREAL *d, int incd, EREAL *a,
                       int lda );
void BLAS_csb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_ysb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_zsb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_ysb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );

void BLAS_ssp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const float *d, int incd, float *ap );
void BLAS_esp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const EREAL *d, int incd, EREAL *ap );
void BLAS_dsp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const double *d, int incd, double *ap );
void BLAS_esp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const EREAL *d, int incd, EREAL *ap );
void BLAS_csp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_ysp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_zsp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_ysp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );

void BLAS_che_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_yhe_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_zhe_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );
void BLAS_yhe_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *a, int lda );

void BLAS_chb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_yhb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_zhb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );
void BLAS_yhb_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, int k, const void *d, int incd, void *a,
                       int lda );

void BLAS_chp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_yhp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_zhp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );
void BLAS_yhp_lrscale( enum blas_order_type order, enum blas_uplo_type uplo,
                       int n, const void *d, int incd, void *ap );

void BLAS_sge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const float *b, int ldb, const float *d,
                              int incd, float *a, int lda );
void BLAS_ege_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const EREAL *b, int ldb, const EREAL *d,
                              int incd, EREAL *a, int lda );
void BLAS_dge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const double *b, int ldb, const double *d,
                              int incd, double *a, int lda );
void BLAS_ege_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const EREAL *b, int ldb, const EREAL *d,
                              int incd, EREAL *a, int lda );
void BLAS_cge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const void *b, int ldb, const void *d,
                              int incd, void *a, int lda );
void BLAS_yge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const void *b, int ldb, const void *d,
                              int incd, void *a, int lda );
void BLAS_zge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const void *b, int ldb, const void *d,
                              int incd, void *a, int lda );
void BLAS_yge_diag_scale_acc( enum blas_order_type order, int m, int n,
                              const void *b, int ldb, const void *d,
                              int incd, void *a, int lda );

void BLAS_sgb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const float *b, int ldb,
                              const float *d, int incd, float *a, int lda );
void BLAS_egb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const EREAL *b, int ldb,
                              const EREAL *d, int incd, EREAL *a, int lda );
void BLAS_dgb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const double *b, int ldb,
                              const double *d, int incd, double *a,
                              int lda );
void BLAS_egb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const EREAL *b, int ldb,
                              const EREAL *d, int incd, EREAL *a,
                              int lda );
void BLAS_cgb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const void *b, int ldb,
                              const void *d, int incd, void *a, int lda );
void BLAS_ygb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const void *b, int ldb,
                              const void *d, int incd, void *a, int lda );
void BLAS_zgb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const void *b, int ldb,
                              const void *d, int incd, void *a, int lda );
void BLAS_ygb_diag_scale_acc( enum blas_order_type order, int m, int n,
                              int kl, int ku, const void *b, int ldb,
                              const void *d, int incd, void *a, int lda );

void BLAS_sge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, float alpha, const float *a, int lda,
                   float beta, float *b, int ldb );
void BLAS_ege_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, EREAL alpha, const EREAL *a, int lda,
                   EREAL beta, EREAL *b, int ldb );
void BLAS_dge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, double alpha, const double *a, int lda,
                   double beta, double *b, int ldb );
void BLAS_ege_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, EREAL alpha, const EREAL *a, int lda,
                   EREAL beta, EREAL *b, int ldb );
void BLAS_cge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_yge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_zge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_yge_acc( enum blas_order_type order, enum blas_trans_type trans,
                   int m, int n, const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );

void BLAS_ssy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, float alpha,
                   const float *a, int lda, float beta, float *b, int ldb );
void BLAS_esy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b, int ldb );
void BLAS_dsy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, double alpha,
                   const double *a, int lda, double beta, double *b,
                   int ldb );
void BLAS_esy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b,
                   int ldb );
void BLAS_csy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ysy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_zsy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ysy_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );

void BLAS_ssb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k, float alpha,
                   const float *a, int lda, float beta, float *b, int ldb );
void BLAS_esb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b, int ldb );
void BLAS_dsb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k, double alpha,
                   const double *a, int lda, double beta, double *b,
                   int ldb );
void BLAS_esb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b,
                   int ldb );
void BLAS_csb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_ysb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_zsb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_ysb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, int k,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );

void BLAS_ssp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, float alpha,
                   const float *ap, float beta, float *bp );
void BLAS_esp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, EREAL *bp );
void BLAS_dsp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, double alpha,
                   const double *ap, double beta, double *bp );
void BLAS_esp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, EREAL *bp );
void BLAS_csp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_ysp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_zsp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_ysp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_trans_type trans, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );

void BLAS_sgb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   float alpha, const float *a, int lda, float beta,
                   float *b, int ldb );
void BLAS_egb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   EREAL alpha, const EREAL *a, int lda, EREAL beta,
                   EREAL *b, int ldb );
void BLAS_dgb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   double alpha, const double *a, int lda, double beta,
                   double *b, int ldb );
void BLAS_egb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   EREAL alpha, const EREAL *a, int lda, EREAL beta,
                   EREAL *b, int ldb );
void BLAS_cgb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_ygb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_zgb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );
void BLAS_ygb_acc( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, void *b, int ldb );

void BLAS_str_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, float alpha,
                   const float *a, int lda, float beta, float *b, int ldb );
void BLAS_etr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b, int ldb );
void BLAS_dtr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, double alpha,
                   const double *a, int lda, double beta, double *b,
                   int ldb );
void BLAS_etr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b,
                   int ldb );
void BLAS_ctr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b, 
                   int ldb );
void BLAS_ytr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b, 
                   int ldb );
void BLAS_ztr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ytr_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );

void BLAS_stb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, float alpha,
                   const float *a, int lda, float beta, float *b, int ldb );
void BLAS_etb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b, int ldb );
void BLAS_dtb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, double alpha,
                   const double *a, int lda, double beta, double *b,
                   int ldb );
void BLAS_etb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, EREAL *b,
                   int ldb );
void BLAS_ctb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ytb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ztb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );
void BLAS_ytb_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, void *b,
                   int ldb );

void BLAS_stp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, float alpha,
                   const float *ap, float beta, float *bp );
void BLAS_etp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, EREAL *bp );
void BLAS_dtp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, double alpha,
                   const double *ap, double beta, double *bp );
void BLAS_etp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, EREAL *bp );
void BLAS_ctp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_ytp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_ztp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );
void BLAS_ytp_acc( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, void *bp );

void BLAS_sge_add( enum blas_order_type order, int m, int n, float alpha,
                   const float *a, int lda, float beta, const float *b,
                   int ldb, float *c, int ldc );
void BLAS_ege_add( enum blas_order_type order, int m, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_dge_add( enum blas_order_type order, int m, int n, double alpha,
                   const double *a, int lda, double beta, const double *b,
                   int ldb, double *c, int ldc );
void BLAS_ege_add( enum blas_order_type order, int m, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_cge_add( enum blas_order_type order, int m, int n,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_yge_add( enum blas_order_type order, int m, int n,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_zge_add( enum blas_order_type order, int m, int n,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_yge_add( enum blas_order_type order, int m, int n,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b,
                   int ldb, void *c, int ldc );

void BLAS_sgb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   float alpha, const float *a, int lda, float beta,
                   const float *b, int ldb, float *c, int ldc );
void BLAS_egb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   EREAL alpha, const EREAL *a, int lda, EREAL beta,
                   const EREAL *b, int ldb, EREAL *c, int ldc );
void BLAS_dgb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   double alpha, const double *a, int lda, double beta,
                   const double *b, int ldb, double *c, int ldc );
void BLAS_egb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   EREAL alpha, const EREAL *a, int lda, EREAL beta,
                   const EREAL *b, int ldb, EREAL *c, int ldc );
void BLAS_cgb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ygb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_zgb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ygb_add( enum blas_order_type order, int m, int n, int kl, int ku,
                   const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );

void BLAS_ssy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, float alpha, const float *a, int lda, float beta,
                   const float *b, int ldb, float *c, int ldc );
void BLAS_esy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, EREAL alpha, const EREAL *a, int lda, EREAL beta,
                   const EREAL *b, int ldb, EREAL *c, int ldc );
void BLAS_dsy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, double alpha, const double *a, int lda,
                   double beta, const double *b, int ldb, double *c,
                   int ldc );
void BLAS_esy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, EREAL alpha, const EREAL *a, int lda,
                   EREAL beta, const EREAL *b, int ldb, EREAL *c,
                   int ldc );
void BLAS_csy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ysy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_zsy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ysy_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );

void BLAS_ssb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, float alpha, const float *a, int lda,
                   float beta, const float *b, int ldb, float *c, int ldc );
void BLAS_esb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, EREAL alpha, const EREAL *a, int lda,
                   EREAL beta, const EREAL *b, int ldb, EREAL *c, int ldc );
void BLAS_dsb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, double alpha, const double *a, int lda,
                   double beta, const double *b, int ldb, double *c, 
                   int ldc );
void BLAS_esb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, EREAL alpha, const EREAL *a, int lda,
                   EREAL beta, const EREAL *b, int ldb, EREAL *c, 
                   int ldc );
void BLAS_csb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ysb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_zsb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );
void BLAS_ysb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, int k, const void *alpha, const void *a, int lda,
                   const void *beta, const void *b, int ldb, void *c,
                   int ldc );

void BLAS_ssp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, float alpha, const float *ap, float beta,
                   const float *bp, float *cp );
void BLAS_esp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, EREAL alpha, const EREAL *ap, EREAL beta,
                   const EREAL *bp, EREAL *cp );
void BLAS_dsp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, double alpha, const double *ap, double beta,
                   const double *bp, double *cp );
void BLAS_esp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, EREAL alpha, const EREAL *ap, EREAL beta,
                   const EREAL *bp, EREAL *cp );
void BLAS_csp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *ap,
                   const void *beta, const void *bp, void *cp );
void BLAS_ysp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *ap,
                   const void *beta, const void *bp, void *cp );
void BLAS_zsp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *ap,
                   const void *beta, const void *bp, void *cp );
void BLAS_ysp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   int n, const void *alpha, const void *ap,
                   const void *beta, const void *bp, void *cp );

void BLAS_str_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, float alpha,
                   const float *a, int lda, float beta, const float *b,
                   int ldb, float *c, int ldc );
void BLAS_etr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_dtr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, double alpha,
                   const double *a, int lda, double beta, const double *b,
                   int ldb, double *c, int ldc );
void BLAS_etr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_ctr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ytr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ztr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ytr_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );

void BLAS_stb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, float alpha,
                   const float *a, int lda, float beta, const float *b,
                   int ldb, float *c, int ldc );
void BLAS_etb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_dtb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, double alpha,
                   const double *a, int lda, double beta, const double *b,
                   int ldb, double *c, int ldc );
void BLAS_etb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, EREAL alpha,
                   const EREAL *a, int lda, EREAL beta, const EREAL *b,
                   int ldb, EREAL *c, int ldc );
void BLAS_ctb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ytb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ztb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );
void BLAS_ytb_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, int k, const void *alpha,
                   const void *a, int lda, const void *beta, const void *b,
                   int ldb, void *c, int ldc );

void BLAS_stp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, float alpha,
                   const float *ap, float beta, const float *bp,
                   float *cp );
void BLAS_etp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, const EREAL *bp,
                   EREAL *cp );
void BLAS_dtp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, double alpha,
                   const double *ap, double beta, const double *bp,
                   double *cp );
void BLAS_etp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, EREAL alpha,
                   const EREAL *ap, EREAL beta, const EREAL *bp,
                   EREAL *cp );
void BLAS_ctp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, const void *bp,
                   void *cp );
void BLAS_ytp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, const void *bp,
                   void *cp );
void BLAS_ztp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, const void *bp,
                   void *cp );
void BLAS_ytp_add( enum blas_order_type order, enum blas_uplo_type uplo,
                   enum blas_diag_type diag, int n, const void *alpha,
                   const void *ap, const void *beta, const void *bp,
                   void *cp );

  /* Matrix-Matrix Operations */

void BLAS_sgemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 float alpha, const float *a, int lda, const float *b,
                 int ldb, float beta, float *c, int ldc );
void BLAS_egemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 EREAL alpha, const EREAL *a, int lda, const EREAL *b,
                 int ldb, EREAL beta, EREAL *c, int ldc );
void BLAS_dgemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 double alpha, const double *a, int lda, const double *b,
                 int ldb, double beta, double *c, int ldc );
void BLAS_egemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 EREAL alpha, const EREAL *a, int lda, const EREAL *b,
                 int ldb, EREAL beta, EREAL *c, int ldc );
void BLAS_cgemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 const void *alpha, const void *a, int lda, const void *b,
                 int ldb, const void *beta, void *c, int ldc );
void BLAS_ygemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 const void *alpha, const void *a, int lda, const void *b,
                 int ldb, const void *beta, void *c, int ldc );
void BLAS_zgemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 const void *alpha, const void *a, int lda, const void *b,
                 int ldb, const void *beta, void *c, int ldc );
void BLAS_ygemm( enum blas_order_type order, enum blas_trans_type transa,
                 enum blas_trans_type transb, int m, int n, int k,
                 const void *alpha, const void *a, int lda, const void *b,
                 int ldb, const void *beta, void *c, int ldc );

void BLAS_ssymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, float alpha,
                 const float *a, int lda, const float *b, int ldb,
                 float beta, float *c, int ldc );
void BLAS_esymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, EREAL alpha,
                 const EREAL *a, int lda, const EREAL *b, int ldb,
                 EREAL beta, EREAL *c, int ldc );
void BLAS_dsymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, double alpha,
                 const double *a, int lda, const double *b, int ldb,
                 double beta, double *c, int ldc );
void BLAS_esymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, EREAL alpha,
                 const EREAL *a, int lda, const EREAL *b, int ldb,
                 EREAL beta, EREAL *c, int ldc );
void BLAS_csymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_ysymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_zsymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_ysymm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );

void BLAS_chemm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_yhemm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_zhemm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );
void BLAS_yhemm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, int m, int n, const void *alpha,
                 const void *a, int lda, const void *b, int ldb,
                 const void *beta, void *c, int ldc );

void BLAS_strmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, float alpha,
                 const float *t, int ldt, float *b, int ldb );
void BLAS_etrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, EREAL alpha,
                 const EREAL *t, int ldt, EREAL *b, int ldb );
void BLAS_dtrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, double alpha,
                 const double *t, int ldt, double *b, int ldb );
void BLAS_etrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, EREAL alpha,
                 const EREAL *t, int ldt, EREAL *b, int ldb );
void BLAS_ctrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ytrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ztrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ytrmm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );

void BLAS_strsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, float alpha,
                 const float *t, int ldt, float *b, int ldb );
void BLAS_etrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, EREAL alpha,
                 const EREAL *t, int ldt, EREAL *b, int ldb );
void BLAS_dtrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, double alpha,
                 const double *t, int ldt, double *b, int ldb );
void BLAS_etrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, EREAL alpha,
                 const EREAL *t, int ldt, EREAL *b, int ldb );
void BLAS_ctrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ytrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ztrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );
void BLAS_ytrsm( enum blas_order_type order, enum blas_side_type side,
                 enum blas_uplo_type uplo, enum blas_trans_type transt,
                 enum blas_diag_type diag, int m, int n, const void *alpha,
                 const void *t, int ldt, void *b, int ldb );

void BLAS_ssyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, float alpha,
                 const float *a, int lda, float beta, float *c, int ldc );
void BLAS_esyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, EREAL alpha,
                 const EREAL *a, int lda, EREAL beta, EREAL *c, int ldc );
void BLAS_dsyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, double alpha,
                 const double *a, int lda, double beta, double *c, int ldc );
void BLAS_esyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, EREAL alpha,
                 const EREAL *a, int lda, EREAL beta, EREAL *c, int ldc );
void BLAS_csyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, const void *alpha,
                 const void *a, int lda, const void *beta, void *c,
                 int ldc );
void BLAS_ysyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, const void *alpha,
                 const void *a, int lda, const void *beta, void *c,
                 int ldc );
void BLAS_zsyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, const void *alpha,
                 const void *a, int lda, const void *beta, void *c,
                 int ldc );
void BLAS_ysyrk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, const void *alpha,
                 const void *a, int lda, const void *beta, void *c,
                 int ldc );

void BLAS_cherk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, float alpha,
                 const void *a, int lda, float beta, void *c, int ldc );
void BLAS_yherk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, EREAL alpha,
                 const void *a, int lda, EREAL beta, void *c, int ldc );
void BLAS_zherk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, double alpha,
                 const void *a, int lda, double beta, void *c, int ldc );
void BLAS_yherk( enum blas_order_type order, enum blas_uplo_type uplo,
                 enum blas_trans_type trans, int n, int k, EREAL alpha,
                 const void *a, int lda, EREAL beta, void *c, int ldc );

void BLAS_ssy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          float alpha, const float *a, int lda, 
                          const float *d, const float *e, float beta,
                          float *c, int ldc );
void BLAS_esy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          EREAL alpha, const EREAL *a, int lda, 
                          const EREAL *d, const EREAL *e, EREAL beta,
                          EREAL *c, int ldc );
void BLAS_dsy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          double alpha, const double *a, int lda,
                          const double *d, const double *e, double beta,
                          double *c, int ldc );
void BLAS_esy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          EREAL alpha, const EREAL *a, int lda,
                          const EREAL *d, const EREAL *e, EREAL beta,
                          EREAL *c, int ldc );
void BLAS_csy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          const void *alpha, const void *a, int lda,
                          const void *d, const void *e, const void *beta,
                          void *c, int ldc );
void BLAS_ysy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          const void *alpha, const void *a, int lda,
                          const void *d, const void *e, const void *beta,
                          void *c, int ldc );
void BLAS_zsy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          const void *alpha, const void *a, int lda,
                          const void *d, const void *e, const void *beta,
                          void *c, int ldc );
void BLAS_ysy_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          const void *alpha, const void *a, int lda,
                          const void *d, const void *e, const void *beta,
                          void *c, int ldc );

void BLAS_che_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          float alpha, const void *a, int lda,
                          const void *d, const void *e, float beta,
                          void *c, int ldc );
void BLAS_yhe_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          EREAL alpha, const void *a, int lda,
                          const void *d, const void *e, EREAL beta,
                          void *c, int ldc );
void BLAS_zhe_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          double alpha, const void *a, int lda,
                          const void *d, const void *e, double beta,
                          void *c, int ldc );
void BLAS_yhe_tridiag_rk( enum blas_order_type order,
                          enum blas_uplo_type uplo,
                          enum blas_trans_type trans, int n, int k,
                          EREAL alpha, const void *a, int lda,
                          const void *d, const void *e, EREAL beta,
                          void *c, int ldc );

void BLAS_ssyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, float alpha,
                  const float *a, int lda, const float *b, int ldb,
                  float beta, float *c, int ldc );
void BLAS_esyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, EREAL alpha,
                  const EREAL *a, int lda, const EREAL *b, int ldb,
                  EREAL beta, EREAL *c, int ldc );
void BLAS_dsyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, double alpha,
                  const double *a, int lda, const double *b, int ldb,
                  double beta, double *c, int ldc );
void BLAS_esyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, EREAL alpha,
                  const EREAL *a, int lda, const EREAL *b, int ldb,
                  EREAL beta, EREAL *c, int ldc );
void BLAS_csyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *a, int lda, const void *b,
                  int ldb, const void *beta, void *c, int ldc );
void BLAS_ysyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *a, int lda, const void *b,
                  int ldb, const void *beta, void *c, int ldc );
void BLAS_zsyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, 
                  const void *alpha, const void *a, int lda, const void *b,
                  int ldb, const void *beta, void *c, int ldc );
void BLAS_ysyr2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k, 
                  const void *alpha, const void *a, int lda, const void *b,
                  int ldb, const void *beta, void *c, int ldc );

void BLAS_cher2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *A, int lda, const void *b,
                  int ldb, float beta, void *c, int ldc );
void BLAS_yher2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *A, int lda, const void *b,
                  int ldb, EREAL beta, void *c, int ldc );
void BLAS_zher2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *A, int lda, const void *b,
                  int ldb, double beta, void *c, int ldc );
void BLAS_yher2k( enum blas_order_type order, enum blas_uplo_type uplo,
                  enum blas_trans_type trans, int n, int k,
                  const void *alpha, const void *A, int lda, const void *b,
                  int ldb, EREAL beta, void *c, int ldc );

void BLAS_ssy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           float alpha, const float *a, int lda,
                           const float *d, const float *e, const float *b,
                           int ldb, float beta, float *c, int ldc );
void BLAS_esy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           EREAL alpha, const EREAL *a, int lda,
                           const EREAL *d, const EREAL *e, const EREAL *b,
                           int ldb, EREAL beta, EREAL *c, int ldc );
void BLAS_dsy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           double alpha, const double *a, int lda,
                           const double *d, const double *e, const double *b,
                           int ldb, double beta, double *c, int ldc );
void BLAS_esy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           EREAL alpha, const EREAL *a, int lda,
                           const EREAL *d, const EREAL *e, const EREAL *b,
                           int ldb, EREAL beta, EREAL *c, int ldc );
void BLAS_csy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, const void *beta, void *c, int ldc );
void BLAS_ysy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, const void *beta, void *c, int ldc );
void BLAS_zsy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, const void *beta, void *c, int ldc );
void BLAS_ysy_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, const void *beta, void *c, int ldc );

void BLAS_che_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, float beta, void *c, int ldc );
void BLAS_yhe_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, EREAL beta, void *c, int ldc );
void BLAS_zhe_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, double beta, void *c, int ldc );
void BLAS_yhe_tridiag_r2k( enum blas_order_type order,
                           enum blas_uplo_type uplo,
                           enum blas_trans_type trans, int n, int k,
                           const void *alpha, const void *a, int lda,
                           const void *d, const void *e, const void *b,
                           int ldb, EREAL beta, void *c, int ldc );

  /* Data Movement with Matrices */

void BLAS_sge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const float *a, int lda, float *b,
                    int ldb );
void BLAS_ege_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const EREAL *a, int lda, EREAL *b,
                    int ldb );
void BLAS_dge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const double *a, int lda, double *b,
                    int ldb );
void BLAS_ege_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const EREAL *a, int lda, EREAL *b,
                    int ldb );
void BLAS_cge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const void *a, int lda, void *b, int ldb );
void BLAS_yge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const void *a, int lda, void *b, int ldb );
void BLAS_zge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const void *a, int lda, void *b, int ldb );
void BLAS_yge_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, const void *a, int lda, void *b, int ldb );

void BLAS_sgb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const float *a, int lda,
                    float *b, int ldb );
void BLAS_egb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const EREAL *a, int lda,
                    EREAL *b, int ldb );
void BLAS_dgb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const double *a, int lda,
                    double *b, int ldb );
void BLAS_egb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const EREAL *a, int lda,
                    EREAL *b, int ldb );
void BLAS_cgb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    void *b, int ldb );
void BLAS_ygb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    void *b, int ldb );
void BLAS_zgb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    void *b, int ldb );
void BLAS_ygb_copy( enum blas_order_type order, enum blas_trans_type trans,
                    int m, int n, int kl, int ku, const void *a, int lda,
                    void *b, int ldb );

void BLAS_ssy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const float *a, int lda, float *b, int ldb );
void BLAS_esy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const EREAL *a, int lda, EREAL *b, int ldb );
void BLAS_dsy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const double *a, int lda, double *b, int ldb );
void BLAS_esy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const EREAL *a, int lda, EREAL *b, int ldb );
void BLAS_csy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_ysy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_zsy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_ysy_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );

void BLAS_ssb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const float *a, int lda, float *b,
                    int ldb );
void BLAS_esb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const EREAL *a, int lda, EREAL *b,
                    int ldb );
void BLAS_dsb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const double *a, int lda, double *b,
                    int ldb );
void BLAS_esb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const EREAL *a, int lda, EREAL *b,
                    int ldb );
void BLAS_csb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_ysb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_zsb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_ysb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );

void BLAS_ssp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const float *ap, float *bp );
void BLAS_esp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const EREAL *ap, EREAL *bp );
void BLAS_dsp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const double *ap, double *bp );
void BLAS_esp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const EREAL *ap, EREAL *bp );
void BLAS_csp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_ysp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_zsp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_ysp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );

void BLAS_str_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const float *a, int lda, float *b, int ldb );
void BLAS_etr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const EREAL *a, int lda, EREAL *b, int ldb );
void BLAS_dtr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const double *a, int lda, double *b, int ldb );
void BLAS_etr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const EREAL *a, int lda, EREAL *b, int ldb );
void BLAS_ctr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_ytr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_ztr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_ytr_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *a, int lda, void *b, int ldb );

void BLAS_stb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const float *a, int lda, float *b, 
                    int ldb );
void BLAS_etb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const EREAL *a, int lda, EREAL *b, 
                    int ldb );
void BLAS_dtb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const double *a, int lda, double *b,
                    int ldb );
void BLAS_etb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const EREAL *a, int lda, EREAL *b,
                    int ldb );
void BLAS_ctb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_ytb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_ztb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_ytb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, int k, const void *a, int lda, void *b, int ldb );

void BLAS_stp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const float *ap, float *bp );
void BLAS_etp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const EREAL *ap, EREAL *bp );
void BLAS_dtp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const double *ap, double *bp );
void BLAS_etp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const EREAL *ap, EREAL *bp );
void BLAS_ctp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *ap, void *bp );
void BLAS_ytp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *ap, void *bp );
void BLAS_ztp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *ap, void *bp );
void BLAS_ytp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    enum blas_trans_type trans, enum blas_diag_type diag,
                    int n, const void *ap, void *bp );

void BLAS_che_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_yhe_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_zhe_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );
void BLAS_yhe_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *a, int lda, void *b, int ldb );

void BLAS_chb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_yhb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_zhb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );
void BLAS_yhb_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, int k, const void *a, int lda, void *b, int ldb );

void BLAS_chp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_yhp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_zhp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );
void BLAS_yhp_copy( enum blas_order_type order, enum blas_uplo_type uplo,
                    int n, const void *ap, void *bp );

void BLAS_sge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, float *a, int lda );
void BLAS_ege_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, EREAL *a, int lda );
void BLAS_dge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, double *a, int lda );
void BLAS_ege_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, EREAL *a, int lda );
void BLAS_cge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, void *a, int lda );
void BLAS_yge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, void *a, int lda );
void BLAS_zge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, void *a, int lda );
void BLAS_yge_trans( enum blas_order_type order, enum blas_conj_type conj,
                     int n, void *a, int lda );

void BLAS_sge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, float *a,
                       int lda );
void BLAS_ege_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, EREAL *a,
                       int lda );
void BLAS_dge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, double *a,
                       int lda );
void BLAS_ege_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, EREAL *a,
                       int lda );
void BLAS_cge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, void *a,
                       int lda );
void BLAS_yge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, void *a,
                       int lda );
void BLAS_zge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, void *a,
                       int lda );
void BLAS_yge_permute( enum blas_order_type order, enum blas_side_type side,
                       int m, int n, const int *p, int incp, void *a,
                       int lda );

  /* Environmental Enquiry */

float c_sfpinfo( enum blas_cmach_type cmach );
EREAL c_efpinfo( enum blas_cmach_type cmach );
double c_dfpinfo( enum blas_cmach_type cmach );
EREAL c_efpinfo( enum blas_cmach_type cmach );

#endif
   /* BLAS_DENSE_PROTO_H */
