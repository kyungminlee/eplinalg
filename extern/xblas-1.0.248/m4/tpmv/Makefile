include ../../make.conf
include ../../$(MAKEINC)


SRC_PATH = ../../src/tpmv
M4_PATH = -I. -I..

TPMV_SRCS = \
	$(SRC_PATH)/BLAS_dtpmv_s.c\
	$(SRC_PATH)/BLAS_ztpmv_c.c\
	$(SRC_PATH)/BLAS_ctpmv_s.c\
	$(SRC_PATH)/BLAS_ztpmv_d.c\
	$(SRC_PATH)/BLAS_stpmv_x.c\
	$(SRC_PATH)/BLAS_dtpmv_x.c\
	$(SRC_PATH)/BLAS_ctpmv_x.c\
	$(SRC_PATH)/BLAS_ztpmv_x.c\
	$(SRC_PATH)/BLAS_dtpmv_s_x.c\
	$(SRC_PATH)/BLAS_ztpmv_c_x.c\
	$(SRC_PATH)/BLAS_ctpmv_s_x.c\
	$(SRC_PATH)/BLAS_ztpmv_d_x.c \
	\
	$(SRC_PATH)/BLAS_dtpmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztpmv_c-f2c.c\
	$(SRC_PATH)/BLAS_ctpmv_s-f2c.c\
	$(SRC_PATH)/BLAS_ztpmv_d-f2c.c\
	$(SRC_PATH)/BLAS_stpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ctpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_ztpmv_x-f2c.c\
	$(SRC_PATH)/BLAS_dtpmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztpmv_c_x-f2c.c\
	$(SRC_PATH)/BLAS_ctpmv_s_x-f2c.c\
	$(SRC_PATH)/BLAS_ztpmv_d_x-f2c.c

tpmv-source: $(TPMV_SRCS)

tpmv-header: tpmv.h

$(SRC_PATH)/BLAS_dtpmv_s.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(d, s)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtpmv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_s.c && \
    mv tmp.BLAS_dtpmv_s.c $@ && rm -f tmp.BLAS_dtpmv_s.c*

$(SRC_PATH)/BLAS_ztpmv_c.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(z, c)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztpmv_c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_c.c && \
    mv tmp.BLAS_ztpmv_c.c $@ && rm -f tmp.BLAS_ztpmv_c.c*

$(SRC_PATH)/BLAS_ctpmv_s.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(c, s)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctpmv_s.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_s.c && \
    mv tmp.BLAS_ctpmv_s.c $@ && rm -f tmp.BLAS_ctpmv_s.c*

$(SRC_PATH)/BLAS_ztpmv_d.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(z, d)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztpmv_d.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_d.c && \
    mv tmp.BLAS_ztpmv_d.c $@ && rm -f tmp.BLAS_ztpmv_d.c*

$(SRC_PATH)/BLAS_stpmv_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(s, s, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_stpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stpmv_x.c && \
    mv tmp.BLAS_stpmv_x.c $@ && rm -f tmp.BLAS_stpmv_x.c*

$(SRC_PATH)/BLAS_dtpmv_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(d, d, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_x.c && \
    mv tmp.BLAS_dtpmv_x.c $@ && rm -f tmp.BLAS_dtpmv_x.c*

$(SRC_PATH)/BLAS_ctpmv_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(c, c, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_x.c && \
    mv tmp.BLAS_ctpmv_x.c $@ && rm -f tmp.BLAS_ctpmv_x.c*

$(SRC_PATH)/BLAS_ztpmv_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(z, z, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztpmv_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_x.c && \
    mv tmp.BLAS_ztpmv_x.c $@ && rm -f tmp.BLAS_ztpmv_x.c*

$(SRC_PATH)/BLAS_dtpmv_s_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(d, s, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_dtpmv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_s_x.c && \
    mv tmp.BLAS_dtpmv_s_x.c $@ && rm -f tmp.BLAS_dtpmv_s_x.c*

$(SRC_PATH)/BLAS_ztpmv_c_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(z, c, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztpmv_c_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_c_x.c && \
    mv tmp.BLAS_ztpmv_c_x.c $@ && rm -f tmp.BLAS_ztpmv_c_x.c*

$(SRC_PATH)/BLAS_ctpmv_s_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(c, s, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ctpmv_s_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_s_x.c && \
    mv tmp.BLAS_ctpmv_s_x.c $@ && rm -f tmp.BLAS_ctpmv_s_x.c*

$(SRC_PATH)/BLAS_ztpmv_d_x.c: tpmv.m4 tpmv-common.m4 ../cblas.m4
	echo "TPMV(z, d, _x)" | cat tpmv.m4 - | \
    $(M4) $(M4_OPTS) $(M4_PATH) - > tmp.BLAS_ztpmv_d_x.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_d_x.c && \
    mv tmp.BLAS_ztpmv_d_x.c $@ && rm -f tmp.BLAS_ztpmv_d_x.c*



$(SRC_PATH)/BLAS_dtpmv_s-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="d, s" ../f2c-bridge.m4 > tmp.BLAS_dtpmv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_s-f2c.c && \
    mv tmp.BLAS_dtpmv_s-f2c.c $@ && rm -f tmp.BLAS_dtpmv_s-f2c.c*

$(SRC_PATH)/BLAS_ztpmv_c-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="z, c" ../f2c-bridge.m4 > tmp.BLAS_ztpmv_c-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_c-f2c.c && \
    mv tmp.BLAS_ztpmv_c-f2c.c $@ && rm -f tmp.BLAS_ztpmv_c-f2c.c*

$(SRC_PATH)/BLAS_ctpmv_s-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="c, s" ../f2c-bridge.m4 > tmp.BLAS_ctpmv_s-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_s-f2c.c && \
    mv tmp.BLAS_ctpmv_s-f2c.c $@ && rm -f tmp.BLAS_ctpmv_s-f2c.c*

$(SRC_PATH)/BLAS_ztpmv_d-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="z, d" ../f2c-bridge.m4 > tmp.BLAS_ztpmv_d-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_d-f2c.c && \
    mv tmp.BLAS_ztpmv_d-f2c.c $@ && rm -f tmp.BLAS_ztpmv_d-f2c.c*

$(SRC_PATH)/BLAS_stpmv_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="s, s, _x" ../f2c-bridge.m4 > tmp.BLAS_stpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_stpmv_x-f2c.c && \
    mv tmp.BLAS_stpmv_x-f2c.c $@ && rm -f tmp.BLAS_stpmv_x-f2c.c*

$(SRC_PATH)/BLAS_dtpmv_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="d, d, _x" ../f2c-bridge.m4 > tmp.BLAS_dtpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_x-f2c.c && \
    mv tmp.BLAS_dtpmv_x-f2c.c $@ && rm -f tmp.BLAS_dtpmv_x-f2c.c*

$(SRC_PATH)/BLAS_ctpmv_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="c, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ctpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_x-f2c.c && \
    mv tmp.BLAS_ctpmv_x-f2c.c $@ && rm -f tmp.BLAS_ctpmv_x-f2c.c*

$(SRC_PATH)/BLAS_ztpmv_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="z, z, _x" ../f2c-bridge.m4 > tmp.BLAS_ztpmv_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_x-f2c.c && \
    mv tmp.BLAS_ztpmv_x-f2c.c $@ && rm -f tmp.BLAS_ztpmv_x-f2c.c*

$(SRC_PATH)/BLAS_dtpmv_s_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="d, s, _x" ../f2c-bridge.m4 > tmp.BLAS_dtpmv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_dtpmv_s_x-f2c.c && \
    mv tmp.BLAS_dtpmv_s_x-f2c.c $@ && rm -f tmp.BLAS_dtpmv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztpmv_c_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="z, c, _x" ../f2c-bridge.m4 > tmp.BLAS_ztpmv_c_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_c_x-f2c.c && \
    mv tmp.BLAS_ztpmv_c_x-f2c.c $@ && rm -f tmp.BLAS_ztpmv_c_x-f2c.c*

$(SRC_PATH)/BLAS_ctpmv_s_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="c, s, _x" ../f2c-bridge.m4 > tmp.BLAS_ctpmv_s_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ctpmv_s_x-f2c.c && \
    mv tmp.BLAS_ctpmv_s_x-f2c.c $@ && rm -f tmp.BLAS_ctpmv_s_x-f2c.c*

$(SRC_PATH)/BLAS_ztpmv_d_x-f2c.c: ../f2c-bridge.m4 tpmv-common.m4 ../cblas.m4
	$(M4) $(M4_OPTS) $(M4_PATH) -Droutine=tpmv -DRARG="z, d, _x" ../f2c-bridge.m4 > tmp.BLAS_ztpmv_d_x-f2c.c && \
    $(INDENT) $(INDENT_OPTS) tmp.BLAS_ztpmv_d_x-f2c.c && \
    mv tmp.BLAS_ztpmv_d_x-f2c.c $@ && rm -f tmp.BLAS_ztpmv_d_x-f2c.c*


tpmv.h: ../cblas.m4 ../header.template tpmv-common.m4
	sed -e "s/routine/tpmv/g" -e "s/ROUTINE/TPMV/g" \
  ../header.template | $(M4) $(M4_OPTS) $(M4_PATH) - >$@

source-clean:
	rm -f $(TPMV_SRCS) tpmv.h

.PHONY: tpmv-source tpmv-header source-clean

