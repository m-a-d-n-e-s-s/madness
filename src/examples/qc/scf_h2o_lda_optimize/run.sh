#!/bin/sh
exec ${MADQC:-madqc} --optimize --wf=scf scf_h2o_lda_optimize.in
