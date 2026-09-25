#!/bin/sh
exec ${MADQC:-madqc} --wf=response response_h2o_lda_alpha.in
