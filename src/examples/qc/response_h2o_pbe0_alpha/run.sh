#!/bin/sh
exec ${MADQC:-madqc} --wf=response response_h2o_pbe0_alpha.in
