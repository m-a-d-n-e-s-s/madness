#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_h2o_hf.in
