#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_lih_pbe_d3.in
