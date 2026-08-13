#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_lih_gopt.in
