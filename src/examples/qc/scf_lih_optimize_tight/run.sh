#!/bin/sh
exec ${MADQC:-madqc} --optimize --wf=scf scf_lih_optimize_tight.in
