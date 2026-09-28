#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_lih_pcm_water.in
