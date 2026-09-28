#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_li_pcm_born.in
