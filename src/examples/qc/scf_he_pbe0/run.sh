#!/bin/sh
exec ${MADQC:-madqc} --wf=scf scf_he_pbe0.in
