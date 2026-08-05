#!/bin/sh
exec ${MADQC:-madqc} --optimize --wf=nemo nemo_lih_optimize.in
