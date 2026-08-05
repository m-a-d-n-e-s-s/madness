#!/bin/sh
exec ${MADQC:-madqc} --wf=optimize optimize_lih_moldft.in
