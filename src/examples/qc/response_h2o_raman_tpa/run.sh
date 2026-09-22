#!/bin/sh
exec ${MADQC:-madqc} --wf=response response_h2o_raman_tpa.in
