#!/bin/sh
exec ${MADQC:-madqc} --wf=response response_h2o_tpa_seeded.in
