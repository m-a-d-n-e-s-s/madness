#! /bin/sh

set -e

export ACLOCAL_PATH=/usr/share/aclocal

aclocal -I ./config
autoconf
autoheader
if hash glibtoolize 2>/dev/null; then
  glibtoolize
else
  libtoolize
fi
automake --add-missing
