#! /bin/sh

set -e

aclocal -I ./config
autoconf
autoheader
if hash glibtoolize 2>/dev/null; then
  glibtoolize
else
  libtoolize
fi
automake --add-missing
