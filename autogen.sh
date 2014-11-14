#! /bin/sh

set -e

aclocal -I ./config
autoconf
autoheader
automake --add-missing
