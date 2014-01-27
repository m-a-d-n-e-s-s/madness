#! /bin/sh

aclocal -I ./config
autoconf
autoheader
automake --add-missing
