AC_DEFUN([ACX_MAC], [
        # If on a MAC
        #   - kiss steve jobs' ass
        #   - get an iphone
        #   - get another mac
        ON_A_MAC="no"
        uname -a | grep -iq Darwin
        if test $? = 0; then
            ITS_A_MAC="yes"
            ON_A_MAC="yes"
            LDFLAGS="-Wl,-no_pie $LDFLAGS"
            AC_MSG_NOTICE([You are building on a mac ... now tell ten of your friends.])
            AC_DEFINE(ON_A_MAC,[1],[Set if building on a mac])
        fi
])




