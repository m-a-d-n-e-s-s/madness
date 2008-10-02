AC_DEFUN([ACX_MAC], [
        # If on a MAC
        #   - tear hair out
        #   - stomp on ipod
        #   - toss iphone in toilet
        #   - call george
        uname -a | grep -iq Darwin
        ON_A_MAC="no"
        if test $? = 0; then
            ON_A_MAC="yes"
            echo "Sorry ... you are building on a mac"
            AC_DEFINE(ITS_A_MAC,[1],[Set if building on a mac])
        fi
])




