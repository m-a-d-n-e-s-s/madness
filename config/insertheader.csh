#!/bin/csh

set files = (`find . \( -name "*.h" -o -name "*.cc" -o -name "*.c" -o -name "*.hpp" -o -name "*.cpp" \)`)
set files = (`grep -L 'This program is free software' $files | grep -v mainpage | grep -v lookup3 `)

foreach file ($files)
  echo Processing $file
  cat config/HEADER $file > tmptmptmp
  mv tmptmptmp $file
end


