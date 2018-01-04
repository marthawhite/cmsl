#!/bin/sh
#Check if file begins with blank line, if not, insert a blank line

echo 'Outputting results to test_output'

echo 'Results for test on' $(date +"%m-%d-%Y") > test_output

for f in *.m;
do
  firstline="`head -1 "${f}"`"

  if [ "${firstline}" = "%TEST" ]; then
    echo "+++++++++++++++++++++++++ Running test for ${f}"
    echo "+++++++++++++++++++++++++ Running test for ${f}" >> test_output
    matlab -nodisplay < ${f} >> test_output
    echo "+++++++++++++++++++++++++ Done test for ${f}" >> test_output   
  fi

done
