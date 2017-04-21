#!/bin/bash
# this script run the test and save the results as new references. To use when something new comes out. Be careful though
echo 'Running the scripts...'
./test.sh >log 2>err
# here you can put other test scripts
#./test2.sh >log 2>err
echo 'Done.'

for file in *output
do
  mv "${file}" "reference/${file}"
done

