#! /bin/bash
# Run this script to run the tests and check if the results are equal to the ones stored in reference/
# log is stored in report.txt
{
rm -f *.output > /dev/null 2>/dev/null
echo "Running the test scripts..."
time ./test.sh >log 2>err
# here you can put more scripts, of course...
#rm -f out2-*.dat > /dev/null 2>/dev/null
#echo "Running the density scripts..."
#time ./test2.sh >>log 2>>err

echo "Done."

if ls reference/* > /dev/null
then
for file in reference/* ;
do
  new="${file:10}"
  echo $new
  if test -f "$new" ; then
    out="$(diff "$file" "$new")"
  test -n "$out" && {
      echo FAILURE
      echo "Diff for ${file}:"
      diff "${file}" "$new"
    }
  else
    echo FAILURE
    echo FILE $new does not exist
  fi
done
else
    echo WARNING
    echo no file has been checked
fi

} | tee report.txt
