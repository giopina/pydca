{
echo 'running octave script'
cat <<EOF >tmp_script.m
dca('stripped_286.afa_filtered200','test.out')
EOF
octave tmp_script.m >log 2>err

echo 'DCA done!'

new=test.out
file=reference/stripped_286.afa_filtered200.output

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
} | tee report.txt
