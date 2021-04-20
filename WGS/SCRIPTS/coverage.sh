reference="$1"

consensus="$2"

sed -e '/>/! s/[nN]//g' "$consensus" > "$consensus".tmp

sed -e '/>/! s/[nN]//g' "$reference" > "$reference".tmp

ref=`fastalength "$reference".tmp | awk '{print $1}'`

fastalength "$consensus".tmp | awk '{print $2, $1/'"$ref"'*100}'

rm "$consensus".tmp "$reference".tmp