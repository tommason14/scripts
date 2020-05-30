#!/usr/bin/env sh

replace_comma_with_ampersand(){
sed 's/,/ \& /g'
}

add_newline(){
sed 's/$/ \\\\/'
}

data="$1"
header=$(head -1 "$data")
echo $header
num_columns=$(echo $header | awk -F"," '{print NF}')
longtable=$(printf "|c%.0s" $(seq $num_columns) && printf "|")
formatted_header=$(echo $header | sed 's/\,/\ \&\ /g')

base=$(basename $data)

echo "\\begin{longtable}{$longtable}"
echo "\\caption{$base%.*}"
echo "\hline"
echo "$formatted_header \\"
echo "\hline"
echo "\endfirsthead"
echo "\hline"
echo "$formatted_header \\"
echo "\hline"
echo "\endhead"
echo "\hline \multicolumn{$num_columns}{r}{\textit{Continued on next page}} \\"
echo "\endfoot"
printf "  %s\n" $(tail -n +2 "$data") | replace_comma_with_ampersand | add_newline
# echo "\label{table:"$data"}"
echo "\end{table}"
