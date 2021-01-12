#!/usr/bin/env bash                                                                                
                                                                                                   
# using Ovito fails with drudes, but vmd works
                                                                                                   
[[ $# -eq 0 ]] && echo "Syntax: $(basename $0) data.lmps" && exit 1                                  
                                                                                                   
data="$1"                                                                                          
[[ ! -f "$data" ]] && echo "Error: $data doesn't exist" && exit 1                                    
[[ $# -gt 1 ]] && output="$2" || output=$(echo $1 | awk -F"." '{OFS="."; $NF="unwrapped.xyz"; print $0}')

if [[ "$USER" = "tommason" ]]
then                                                                                                  
  vmd="/Applications/VMD 1.9.4a38.app/Contents/vmd/vmd_MACOSXX86_64"                                  
  sed='gsed'
elif [[ "$USER" = "tmas0023" ]]
then                                                                                                  
  vmd="/Applications/VMD 1.9.3.app/Contents/vmd/vmd_MACOSXX86"                                        
  sed='gsed'
else                                                                                                  
  vmd="vmd" # module loaded                                                                           
  sed='sed'
fi                                                                                                    
                                                                                                      
[[ -f tmpvmd ]] && rm tmpvmd                                                                            
cat << EOF >> tmpvmd
package require topotools
topo readlammpsdata "$data"
set sel [atomselect top all]
\$sel writexyz "$output"
exit
EOF
                                                                                                      
"$vmd" -e tmpvmd -dispdev none
rm tmpvmd

atoms=($(sed -n '/Masses/,/^\s*[A-Z]/p' "$data" | grep '^\s*[0-9]' | awk '{if($NF=="DP"){print $1,"DP"} else {print $1,$4}}'))

while read num element
  do 
    $sed -i "s/^\s\+$num/$element/" "$output"
  done < <(printf "%s %s\n" $(echo ${atoms[@]} | xargs -n2))
