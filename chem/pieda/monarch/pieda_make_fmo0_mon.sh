#!/usr/bin/env bash

files="$(find . -name "spec.log")"

cwd="$(pwd)"
for f in $files
do
  echo "$f"
  base="$(basename "$f")"
  base="${base%.*}" # lose extension
  dir="$(dirname "$f")"
  cd "$dir"
  # Max AOs per n-mer:      560     1120     1575
  # length of line can vary (FMO2 v FMO3), so split on : then on whitespace
  max_ao="$(grep "Max AOs per" "$base.log")"
  max_ao="$(echo "${max_ao#*:}" | awk '{print $1}')"
  [[ ! -d "fmo0" ]]  && mkdir "fmo0"
  cp "$base.inp" "$base.job" "fmo0"
  # change runtyp, lose NBODY=...
  sed -i 's/RUNTYP=ENERGY/RUNTYP=FMO0/;s/NBODY=[0-9]//' "fmo0/$base.inp"
  # change fmoprp line- assumes one is already present; if not this script will fail!
  sed -i "s/\$FMOPRP.*/\$FMOPRP MAXAOC=$max_ao MODORB=3 MAXIT=200 \$END/" "fmo0/$base.inp"
  # remove $GDDI line
  sed -i "/\$GDDI/d" "fmo0/$base.inp"
  # one node, normal queue
  # 1 hr walltime
  sed -i "s/\#SBATCH --time=.*/\#SBATCH --time=1:00:00/" "fmo0/$base.job"
  # retain temporary files!
  sed -i "s|rungms.monarch|~/gamess16-srs/rungms.monarch.keep_files|" "fmo0/$base.job"

  cd "$cwd"
done
