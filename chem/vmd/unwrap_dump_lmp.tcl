package require pbctools
pbc unwrap

set nframes [molinfo top get numframes]
set all [atomselect top all]
set natoms [$all num]

set out [open "dump.lmp" w]
for {set i 0} {$i < $nframes} {incr i} {
animate goto $i
puts $out "ITEM: TIMESTEP"
puts $out "$i"
puts $out "ITEM: NUMBER OF ATOMS"
puts $out "$natoms"
puts $out "ITEM: BOX BOUNDS pp pp pp"
lassign [measure minmax $all] min max
lassign $min xlo ylo zlo
lassign $max xhi yhi zhi
puts $out [format "%.6f %.6f" $xlo $xhi]
puts $out [format "%.6f %.6f" $ylo $yhi]
puts $out [format "%.6f %.6f" $zlo $zhi]
puts $out "ITEM: ATOMS element xu yu zu"
foreach adat [$all get {name x y z}] {
  lassign $adat name x y z
  puts $out [format "%s %.6f %.6f %.6f" $name $x $y $z]
}
}
close $out
exit
