package require pbctools

puts "Unwrapping..."
pbc join residue -first 0 -last 0
pbc unwrap
pbc wrap -center com -centersel "resname pol" -compound residue -all

proc minmax_z {sel} {
 set coords [$sel get {x y z}]
  set coord [lvarpop coords]
  lassign $coord minx miny minz
  lassign $coord maxx maxy maxz
  foreach coord $coords {
    lassign $coord x y z
    if {$x < $minx} {set minx $x} else {if {$x > $maxx} {set maxx $x}}
    if {$y < $miny} {set miny $y} else {if {$y > $maxy} {set maxy $y}}
    if {$z < $minz} {set minz $z} else {if {$z > $maxz} {set maxz $z}}
  }
  # return [list [list $minx $miny $minz] [list $maxx $maxy $maxz]]
  return [list $minz $maxz]
}

set nframes [molinfo top get numframes]
# polymer dimensions are fixed - so compute once
set pol [atomselect top "resname pol"]
lassign [minmax_z $pol] minz maxz

puts "Counting ions..."

set file [open "ion_counts.csv" "w"]
puts $file "frame,cl_salt,cl_mem,cl_fresh,na_salt,na_mem,na_fresh"
for {set i 0} {$i < $nframes} {incr i} {
  animate goto $i
  set cl_in_saltwater [[atomselect top "resname CL and z < $minz"] num]
  set cl_in_membrane [[atomselect top "resname CL and z > $minz and z < $maxz"] num]
  set cl_in_freshwater [[atomselect top "resname CL and z > $maxz"] num]
  set na_in_saltwater [[atomselect top "resname NA and z < $minz"] num]
  set na_in_membrane [[atomselect top "resname NA and z > $minz and z < $maxz"] num]
  set na_in_freshwater [[atomselect top "resname NA and z > $maxz"] num]
  puts $file "$i,$cl_in_saltwater,$cl_in_membrane,$cl_in_freshwater,$na_in_saltwater,$na_in_membrane,$na_in_freshwater"
  # Print out percentage completion of nframes every 100 frames
  if {$i % 100 == 0} {
    set percent [expr {100*double($i)/$nframes}]
    puts [format "frame $i of $nframes (%.2f %%)" $percent]
  }
}
exit
