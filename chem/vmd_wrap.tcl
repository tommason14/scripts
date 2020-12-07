package require pbctools
set numFrames [molinfo top get numframes]
for  {set i 1} {$i <= $numFrames} {incr i} {
animate goto $i
pbc set [pbc get -now]
pbc wrap    
}
