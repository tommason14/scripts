set res [lindex $argv 0]
mol addrep 0
mol modselect 0 0 not resname $res
mol modselect 1 0 resname $res
mol modstyle 1 0 VDW 1.000000 12.000000
