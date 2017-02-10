# Purpose: Create subdirectories and generate input pdbs with varying dihedral angle.
# Usage: "vmd -dispdev none -e setupDihedrals.tcl"

# Before using, check the bottom to:
#    - edit home and final directories and PDB file names
#    - edit atom numbers, subtracting one (VMD starts with 0)
# To do: make the above ^ command line arguments

# Adapted from:
# https://github.com/Eigenstate/vmd-python/blob/master/vmd/plugins/molefacture/molefacture_builder.tcl

proc set_dihedral {molid movesel newval ind1 ind2 ind3 ind4} {
  set tmpmolid $molid
  set dihedral [measure dihed [list [list $ind1 $tmpmolid] [list $ind2 $tmpmolid] [list $ind3 $tmpmolid] [list $ind4 $tmpmolid] ]]

  ### Validate input.
  if {$ind1 == "" || $ind2 == "" || $ind3 == "" || $ind4 == ""} {
    puts "WARNING: Couldn't find one or more atoms for set_dihedral"
    return
  }

  ### Set both atoms of the central bond.
  set bsel1 [atomselect $tmpmolid "index $ind2"]
  set bsel2 [atomselect $tmpmolid "index $ind3"]

  ### Get degree to move by, and the move matrix.
  set delta [expr $dihedral - $newval]
  set mat [trans bond [lindex [$bsel1 get {x y z}] 0] [lindex [$bsel2 get {x y z}] 0] $delta deg]

  ### Move the selection by that matrix.
  $movesel move $mat
  $bsel1 delete
  $bsel2 delete

  ### User information.
  set dihedral [measure dihed [list [list $ind1 $tmpmolid] [list $ind2 $tmpmolid] [list $ind3 $tmpmolid] [list $ind4 $tmpmolid] ]]
  puts $dihedral
}

# ======================================================================

mol new "Butane.pdb" type {pdb} first 0 last -1 step 1 waitfor -1
set butane [atomselect 0 "index 0 1 2 3 4 16 17"]

for {set angle 0} {$angle < 360} {incr angle 1} {
    set_dihedral 0 $guanidino $angle 1 5 6 7
    mkdir ../angles-qm/$angle
    [atomselect top "all"] writepdb ../angles-qm/$angle/gbi2-$angle.pdb
}
