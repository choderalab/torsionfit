package require psfgen
topology ../data/charmm_ff/top_all36_cgenff.rtf

mol load pdb butane.pdb
set BUTA [atomselect top all]

segment BUTA {
pdb butane.pdb
first NONE
last NONE
}
coordpdb butane.pdb BUTA

regenerate dihedrals angles
guesscoord

writepsf butane.psf
writepdb butane_new.pdb
