# Scripts

Scripts organised into various categories.

## Installation

Move this repo to `~/.local/scripts` and the following command to your shell rc file:
```
export PATH="$(find "$HOME/.local/scripts" -type d | grep -v "^.$\|.git\|pycache" | tr '\n' ':' | sed 's/:$//'):$PATH"
```

## Chem scripts 

- Faciltating LAMMPS, Gromacs, OpenMM file creation for classical and
  polarisable molecular dynamics 
- Analysis scripts, written to be file agnostic through the MDAnalysis
  and MDTraj python libraries, with specific scripts for desalination-type
  simulations 
- Scripts to process FMO PIEDA and Psi4 SAPT calculations, along with
  polarisabilities for polarisable MD
- quick plotting via gnuplot 
- various other utilities related to comp chem, for example converting
  GAMESS output into a molden-compatible format

## Supercomps 

- Scripts that work with the SLURM and PBS job schedulers

## Utils 

- Various shell scripts / python utility functions
- Unicode character selection via fzf (useful for Alacritty)
