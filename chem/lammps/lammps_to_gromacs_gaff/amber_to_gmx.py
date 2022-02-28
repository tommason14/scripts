#!/usr/bin/env python3
import parmed as pmd
p = pmd.load_file('polymer.prmtop', xyz='polymer.pdb')

pdb = pmd.load_file('polymer.pdb')
for old, new in zip(p.atoms, pdb.atoms):
    old.name = new.name

p.save('polymer.gro')
p.save('polymer.top')
