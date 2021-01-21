#!/usr/bin/env bash

autochem -ro energies.csv # collects all energies
# now remove fmo0/free-state-polarisation/full-pieda lines
head -1 energies.csv > fmo3.csv && grep 'cluster-[0-9]\+,' energies.csv >> fmo3.csv
