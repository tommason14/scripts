#!/usr/bin/env python3

from autochem import OrcaResults
import glob

log=glob.glob('*log')[0]

log = OrcaResults(log)

with open('uv_vis.csv', 'w') as f:
    f.write('Config,'
            'Root,'
            'Iteration,'
            'Oscillator Strength (eV),'
            'Wavelength (nm),'
            'Intensity (au)\n')
    for num, val in enumerate(
        zip(log.td_dft_wavelengths, 
            log.td_dft_intensities,
            log.td_dft_transition_energies), 
        1):
        w, i, o = val
        f.write(f'{log.log},{num},init,{o},{w},{i}\n')
