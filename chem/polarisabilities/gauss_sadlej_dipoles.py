#!/usr/bin/env python3
from autochem import Settings, GaussJob, Molecule
import os
from glob import glob

cwd = os.getcwd()

sett = Settings()
sett.input.method = 'M062X'
sett.input.basis = 'Gen'
sett.input.polar = 'Dipole'  # print dipoles, might not be needed...
sett.meta.ncpus = 48
sett.meta.mem = '160gb'
sett.meta.nodemem = '192gb'
sett.meta.time = '1:00:00'
sett.meta.partition = 'normal'
sett.meta.jobfs = '200gb'

# basis sets from https://www.basissetexchange.org/

bsets = {}
bsets['CHO'] = [
    "! Sadlej-pVTZ for H,C,O\n"
    "H     0\n"
    "S    4   1.00\n"
    "     33.685014                .006068\n"
    "      5.094788                .045316\n"
    "      1.158786                .202846\n"
    "      0.325840                .503709\n"
    "S    1   1.00\n"
    "      0.102741               1.00000\n"
    "S    1   1.00\n"
    "      0.032400               1.00000\n"
    "P    2   1.00\n"
    "      1.1588                  .188440\n"
    "      0.3258                  .882420\n"
    "P    2   1.00\n"
    "      0.1027                  .117800\n"
    "      0.0324                  .004200\n"
    "****\n"
    "C     0\n"
    "S    5   1.00\n"
    "   5240.6353                  .000937\n"
    "    782.20480                 .007228\n"
    "    178.35083                 .036344\n"
    "     50.815942                .130600\n"
    "     16.823562                .318931\n"
    "S    2   1.00\n"
    "      6.175776                .438742\n"
    "      2.418049                .214974\n"
    "S    1   1.00\n"
    "      0.511900               1.\n"
    "S    1   1.00\n"
    "      0.156590               1.\n"
    "S    1   1.00\n"
    "      0.047900               1.\n"
    "P    4   1.00\n"
    "     18.841800                .013887\n"
    "      4.159240                .086279\n"
    "      1.206710                .288744\n"
    "      0.385540                .499411\n"
    "P    1   1.00\n"
    "      0.121940               1.\n"
    "P    1   1.00\n"
    "      0.038568               1.\n"
    "D    2   1.00\n"
    "      1.206710                .262852\n"
    "      0.385540                .804300\n"
    "D    2   1.00\n"
    "      0.121940                .653500\n"
    "      0.038658                .863600\n"
    "****\n"
    "O     0\n"
    "S    5   1.00\n"
    "  10662.285                   .000799\n"
    "   1599.7097                  .006153\n"
    "    364.72526                 .031157\n"
    "    103.65179                 .115596\n"
    "     33.905805                .301552\n"
    "S    2   1.00\n"
    "     12.287469                .444870\n"
    "      4.756805                .243172\n"
    "S    1   1.00\n"
    "      1.004271               1.\n"
    "S    1   1.00\n"
    "      0.300686               1.\n"
    "S    1   1.00\n"
    "      0.09003                1.\n"
    "P    4   1.00\n"
    "     34.856463                .015648\n"
    "      7.843131                .098197\n"
    "      2.306249                .307768\n"
    "      0.723164                .492470\n"
    "P    1   1.00\n"
    "      0.214882               1.\n"
    "P    1   1.00\n"
    "      0.06385                1.\n"
    "D    2   1.00\n"
    "      2.3062                  .20270\n"
    "      0.7232                  .5791\n"
    "D    2   1.00\n"
    "      0.2149                  .78545\n"
    "      0.0639                  .53387\n"
    "****\n"
]
bsets['CHNO'] = [
    "! Sadlej-pVTZ for H,C,N,O\n"
    "H     0\n"
    "S    4   1.00\n"
    "     33.685014                .006068\n"
    "      5.094788                .045316\n"
    "      1.158786                .202846\n"
    "      0.325840                .503709\n"
    "S    1   1.00\n"
    "      0.102741               1.00000\n"
    "S    1   1.00\n"
    "      0.032400               1.00000\n"
    "P    2   1.00\n"
    "      1.1588                  .188440\n"
    "      0.3258                  .882420\n"
    "P    2   1.00\n"
    "      0.1027                  .117800\n"
    "      0.0324                  .004200\n"
    "****\n"
    "C     0\n"
    "S    5   1.00\n"
    "   5240.6353                  .000937\n"
    "    782.20480                 .007228\n"
    "    178.35083                 .036344\n"
    "     50.815942                .130600\n"
    "     16.823562                .318931\n"
    "S    2   1.00\n"
    "      6.175776                .438742\n"
    "      2.418049                .214974\n"
    "S    1   1.00\n"
    "      0.511900               1.\n"
    "S    1   1.00\n"
    "      0.156590               1.\n"
    "S    1   1.00\n"
    "      0.047900               1.\n"
    "P    4   1.00\n"
    "     18.841800                .013887\n"
    "      4.159240                .086279\n"
    "      1.206710                .288744\n"
    "      0.385540                .499411\n"
    "P    1   1.00\n"
    "      0.121940               1.\n"
    "P    1   1.00\n"
    "      0.038568               1.\n"
    "D    2   1.00\n"
    "      1.206710                .262852\n"
    "      0.385540                .804300\n"
    "D    2   1.00\n"
    "      0.121940                .653500\n"
    "      0.038658                .863600\n"
    "****\n"
    "N     0\n"
    "S    5   1.00\n"
    "   8104.0716                  .000802\n"
    "   1216.0215                  .006174\n"
    "    277.23428                 .031233\n"
    "     76.904023                .115198\n"
    "     25.874419                .296951\n"
    "S    2   1.00\n"
    "      9.346767                .447349\n"
    "      3.579794                .245003\n"
    "S    1   1.00\n"
    "      0.739610               1.\n"
    "S    1   1.00\n"
    "      0.222617               1.\n"
    "S    1   1.00\n"
    "      0.067006               1.\n"
    "P    4   1.00\n"
    "     26.868987                .014478\n"
    "      5.991227                .091156\n"
    "      1.750842                .297420\n"
    "      0.560511                .493796\n"
    "P    1   1.00\n"
    "      0.175948               1.\n"
    "P    1   1.00\n"
    "      0.055231               1.\n"
    "D    2   1.00\n"
    "      1.750842                .224774\n"
    "      0.560511                .659562\n"
    "D    2   1.00\n"
    "      0.175948                .871355\n"
    "      0.055231                .704217\n"
    "****\n"
    "O     0\n"
    "S    5   1.00\n"
    "  10662.285                   .000799\n"
    "   1599.7097                  .006153\n"
    "    364.72526                 .031157\n"
    "    103.65179                 .115596\n"
    "     33.905805                .301552\n"
    "S    2   1.00\n"
    "     12.287469                .444870\n"
    "      4.756805                .243172\n"
    "S    1   1.00\n"
    "      1.004271               1.\n"
    "S    1   1.00\n"
    "      0.300686               1.\n"
    "S    1   1.00\n"
    "      0.09003                1.\n"
    "P    4   1.00\n"
    "     34.856463                .015648\n"
    "      7.843131                .098197\n"
    "      2.306249                .307768\n"
    "      0.723164                .492470\n"
    "P    1   1.00\n"
    "      0.214882               1.\n"
    "P    1   1.00\n"
    "      0.06385                1.\n"
    "D    2   1.00\n"
    "      2.3062                  .20270\n"
    "      0.7232                  .5791\n"
    "D    2   1.00\n"
    "      0.2149                  .78545\n"
    "      0.0639                  .53387\n"
    "****\n"
]
bsets['CHOS'] = [
    "! Sadlej-pVTZ for H,C,O,S\n"
    "H     0\n"
    "S    4   1.00\n"
    "     33.685014                .006068\n"
    "      5.094788                .045316\n"
    "      1.158786                .202846\n"
    "      0.325840                .503709\n"
    "S    1   1.00\n"
    "      0.102741               1.00000\n"
    "S    1   1.00\n"
    "      0.032400               1.00000\n"
    "P    2   1.00\n"
    "      1.1588                  .188440\n"
    "      0.3258                  .882420\n"
    "P    2   1.00\n"
    "      0.1027                  .117800\n"
    "      0.0324                  .004200\n"
    "****\n"
    "C     0\n"
    "S    5   1.00\n"
    "   5240.6353                  .000937\n"
    "    782.20480                 .007228\n"
    "    178.35083                 .036344\n"
    "     50.815942                .130600\n"
    "     16.823562                .318931\n"
    "S    2   1.00\n"
    "      6.175776                .438742\n"
    "      2.418049                .214974\n"
    "S    1   1.00\n"
    "      0.511900               1.\n"
    "S    1   1.00\n"
    "      0.156590               1.\n"
    "S    1   1.00\n"
    "      0.047900               1.\n"
    "P    4   1.00\n"
    "     18.841800                .013887\n"
    "      4.159240                .086279\n"
    "      1.206710                .288744\n"
    "      0.385540                .499411\n"
    "P    1   1.00\n"
    "      0.121940               1.\n"
    "P    1   1.00\n"
    "      0.038568               1.\n"
    "D    2   1.00\n"
    "      1.206710                .262852\n"
    "      0.385540                .804300\n"
    "D    2   1.00\n"
    "      0.121940                .653500\n"
    "      0.038658                .863600\n"
    "****\n"
    "O     0\n"
    "S    5   1.00\n"
    "  10662.285                   .000799\n"
    "   1599.7097                  .006153\n"
    "    364.72526                 .031157\n"
    "    103.65179                 .115596\n"
    "     33.905805                .301552\n"
    "S    2   1.00\n"
    "     12.287469                .444870\n"
    "      4.756805                .243172\n"
    "S    1   1.00\n"
    "      1.004271               1.\n"
    "S    1   1.00\n"
    "      0.300686               1.\n"
    "S    1   1.00\n"
    "      0.09003                1.\n"
    "P    4   1.00\n"
    "     34.856463                .015648\n"
    "      7.843131                .098197\n"
    "      2.306249                .307768\n"
    "      0.723164                .492470\n"
    "P    1   1.00\n"
    "      0.214882               1.\n"
    "P    1   1.00\n"
    "      0.06385                1.\n"
    "D    2   1.00\n"
    "      2.3062                  .20270\n"
    "      0.7232                  .5791\n"
    "D    2   1.00\n"
    "      0.2149                  .78545\n"
    "      0.0639                  .53387\n"
    "****\n"
    "S     0\n"
    "S    6   1.00\n"
    "  99413.4                     .000742\n"
    "  13961.7                     .005790\n"
    "   3169.9                     .029941\n"
    "    902.46                    .118975\n"
    "    297.16                    .368290\n"
    "    108.702                   .577489\n"
    "S    3   1.00\n"
    "    108.702                   .142932\n"
    "     43.155                   .624649\n"
    "     18.108                   .283400\n"
    "S    1   1.00\n"
    "      5.5705                 1.\n"
    "S    1   1.00\n"
    "      2.1427                 1.\n"
    "S    1   1.00\n"
    "      0.4340                 1.\n"
    "S    1   1.00\n"
    "      0.1570                 1.\n"
    "S    1   1.00\n"
    "      0.0568                 1.\n"
    "P    6   1.00\n"
    "    495.04                    .003171\n"
    "    117.22                    .024643\n"
    "     37.507                   .107845\n"
    "     13.910                   .288489\n"
    "      5.5045                  .448108\n"
    "      2.2433                  .320517\n"
    "P    1   1.00\n"
    "      0.7762                 1.\n"
    "P    1   1.00\n"
    "      0.2919                 1.\n"
    "P    1   1.00\n"
    "      0.1029                 1.\n"
    "P    1   1.00\n"
    "      0.0363                 1.\n"
    "D    2   1.00\n"
    "      0.7762                  .2915\n"
    "      0.2919                  .8966\n"
    "D    2   1.00\n"
    "      0.1029                  .1138\n"
    "      0.0363                  .0568\n"
    "****\n"
]
bsets['HOP'] = [
    "! Sadlej-pVTZ for H,O,P\n"
    "H     0\n"
    "S    4   1.00\n"
    "     33.685014                .006068\n"
    "      5.094788                .045316\n"
    "      1.158786                .202846\n"
    "      0.325840                .503709\n"
    "S    1   1.00\n"
    "      0.102741               1.00000\n"
    "S    1   1.00\n"
    "      0.032400               1.00000\n"
    "P    2   1.00\n"
    "      1.1588                  .188440\n"
    "      0.3258                  .882420\n"
    "P    2   1.00\n"
    "      0.1027                  .117800\n"
    "      0.0324                  .004200\n"
    "****\n"
    "O     0\n"
    "S    5   1.00\n"
    "  10662.285                   .000799\n"
    "   1599.7097                  .006153\n"
    "    364.72526                 .031157\n"
    "    103.65179                 .115596\n"
    "     33.905805                .301552\n"
    "S    2   1.00\n"
    "     12.287469                .444870\n"
    "      4.756805                .243172\n"
    "S    1   1.00\n"
    "      1.004271               1.\n"
    "S    1   1.00\n"
    "      0.300686               1.\n"
    "S    1   1.00\n"
    "      0.09003                1.\n"
    "P    4   1.00\n"
    "     34.856463                .015648\n"
    "      7.843131                .098197\n"
    "      2.306249                .307768\n"
    "      0.723164                .492470\n"
    "P    1   1.00\n"
    "      0.214882               1.\n"
    "P    1   1.00\n"
    "      0.06385                1.\n"
    "D    2   1.00\n"
    "      2.3062                  .20270\n"
    "      0.7232                  .5791\n"
    "D    2   1.00\n"
    "      0.2149                  .78545\n"
    "      0.0639                  .53387\n"
    "****\n"
    "P     0\n"
    "S    6   1.00\n"
    "  77492.43                    .000787\n"
    "  11605.79                    .006109\n"
    "   2645.96                    .031368\n"
    "    754.98                    .124291\n"
    "    248.75                    .380682\n"
    "     91.157                   .559989\n"
    "S    3   1.00\n"
    "     91.157                   .163997\n"
    "     36.226                   .625927\n"
    "     15.211                   .262211\n"
    "S    1   1.00\n"
    "      4.7138                 1.\n"
    "S    1   1.00\n"
    "      1.7827                 1.\n"
    "S    1   1.00\n"
    "      0.3425                 1.\n"
    "S    1   1.00\n"
    "      0.1246                 1.\n"
    "S    1   1.00\n"
    "      0.0453                 1.\n"
    "P    6   1.00\n"
    "    384.84                    .003765\n"
    "     90.552                   .029010\n"
    "     28.806                   .122953\n"
    "     10.688                   .306470\n"
    "      4.2521                  .441349\n"
    "      1.7405                  .295817\n"
    "P    1   1.00\n"
    "      0.5979                 1.\n"
    "P    1   1.00\n"
    "      0.2292                 1.\n"
    "P    1   1.00\n"
    "      0.0838                 1.\n"
    "P    1   1.00\n"
    "      0.0306                 1.\n"
    "D    2   1.00\n"
    "      0.5979                  .3093\n"
    "      0.2292                  .9715\n"
    "D    2   1.00\n"
    "      0.0838                  .1278\n"
    "      0.0306                  .0774\n"
    "****\n"
]

bsets['CHOP'] = [
    "! Sadlej-pVTZ for H,C,O,P\n"
    "H     0\n"
    "S    4   1.00\n"
    "     33.685014                .006068\n"
    "      5.094788                .045316\n"
    "      1.158786                .202846\n"
    "      0.325840                .503709\n"
    "S    1   1.00\n"
    "      0.102741               1.00000\n"
    "S    1   1.00\n"
    "      0.032400               1.00000\n"
    "P    2   1.00\n"
    "      1.1588                  .188440\n"
    "      0.3258                  .882420\n"
    "P    2   1.00\n"
    "      0.1027                  .117800\n"
    "      0.0324                  .004200\n"
    "****\n"
    "C     0\n"
    "S    5   1.00\n"
    "   5240.6353                  .000937\n"
    "    782.20480                 .007228\n"
    "    178.35083                 .036344\n"
    "     50.815942                .130600\n"
    "     16.823562                .318931\n"
    "S    2   1.00\n"
    "      6.175776                .438742\n"
    "      2.418049                .214974\n"
    "S    1   1.00\n"
    "      0.511900               1.\n"
    "S    1   1.00\n"
    "      0.156590               1.\n"
    "S    1   1.00\n"
    "      0.047900               1.\n"
    "P    4   1.00\n"
    "     18.841800                .013887\n"
    "      4.159240                .086279\n"
    "      1.206710                .288744\n"
    "      0.385540                .499411\n"
    "P    1   1.00\n"
    "      0.121940               1.\n"
    "P    1   1.00\n"
    "      0.038568               1.\n"
    "D    2   1.00\n"
    "      1.206710                .262852\n"
    "      0.385540                .804300\n"
    "D    2   1.00\n"
    "      0.121940                .653500\n"
    "      0.038658                .863600\n"
    "****\n"
    "O     0\n"
    "S    5   1.00\n"
    "  10662.285                   .000799\n"
    "   1599.7097                  .006153\n"
    "    364.72526                 .031157\n"
    "    103.65179                 .115596\n"
    "     33.905805                .301552\n"
    "S    2   1.00\n"
    "     12.287469                .444870\n"
    "      4.756805                .243172\n"
    "S    1   1.00\n"
    "      1.004271               1.\n"
    "S    1   1.00\n"
    "      0.300686               1.\n"
    "S    1   1.00\n"
    "      0.09003                1.\n"
    "P    4   1.00\n"
    "     34.856463                .015648\n"
    "      7.843131                .098197\n"
    "      2.306249                .307768\n"
    "      0.723164                .492470\n"
    "P    1   1.00\n"
    "      0.214882               1.\n"
    "P    1   1.00\n"
    "      0.06385                1.\n"
    "D    2   1.00\n"
    "      2.3062                  .20270\n"
    "      0.7232                  .5791\n"
    "D    2   1.00\n"
    "      0.2149                  .78545\n"
    "      0.0639                  .53387\n"
    "****\n"
    "P     0\n"
    "S    6   1.00\n"
    "  77492.43                    .000787\n"
    "  11605.79                    .006109\n"
    "   2645.96                    .031368\n"
    "    754.98                    .124291\n"
    "    248.75                    .380682\n"
    "     91.157                   .559989\n"
    "S    3   1.00\n"
    "     91.157                   .163997\n"
    "     36.226                   .625927\n"
    "     15.211                   .262211\n"
    "S    1   1.00\n"
    "      4.7138                 1.\n"
    "S    1   1.00\n"
    "      1.7827                 1.\n"
    "S    1   1.00\n"
    "      0.3425                 1.\n"
    "S    1   1.00\n"
    "      0.1246                 1.\n"
    "S    1   1.00\n"
    "      0.0453                 1.\n"
    "P    6   1.00\n"
    "    384.84                    .003765\n"
    "     90.552                   .029010\n"
    "     28.806                   .122953\n"
    "     10.688                   .306470\n"
    "      4.2521                  .441349\n"
    "      1.7405                  .295817\n"
    "P    1   1.00\n"
    "      0.5979                 1.\n"
    "P    1   1.00\n"
    "      0.2292                 1.\n"
    "P    1   1.00\n"
    "      0.0838                 1.\n"
    "P    1   1.00\n"
    "      0.0306                 1.\n"
    "D    2   1.00\n"
    "      0.5979                  .3093\n"
    "      0.2292                  .9715\n"
    "D    2   1.00\n"
    "      0.0838                  .1278\n"
    "      0.0306                  .0774\n"
    "****\n"
]

fields = {
    'plusX': 'X+8',
    'minusX': 'X-8',
    'plusY': 'Y+8',
    'minusY': 'Y-8',
    'plusZ': 'Z+8',
    'minusZ': 'Z-8'
}

for xyz in glob('*xyz'):
    xyzdir = xyz.replace('.xyz', '')
    mol = Molecule(xyz)
    atoms = list(set([a.symbol for a in mol.coords]))
    try:
        basis = bsets[''.join(sorted(atoms))] + ['\n']
    except ValueError:
        print(
            'Download the required basis set from basissetexchange.org and add to gauss_sadlej_dipoles.py'
        )
    os.mkdir(xyzdir)
    for name, field in fields.items():
        os.mkdir(f'{xyzdir}/{name}')
        os.chdir(f'{xyzdir}/{name}')
        _sett = sett.copy()
        _sett.input.field = field
        GaussJob(f'../../{xyz}', settings=_sett)  # makes spec.job
        # add basis set
        with open('spec.job') as f:
            job = f.readlines()
        job = job[:-1] + basis + job[-1:]
        with open('spec.job', 'w') as new:
            for line in job:
                new.write(line)
        os.chdir(cwd)
