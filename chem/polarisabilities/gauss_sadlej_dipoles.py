from autochem import Settings, GaussJob
import os
from glob import glob

cwd = os.getcwd()

sett=Settings()
sett.input.method='M062X'
sett.input.basis='Gen'
sett.input.polar='Dipole' # print dipoles, might not be needed...

# basis set from https://www.basissetexchange.org/
# remove nitrogen basis functions if nitrogen isn't present otherwise run will fail
basis = [
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
"****\n",
"\n"
]

sett.meta.ncpus=48
sett.meta.mem='160gb' 
sett.meta.nodemem='192gb'
sett.meta.time='1:00:00'
sett.meta.partition='normal'
sett.meta.jobfs='200gb'

fields = {'plusX': 'X+8', 'minusX': 'X-8', 'plusY': 'Y+8', 'minusY': 'Y-8', 'plusZ': 'Z+8',
'minusZ': 'Z-8'}

for xyz in glob('*xyz'):
    xyzdir = xyz.replace('.xyz', '')
    os.mkdir(xyzdir)
    for name, field in fields.items():
        os.mkdir(f'{xyzdir}/{name}')
        os.chdir(f'{xyzdir}/{name}')
        _sett = sett.copy()
        _sett.input.field = field
        GaussJob(f'../../{xyz}', settings = _sett) # makes spec.job
        # add basis set
        with open('spec.job') as f:
            job = f.readlines()
        job = job[:-1] + basis + job[-1:]
        with open('spec.job', 'w') as new:
            for line in job:
                new.write(line)
        os.chdir(cwd)
