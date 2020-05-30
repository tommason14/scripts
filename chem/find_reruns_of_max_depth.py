#!/usr/bin/env python3

import subprocess
import glob
import shutil
import os
import sys

if '-h' in sys.argv:
    print('To show which files are copied, add --dry to the call')

dry_run=False
if '--dry' in sys.argv:
    dry_run=True

configs = subprocess.check_output('find . -maxdepth 1 -type d | cut -d "/" -f 2', shell=True)
configs = configs.decode('utf-8').split('\n')[:-1]

rm = ('.', 'reruns')
for removal in rm:
    configs.remove(removal)

configs={config:[] for config in configs}

reruns = glob.glob('**/rerun.xyz', recursive=True)

for r in reruns: 
    conf = r.split('/')[0]
    configs[conf].append(r)

newdir=os.path.join(os.getcwd(), 'reruns')
# if not os.path.isdir(newdir):
#     os.mkdir(newdir)

#max depth = max length string
for k, v in configs.items():
    rerun = max(v)
    newname = f'{newdir}/{k}.xyz'
    if dry_run:
        print('-'*60)
        print(f'{rerun}\n↳ {newname}')
    else:
        # print('copied')
        shutil.copy(rerun, newname)

