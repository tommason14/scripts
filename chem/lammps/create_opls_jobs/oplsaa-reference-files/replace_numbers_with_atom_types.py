import re

with open('oplsaa.prm') as f:
    contents = f.readlines()

atoms = {}

for line in contents:
    if re.search('^\s*atom', line):
        line = line.split()
        atoms[line[2]] = line[3]

newfile = []
for line in contents:
    if re.search('^\s*bond', line):
        tmp = line.split()
        tmp[1] = atoms[tmp[1]]
        tmp[2] = atoms[tmp[2]]
        newfile.append('  '.join(tmp) + '\n')
    elif re.search('^\s*angle', line):
        tmp = line.split()
        tmp[1] = atoms[tmp[1]]
        tmp[2] = atoms[tmp[2]]
        tmp[3] = atoms[tmp[3]]
        newfile.append('  '.join(tmp) + '\n')
    elif re.search('^\s*torsion\s', line) or re.search('^\s*#torsion\s', line):
        tmp = line.split()
        # if missing, atom number = 0, replace with X
        tmp[1] = atoms.get(tmp[1], 'X')
        tmp[2] = atoms.get(tmp[2], 'X')
        tmp[3] = atoms.get(tmp[3], 'X')
        tmp[4] = atoms.get(tmp[4], 'X')
        newfile.append('  '.join(tmp) + '\n')
    else:
        newfile.append(line)

with open('modified_oplsaa.prm', 'w') as f:
    for line in newfile:
        f.write(line)
