#!/usr/bin/env python3

bonds = {}
with open('connected.in') as f:
    for line in f:
        if 'bond' in line:
            _, one, two = line.split()
            one, two = int(one), int(two)
            if one not in bonds:
                bonds[one] = [two]
            else:
                bonds[one].append(two)

atoms = {}
found = False
num = 1
with open('polarisability.out') as p:
    for line in p:
        if 'Name' in line:
            found = True
            continue
        if 'Summed' in line:
            break
        if found:
            line = line.split()
            atoms[num] = {'element': line[0], 'pol': float(line[1])}
            num += 1

with open('summed_polarisabilities.out', 'w') as out:
    for number, bonded in bonds.items():
        pol = 0
        if atoms[number]['element'] != 'H':
            pol += atoms[number]['pol']
            for a in bonded:
                if atoms[a]['element'] == 'H':
                    pol += atoms[a]['pol']
            out.write(f'{atoms[number]["element"]}{number} = {pol:.2f}\n')
