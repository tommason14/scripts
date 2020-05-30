from get_coeffs_gaff import getAtomData, getBond, getImproper
from get_coeffs_gaff import getAngle, getDihedral
from get_coeffs_gaff import getAtomPartialCharge
import subprocess as sp
import glob, re, sys

""" Substitute in different labels to *-l.data files and create
topology file with VMD topo. After VMD original names are
substituted back in and coefficients, box size and atom data are
added. This includes Pair Coeffs, Bond Coeffs, Angle Coeffs,
Dihedral Coeffs, Improper Coeffs and Atoms (partial charges)
sections and box dimensions - xlo, xhi, ylo, yhi, zlo, zhi. """

# PATTERN MATCHING
pattern_labelled_xyz = '^\s*[A-Z]*[a-z]?[0-9]*\s+-?[0-9]\.'
pattern_coeffs       = '^# [0-9]+\s+[A-Z]{1,2}[a-z]?[0-9]*'

# REWRITE -l XYZ FILES FOR NAME CHANGES
Files = glob.glob("*-l.xyz")

if len(Files) == 0:
    sys.exit('Cannot find any files for "*-l.xyz".')

# EXPECT SECOND CHARACTER IS LOWER CASE IF PART OF ELEMENT SYMBOL
for File in Files:

    ### SUBSTITUTE IN ALT NAMES TO -l.xyz ----------------------------

    Name = File.replace('-l.xyz', '')

    # NUMBER OF NEW ELEMENTS
    count = 0

    # NEW ELEMENT NAMES
    elemDict = {}

    # NEW LINES
    newLines = []

    lines = open(File, 'r+').readlines()
    for line_no, line in enumerate(lines):

        # MAKE SURE PAST FIRST LINE
        if line_no > 1:

            # START OF LABELLED XYZ
            if re.search(pattern_labelled_xyz, line):

                # CUT AFTER SYMBOL NAME
                sym, coord = line.strip().split(" ", 1)

                # TURN sym INTO A LIST
                list_sym = list(sym)

                # CHECK IF SECOND LETTER OF ELEMENT IN LOWER CASE
                try:
                    if list_sym[1].islower():

                        # BOTH LETTERS ARE ELEMENT
                        newID = list_sym[:1] + f"{count:03d}"

                    else:
                        # ONLY FIRST LETTER IS ELEMENT
                        newID = list_sym[0]  + f"{count:03d}"
                except:
                    # ONLY FIRST LETTER IS ELEMENT
                    newID = list_sym[0]  + f"{count:03d}"


                # CHECK IF COUNT IN DICTIONARY
                if sym in elemDict.keys():
                    # GET newID VALUE
                    newID = elemDict[sym]

                else:
                    # PUT CHANGE IN DICTIONARY {'C1' : 'CZ'}
                    elemDict[sym] = newID
                    count += 1

                # ADD NEW LINE WITH NEW SYMBOLS
                newLines.append(newID + "   " + coord + '\n')

            else:
                newLines.append(line)

        else:
            newLines.append(line)

    open('topo-in.xyz', 'w+').writelines(newLines)

    ##### RUN VMD TOPO ----------------------------------------

    # COMMANDS FOR VMD
    lines ='package require topotools\n'        +\
            'topo retypebonds\n'                +\
            'topo guessangles\n'                +\
            'topo guessdihedrals\n'             +\
            'topo guessimpropers\n'             +\
            'topo writelammpsdata topo.out\n'   +\
            'exit'

    open('tempfile', 'w+').write(lines)

    # OPEN VMD WITH EDITTED XYZ AND COMMAND FILE
    cmd = "vmd -dispdev none -m topo-in.xyz -e tempfile"
    sp.check_output(cmd, shell=True)

    ### SUBSTITUTE IN ORIG NAMES ------------------------------

    newLines = []

    # SWITCH VALUES AND KEYS IN DICT FOR EASY ACCESS
    elemDict = {y: x for x, y in elemDict.items()}

    # FOR EACH LINE IN TOPO FILE
    lines = open('topo.out', 'r+').readlines()
    for line in lines:

        # FOR EACH NEW SYMBOL TO REPLACE
        for newID in elemDict.keys():

            # IF SYMBOL AFTER A HASH
            if re.search('#.*'+newID, line):
                # GET ORIG NAME
                sym = elemDict[newID]

                # REPLACE NEW WITH OLD STRING
                line = line.replace(newID, sym)

        newLines.append(line)

    ### EDIT DATA FILE ----------------------------------------

    # READ LINES
    lines = newLines[:]

    # SWITCHES - TURN ON WHEN FIND IN FILE
    pcoef = False   # Pair Coeffs
    bcoef = False   # Bond Coeffs
    acoef = False   # Angle Coeffs
    dcoef = False   # Dihedral Coeffs
    icoef = False   # Improper Coeffs
    atoms = False   # Atoms

    # SAVE X, Y, Z VALUES
    xvals = []
    yvals = []
    zvals = []

    # SAVE PARTIAL CHARGES
    pcharges = []

    # START EDITTING LINES
    newLines = []

    for line in lines:

        # PAIR COEFFS SWITCH
        if 'Pair Coeffs' in line:
            pcoef = True
            newLines.append('Pair Coeffs\n')

        # BOND COEFFS SWITCH
        elif 'Bond Coeffs' in line:
            pcoef = False
            bcoef = True
            newLines.append('Bond Coeffs\n')

        # ANGLE COEFFS SWITCH
        elif 'Angle Coeffs' in line:
            bcoef = False
            acoef = True
            newLines.append('Angle Coeffs\n')

        # DIHEDRAL COEFFS SWITCH
        elif 'Dihedral Coeffs' in line:
            acoef = False
            dcoef = True
            newLines.append('Dihedral Coeffs\n')

        # IMPROPER COEFFS SWITCH
        elif 'Improper Coeffs' in line:
            dcoef = False
            icoef = True
            newLines.append('Improper Coeffs\n')

        # ATOMS PARTIAL CHARGES SWITCH
        elif 'Atoms' in line:
            dcoef = False
            icoef = False
            atoms = True
            newLines.append(line)

        # ADD PAIR COEFFS
        elif pcoef and re.search(pattern_coeffs, line):
            h, num, name = line.split()
            eps, sigma   = getAtomData(name)
            # for lj/charmm/coul/long
            newline      = '{:4} {:>10.3f} {:>9.5f}    # {}\n'.format(num, eps, sigma, name)
            newLines.append(newline)

        # ADD BOND COEFFS
        elif bcoef and re.search(pattern_coeffs, line):
            h, num, atms = line.split()
            atm1, atm2   = atms.split('-')
            kr, Re       = getBond(atm1, atm2)
            newline      = '{:4} {:>9.1f} {:>10.3f}    # {}\n'.format(num, kr, Re, atms)
            newLines.append(newline)

        # ADD ANGLE COEFFS
        elif acoef and re.search(pattern_coeffs, line):
            h, num, atms     = line.split()
            atm1, atm2, atm3 = atms.split('-')
            ka, th           = getAngle(atm1, atm2, atm3)
            newline          = '{:4} {:>9.2f} {:>10.2f}    # {}\n'.format(num, ka, th, atms)
            newLines.append(newline)

        # ADD DIHEDRAL COEFFS
        elif dcoef and re.search(pattern_coeffs, line):
            h, num, atms           = line.split()
            atm1, atm2, atm3, atm4 = atms.split('-')
            args = getDihedral(atm1, atm2, atm3, atm4)
            args = '  '.join(args)
            newline                = '{:4} {}  # {}\n'.format(num, args, atms)
            newLines.append(newline)

        # ADD IMPROPER COEFFS
        elif icoef and re.search(pattern_coeffs, line):
            h, num, atms           = line.split()
            atm1, atm2, atm3, atm4 = atms.split('-')
            v1, v2, v3             = getImproper(atm1, atm2, atm3, atm4)
            args = '  '.join([v1, v2, v3])
            newline                = '{:4} {} # {}\n'.format(num, args, atms)
            newLines.append(newline)

        # ADD PARTIAL CHARGES AND FIND BOX MIN & MAX
        elif atoms and re.search('^[0-9]{0,5} [0-9]{0,5} [0-9]{0,5}', line):
            ord, num, ID, pc, x, y, z, h, name = line.split()

            # SAVE VALUES OF X, Y, Z
            xvals.append(float(x))
            yvals.append(float(y))
            zvals.append(float(z))

            # GET NEW PARTIAL CHARGE
            pc = getAtomPartialCharge(name)
            newline = '{:4} {:3} {:3} {:>6.3f} {:>11.6f} {:>11.6f} {:>11.6f}   # {}\n'\
                .format(ord, num, ID, pc, float(x), float(y), float(z), name)
            newLines.append(newline)
            pcharges.append(pc)

        # FINISHED EDITTING LINES
        elif re.search('^\s*Bonds\s*$', line):
            atoms = False
            newLines.append(line)

        # SAVE UNCHANGED LINES
        else:
            newLines.append(line)

    # SORT X, Y, Z, VALS FOR MIN MAX VALUES
    xvals.sort()
    yvals.sort()
    zvals.sort()

    # SUB IN BOX MAX AND MIN AND REMOVE HASH ON BLANK LINES
    for i in range(len(newLines)):

        if 'xlo' in newLines[i]:
            newLines[i] = '{:6.3f} {:6.3f}  xlo xhi\n'.format(xvals[0] - 1, xvals[-1] + 1)

        elif 'ylo' in newLines[i]:
            newLines[i] = '{:6.3f} {:6.3f}  ylo yhi\n'.format(yvals[0] - 1, yvals[-1] + 1)

        elif 'zlo' in newLines[i]:
            newLines[i] = '{:6.3f} {:6.3f}  zlo zhi\n'.format(zvals[0] - 1, zvals[-1] + 1)

        elif re.search('^\s*#\s*$', newLines[i]):
            newLines[i] = '\n'

    # REMOVE EXCESS FILES
    sp.check_output('rm topo-in.xyz topo.out tempfile', shell=True)

    # add new line below first line- packmol's pack.pl needs it
    # lammps doesn't mind it
    newLines.insert(1, '\n')

    # WRITE .data FILE
    open(Name + '.data', "w+").writelines(newLines)

    # PRINT SUM OF PARTIAL CHARGES
    print(f"{File}   Charge: {sum(pcharges):5.5}")
