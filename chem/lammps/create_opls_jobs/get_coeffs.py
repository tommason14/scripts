import re, sys, itertools, os
kj2kcal = 4.1868
factor  = 2

# GROUPS FOR BONDS, ANGLES, DIHEDRALS, IMPROPERS
def getGroup(Atom, path_to_ff):

    # GET ATOM GROUP FROM ff
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search('^'+Atom+'\s', line):
                group = line.split()[1]
                return group

            # IF GET TO BOND SECTION
            elif re.search('^\s*BONDS', line):
                sys.exit('Could not find atom {}'.format(Atom))


def getAtomData(Atom, path_to_ff):

    # GET ATOM DATA FROM ff
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search('^'+Atom+'\s', line):
                line = line.split()
                pot  = float(line[5])
                pars = float(line[6])/kj2kcal
                return pars, pot

            # IF GET TO BOND SECTION
            elif re.search('^\s*BONDS', line):
                sys.exit('Could not find atom {}'.format(Atom))


def getAtomPartialCharge(Atom, path_to_ff):

    # GET ATOM DATA FROM ff
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF LINE STARTS WITH Atom NAME
            if re.search('^'+Atom+'\s', line):
                line = line.split()
                q    = float(line[3])
                return q

            # IF GET TO BOND SECTION
            elif re.search('^\s*BONDS', line):
                sys.exit('Could not find atom {}'.format(Atom))


def getBond(myAtom1, myAtom2, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF FOUND BOND SECTION
            if re.search('BONDS', line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2
            elif found and re.search('^'+Atom1+'\s*'+Atom2+'\s', line):
                line = line.split()
                Re   = float(line[3])
                kr   = float(line[4])/kj2kcal/factor
                return kr, Re

            # IF LINE STARTS WITH Atom2 Atom1
            elif found and re.search('^'+Atom2+'\s*'+Atom1+'\s', line):
                line = line.split()
                Re   = float(line[3])
                kr   = float(line[4])/kj2kcal/factor
                return kr, Re

            # IF GET TO ANGLE SECTION
            elif re.search('^\s*ANGLES', line):
                sys.exit('Could not find bond {} {}'.format(Atom1, Atom2))


def getAngle(myAtom1, myAtom2, myAtom3, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF FOUND ANGLES SECTION
            if re.search('^\s*ANGLES', line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2 Atom3
            elif found and re.search('^'+Atom1+'\s*'+Atom2+'\s*'+Atom3+'\s', line):
                line = line.split()
                th   = float(line[4])
                ka   = float(line[5])/kj2kcal/factor
                return ka, th

            # IF LINE STARTS WITH Atom3 Atom2 Atom1
            elif found and re.search('^'+Atom3+'\s*'+Atom2+'\s*'+Atom1+'\s', line):
                line = line.split()
                th   = float(line[4])
                ka   = float(line[5])/kj2kcal/factor
                return ka, th

            # IF GET TO DIHEDRALS SECTION
            elif re.search('^\s*DIHEDRALS', line):
                sys.exit('Could not find angle {} {} {}'.format(Atom1, Atom2, Atom3))


def getDihedral(myAtom1, myAtom2, myAtom3, myAtom4, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)
    Atom4 = getGroup(myAtom4, path_to_ff)

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, 'r+') as f:
        for line in f:

            # IF FOUND DIHEDRALS SECTION
            if re.search('^\s*DIHEDRALS', line):
                found = True

            # IF LINE STARTS WITH Atom1 Atom2 Atom3 Atom4
            elif found and re.search('^'+Atom1+'\s*'+Atom2+'\s*'+Atom3+'\s*'+Atom4+'\s', line):
                line = line.split()
                v1   = float(line[5])/kj2kcal
                v2   = float(line[6])/kj2kcal
                v3   = float(line[7])/kj2kcal
                v4   = float(line[8])/kj2kcal
                return v1, v2, v3, v4

            # IF LINE STARTS WITH Atom4 Atom3 Atom2 Atom1
            elif found and re.search('^'+Atom4+'\s*'+Atom3+'\s*'+Atom2+'\s*'+Atom1+'\s', line):
                line = line.split()
                v1   = float(line[5])/kj2kcal
                v2   = float(line[6])/kj2kcal
                v3   = float(line[7])/kj2kcal
                v4   = float(line[8])/kj2kcal
                return v1, v2, v3, v4

            # IF GET TO DIHEDRALS SECTION
            elif re.search('^\s*IMPROPER', line):
                sys.exit('Could not find dihedral {} {} {} {}'.format(Atom1, Atom2, Atom3, Atom4))


def getImproper(myAtom1, myAtom2, myAtom3, myAtom4, path_to_ff):

    # GET FF ATOM NAME
    Atom1 = getGroup(myAtom1, path_to_ff)
    Atom2 = getGroup(myAtom2, path_to_ff)
    Atom3 = getGroup(myAtom3, path_to_ff)
    Atom4 = getGroup(myAtom4, path_to_ff)

    perms = list(itertools.permutations([Atom1, Atom2, Atom3, Atom4], 4))

    # GET BOND DATA FROM ff
    found = False
    with open(path_to_ff, 'r+') as f:
        for line in f:
            # IF FOUND DIHEDRALS SECTION
            if re.search('^\s*IMPROPER', line):
                found = True

            elif found:
                for perm in perms:
                    a1, a2, a3, a4 = perm
                    if re.search('^'+a1+'\s*'+a2+'\s*'+a3+'\s*'+a4+'\s', line):
                        line = line.split()
                        v1   = float(line[5])/kj2kcal/factor
                        v2   = float(line[6])/kj2kcal/factor
                        v3   = float(line[7])/kj2kcal/factor
                        v4   = float(line[8])/kj2kcal/factor
                        return v1, v2, v3, v4

    sys.exit('Could not find improper {} {} {} {}'.format(Atom1, Atom2, Atom3, Atom4))
