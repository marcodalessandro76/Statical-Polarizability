def get_atoms(str):
    atoms = []
    # remove numbers from str
    for c in str:
        if c.isdigit():
            str = str.replace(c,'')
    # split elements to identify the atoms
    while True:
        if len(str)==1:
            atoms.append(str)
            break
        if len(str)==2 and str[1].islower():
            atoms.append(str)
            break
        if str[1].isupper():
            atoms.append(str[:1])
            str = str[1:]
            continue
        if str[1].islower():
            atoms.append(str[:2])
            str = str[2:]
            continue
    # remove duplicates (and loose ordering)
    atoms = list(set(atoms))
    return atoms

def molecule_inlist(mol,psp_list):
    atoms = get_atoms(mol)
    inlist = True
    for a in atoms:
        if a not in psp_list:
            inlist = False
        break
    return inlist
