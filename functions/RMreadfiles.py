import numpy as np
import xml.etree.ElementTree as ET
from itertools import islice

#-------------------------------------------------------------------
#-----Lammps Dump Trajectory----------------------------------------
#-------------------------------------------------------------------

def get_natoms_from_dump(traj):
    """
    Extracts the total number of atoms from a LAMMPS dump trajectory.

    Parameters:
    traj (str): Path to the LAMMPS dump trajectory.

    Returns:
    int: Total number of atoms in the system.
    """
    natoms = np.genfromtxt(traj, dtype=int, skip_header=3, max_rows=1)
    return natoms

def read_dump_frame(traj, natoms, framenum=1, atype=True, mol=True, vel=False):
    """
    Reads a specific frame from a LAMMPS dump trajectory.

    Parameters:
    traj (str): Path to the LAMMPS dump trajectory.
    natoms (int): Total number of atoms.
    framenum (int): Frame number to read (1-indexed).
    atype (bool): If True, reads atom types.
    mol (bool): If True, reads molecule IDs.
    vel (bool): If True, reads velocities.

    Returns:
    tuple: Depending on the parameters, a combination of molecule IDs, atom types, positions, and velocities.
    """
    if atype and mol and vel:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "a_type", "mol", "x","y","z", "vx","vy","vz"])
        mol=frame["mol"]
        a_type=frame["a_type"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)
        vel=np.stack((frame["vx"],frame["vy"],frame["vz"]), axis=1)
        
        return mol, a_type, pos, vel
        
    elif atype and mol:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "a_type", "mol", "x","y","z"])
        mol=frame["mol"]
        a_type=frame["a_type"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)
        return mol, a_type, pos



    elif atype and vel:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "a_type", "x","y","z", "vx","vy","vz"])
        a_type=frame["a_type"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)
        vel=np.stack((frame["vx"],frame["vy"],frame["vz"]), axis=1)
        return mol, a_type, pos, vel

        
    elif mol and vel:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "mol", "x","y","z", "vx","vy","vz"])
        mol=frame["mol"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)
        vel=np.stack((frame["vx"],frame["vy"],frame["vz"]), axis=1)
        return mol, pos, vel

        
    
    elif mol:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "mol", "x","y","z"])
        mol=frame["mol"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)

        return mol, pos
    
    elif atype:
        frame = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+9, max_rows=natoms, names=["id", "a_type", "x","y","z"])
        a_type=frame["a_type"]
        pos=np.stack((frame["x"],frame["y"],frame["z"]), axis=1)
        return a_type, pos


def get_box_from_dump(traj, natoms, framenum=1):
    """
    Extracts the simulation box coordinates and dimensions from a LAMMPS dump trajectory.

    Parameters:
    traj (str): Path to the LAMMPS dump trajectory.
    natoms (int): Total number of atoms.
    framenum (int): Frame number to read (1-indexed).

    Returns:
    tuple: Box coordinates and box dimensions.
    """
    
    box = np.genfromtxt(traj, skip_header=(framenum-1)*(natoms+9)+5, max_rows=3)
    box_coord=box.T
    box_dim=np.diff(box_coord,axis=0)
    return box_coord, box_dim

def get_coordinate_style_from_dump(traj):
    """
    Determines the style of coordinates used in a LAMMPS dump trajectory.

    Parameters:
    traj (str): Path to the LAMMPS dump trajectory.

    Returns:
    list: Style of coordinates.
    """
    
    with open(traj) as f:
        line=f.read().split('\n')[8]
    return line.split()[-3:]

#--------------------------------------------------------------------------
#-------LAMMPS Data file---------------------------------------------------
#--------------------------------------------------------------------------

def get_numbers_from_data(datafile):
    """
    Extracts numbers of atoms, bonds, angles, dihedrals, and impropers from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple (integers): Numbers of atoms, bonds, angles, dihedrals, and impropers.
    """
    natoms=0
    nbonds=0
    nangles=0
    ndihedrals=0
    nimpropers=0
    with open(datafile) as myFile:
        for num, line in enumerate(myFile,1):
            if 'atoms' in line:
                natoms=int(line.strip().split()[0])
            if 'bonds' in line:
                nbonds=int(line.strip().split()[0])
            if 'angles' in line:
                nangles=int(line.strip().split()[0])
            if 'dihedrals' in line:
                ndihedrals=int(line.strip().split()[0])
            if 'impropers' in line:
                nimpropers=int(line.strip().split()[0])
    return natoms, nbonds, nangles, ndihedrals, nimpropers

def get_types_from_data(datafile):
    """
    Extracts types of atoms, bonds, angles, dihedrals, and impropers from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Types of atoms, bonds, angles, dihedrals, and impropers.
    """
    atomtypes=0
    bondtypes=0
    angletypes=0
    dihedraltypes=0
    impropertypes=0
    with open(datafile) as myFile:
        for num, line in enumerate(myFile,1):
            if 'atom types' in line:
                atomtypes=int(line.strip().split()[0])
            if 'bond types' in line:
                bondtypes=int(line.strip().split()[0])
            if 'angle types' in line:
                angletypes=int(line.strip().split()[0])
            if 'dihedral types' in line:
                dihedraltypes=int(line.strip().split()[0])
            if 'improper types' in line:
                impropertypes=int(line.strip().split()[0])
    return atomtypes, bondtypes, angletypes, dihedraltypes, impropertypes 

def get_radiitypes_from_data(datafile):
    """
    Extracts radii types from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    ndarray: Radii types.
    """
    
    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Pair Coeffs' in line:
                lcoeffs=int(num)
    atomtypes, _, _, _,_=get_types_from_data(datafile)
    
    radii=np.genfromtxt(datafile,skip_header=(lcoeffs+1),usecols=(2), max_rows=atomtypes)

    return radii/2




def get_radii_from_data(datafile, a_type_m):
    """
    Retrieves atomic radii from a LAMMPS data file based on atom types.

    Parameters:
    datafile (str): Path to the LAMMPS data file.
    a_type_m (array): Array of atom types.

    Returns:
    ndarray: Array of radii corresponding to the atom types.
    """
    radii = get_radiitypes_from_data(datafile)
    radii_m = np.array(radii[np.array(a_type_m).astype(int) - 1]).flatten()
    return radii_m


def get_masstypes_from_data(datafile):
    """
    Extracts mass types from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    ndarray: Array of atomic masses.
    """
    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Masses' in line:
                lmasses = int(num)
    atomtypes, _, _, _, _ = get_types_from_data(datafile)
    masses = np.genfromtxt(datafile, skip_header=(lmasses + 1), usecols=(1), max_rows=atomtypes)
    return masses

def get_masses_from_data(datafile, a_type_m):
    """
    Retrieves masses from a LAMMPS data file based on atom types.

    Parameters:
    datafile (str): Path to the LAMMPS data file.
    a_type_m (array): Array of atom types.

    Returns:
    ndarray: Array of masses corresponding to the atom types.
    """
    mass = get_masstypes_from_data(datafile)
    masses = np.array(mass[np.array(a_type_m).astype(int) - 1]).flatten()
    return masses

def get_atoms_from_data(datafile):
    """
    Extracts atom information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Arrays of molecule IDs, atom types, and atom charges.
    """
    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Atoms' in line:
                latoms = int(num)
    natoms, _, _, _, _ = get_numbers_from_data(datafile)
    atoms = np.genfromtxt(datafile, skip_header=(latoms + 1), usecols=(0, 1, 2, 3), max_rows=natoms)
    atoms = atoms[atoms[:, 0].argsort()]
    return atoms[:, 1], atoms[:, 2], atoms[:, 3]

def get_bonds_from_data(datafile):
    """
    Extracts bond information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Arrays of bond types and connected atoms.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Bonds' in line:
                lbonds=int(num)
    _, nbonds, _, _,_=get_numbers_from_data(datafile)

    bonds=np.genfromtxt(datafile,skip_header=(lbonds+1),usecols=(0,1,2,3), max_rows=nbonds)
    bonds=bonds[bonds[:,0].argsort()]
    return bonds[:,1], bonds[:,2:4]

    
def get_angles_from_data(datafile):
    """
    Extracts angle information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Arrays of angle type and connected atoms.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Angles' in line:
                langles=int(num)
    _, _, nangles, _,_=get_numbers_from_data(datafile)

    angles=np.genfromtxt(datafile,skip_header=(langles+1),usecols=(0,1,2,3,4), max_rows=nangles)
    angles=angles[angles[:,0].argsort()]

    #angles[:,0] angle id 
    #angles[:,1] angle type
    #angles[:,2:5] 
    
    return angles[:,1], angles[:,2:5]

    

def get_dihedrals_from_data(datafile):
    """
    Extracts dihedral information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Arrays of dihedral types and connected atoms.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Dihedrals' in line:
                ldihedrals=int(num)
    _, _, _, ndihedrals,_=get_numbers_from_data(datafile)

    dihedrals=np.genfromtxt(datafile,skip_header=(ldihedrals+1),usecols=(0,1,2,3,4,5), max_rows=ndihedrals)
    dihedrals=dihedrals[dihedrals[:,0].argsort()]

    #dihedrals[:,0] dihedral id 
    #dihedrals[:,1] dihedral type
    #dihedrals[:,2:6] 
    
    return dihedrals[:,1], dihedrals[:,2:6]

def get_impropers_from_data(datafile):
    """
    Extracts improper dihedral information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Arrays of improper dihedral types and connected atoms.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Impropers' in line:
                limpropers=int(num)
    _, _, _, _,nimpropers=get_numbers_from_data(datafile)

    impropers=np.genfromtxt(datafile,skip_header=(limpropers+1),usecols=(0,1,2,3,4,5), max_rows=nimpropers)
    impropers=impropers[impropers[:,0].argsort()]
    
    return impropers[:,1], impropers[:,2:6]


def get_box_from_data(datafile):
    """
    Extracts simulation box information from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    tuple: Coordinates of the simulation box and its dimensions.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'xlo' in line:
                lxlo=int(num)
            if 'ylo' in line:
                lylo=int(num)
            if 'zlo' in line:
                lzlo=int(num)

  
    xs,xe=np.genfromtxt(datafile,skip_header=(lxlo-1),usecols=(0,1), max_rows=1)
    ys,ye=np.genfromtxt(datafile,skip_header=(lylo-1),usecols=(0,1), max_rows=1)
    zs,ze=np.genfromtxt(datafile,skip_header=(lzlo-1),usecols=(0,1), max_rows=1)

    box_coord=np.array([(xs,ys,zs),(xe,ye,ze)])
    box_dim=np.diff(box_coord,axis=0)
    
    return box_coord, box_dim


def get_coordinates_from_data(datafile):
    """
    Extracts coordinates of atoms from a LAMMPS data file.

    Parameters:
    datafile (str): Path to the LAMMPS data file.

    Returns:
    ndarray: Array of atomic coordinates.
    """

    with open(datafile) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'Atoms' in line:
                latoms=int(num)
    natoms, _, _, _,_=get_numbers_from_data(datafile)

    atoms=np.genfromtxt(datafile,skip_header=(latoms+1),usecols=(0,4,5,6), max_rows=natoms)
    atoms=atoms[atoms[:,0].argsort()]
    
    return atoms[:,1:]





#--------------------------------------------------------------------------
#-------VOTCA Mapping File-------------------------------------------------
#--------------------------------------------------------------------------

def get_mapinfo(mapfile, AAtopo, CGtopo):
    """
    Extracts mapping information from a VOTCA mapping file.

    Parameters:
    mapfile (str): Path to the VOTCA mapping file.
    AAtopo (str): Path to the VOTCA topology file for all-atom representation.
    CGtopo (str): Path to the VOTCA topology file for coarse-grained representation.

    Returns:
    ndarray: Mapping information.
    """

    AAnbeads, AAnames, AAmasses, AAtypes = get_topoinfo(AAtopo)
    CGnbeads, CGnames, CGmasses, CGtypes = get_topoinfo(CGtopo)

    maps=ET.parse(mapfile)
    maps=maps.getroot()
    mapping=np.zeros([AAnbeads,CGnbeads])
    kk=-1
    for child in maps.iter('cg_bead'):
        kk+=1
        lists=child.findtext('beads').split()
        checklist='n'.join(lists)+'n'
        ll=-1
        for bead in AAnames:
            ll+=1
            if (bead+'n' in checklist):
                mapping[ll,kk]=int(1)

    return mapping.astype(int)

    
#--------------------------------------------------------------------------
#-------VOTCA Topology File------------------------------------------------
#--------------------------------------------------------------------------

def get_topoinfo(topofile):
    """
    Extracts topology information from a VOTCA topology file.

    Parameters:
    topofile (str): Path to the VOTCA topology file.

    Returns:
    tuple: Number of beads/atoms, names,masses and types of atoms.
    """

    topo=ET.parse(topofile)
    topo=topo.getroot()
    beads_n=topo.find('./molecules/molecule').get('nbeads')
    beads_n=int(beads_n)
    names_m=[None]*beads_n
    masses_m=np.zeros([beads_n])
    types_m=np.zeros([beads_n])
            
    ii=-1
    for beads in topo.iter('bead'):
        ii+=1
        names_m[ii]=beads.get('name')
        masses_m[ii]=beads.get('mass')
        types_m[ii]=beads.get('type')
    types_m=types_m.astype(int)

    return beads_n, names_m, masses_m, types_m

    
#------------------------------------------------------------------------
#-----Log-File (LAMMPS)---------------------------------------------------
#------------------------------------------------------------------------

def get_log_property(thermo, Property):
    """
    Extracts a specific property from a LAMMPS log file.

    Parameters:
    thermo (str): Path to the LAMMPS log file.
    Property (str): Name of the property to extract.

    Returns:
    ndarray: Extracted property values.
    """
    data=np.genfromtxt(thermo,names=True)
    if Property in data.dtype.names:
        prop=data[Property]
        return prop

    else:
        print("Error: Property not found.")
        print("Valid properties are:")
        for name in data.dtype.names:
            print(name)



#------------------------------------------------------------------------
#-----tabulated Pot-File (LAMMPS)----------------------------------------
#------------------------------------------------------------------------

def get_CG_radius(potfile):
    """
    Determines the coarse-grained radius from a LAMMPS potential file.

    Parameters:
    potfile (str): Path to the LAMMPS potential file.

    Returns:
    float: Coarse-grained radius.
    """

    data=np.genfromtxt(potfile,skip_header=1, usecols=(1,2))
    r=data[:,0]
    pot=data[:,1]

    ii=-1
    for p in pot:
        ii+=1
        if p < 0:
            return r[ii]/2


    print("Error: Radius could not be determined!")


