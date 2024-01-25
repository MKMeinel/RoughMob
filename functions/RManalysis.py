import numpy as np
from scipy.optimize import least_squares

import RMgeometry as geo
import RMreadfiles as rf

def single_bead_roughness(traj, AAdata, CGpotfile, mols_n=500, skipmols_n=0, frames_n=1, scaled_coords=True, use_weights=False, weights=False):
    """
    Calculates the moleuclar roughness difference of molecules that are coarse-grained into a single bead in a trajectory.

    Parameters:
    traj (str): Path to the AA trajectory file.
    AAdata (str): Path to the LAMMPS AA data file.
    CGpotfile (str): Path to the tabulated CG potential file.
    mols_n (int): Number of molecules to process.
    skipmols_n (int): Number of molecules to skip at the start.
    frames_n (int): Number of frames to process.
    scaled_coords (bool): Whether the coordinates in the trajectory are scaled.
    use_weights (bool): If True, uses custom weights.
    weights (array): Custom weights to be used if use_weights is True.

    Returns:
    tuple: moleuclar roughness difference (Ra) and average AA radius (rAA) arrays.
    """
    
    natoms=rf.get_natoms_from_dump(traj) #total number of atoms in trajectory/system
    mol_a, a_type_a, _=rf.read_dump_frame(traj,natoms,1)#,atype=False)

    a_type_m=a_type_a[mol_a==1]
    if use_weights:
        masses_m=weights
    else:
        masses_m=rf.get_masses_from_data(AAdata,a_type_m)
    
    natoms_m=len(a_type_m)
    CG_radius=rf.get_CG_radius(CGpotfile)
    AA_radii_m=rf.get_radii_from_data(AAdata,a_type_m)

    Ra=np.zeros([mols_n, frames_n])
    rAA=np.zeros([mols_n, frames_n])

    for ff in range(frames_n): 

        mol_a, a_type_a, Coord_a=rf.read_dump_frame(traj,natoms,ff+1)#,atype=False)
        _, box_dim=rf.get_box_from_dump(traj, natoms, ff+1)
        
        for ii in range(mols_n):
            jj=int(ii+skipmols_n)
            coord_m=Coord_a[jj*natoms_m:(jj+1)*natoms_m,:]
            coord_m=geo.f_unwrap(coord_m,box_dim,scaled=scaled_coords)
            com, comcoord_m=geo.f_com(coord_m,masses_m,return_coords=True)
            surf_CG=geo.f_sphere_gridpoints(2500,CG_radius)
            surf_AA=geo.f_ray_sphere(comcoord_m,AA_radii_m,np.zeros(3),surf_CG)
            Ra[ii,ff], _, _, rAA[ii,ff]=geo.f_calc_Ra(surf_AA, surf_CG, np.zeros(3))

    return Ra, rAA

def calc_Vout_AA(traj, data, Rout, mols_n=1000, scaled_coords=True, use_weights=False, weights=False, traj_atype=True):
    """
    Calculates the outer volume of molecules in an AA trajectory.

    Parameters:
    traj (str): Path to the AA trajectory file.
    data (str): Path to the LAMMPS AA data file.
    Rout (float): Outer radius for shell volume calculation.
    mols_n (int): Number of molecules to process.
    scaled_coords (bool): Whether the coordinates in the trajectory are scaled.
    use_weights (bool): If True, uses custom weights.
    weights (array): Custom weights to be used if use_weights is True.
    traj_atype (bool): Whether to read atom types from the trajectory.

    Returns:
    tuple: Outer volume and molecular volume.
    """
    
    def f(a):
        if a >= 0:
            return 1
        else:
            return 0
    f = np.vectorize(f)


    natoms=rf.get_natoms_from_dump(traj)
    mol_a, a_type_a,Coord_a=rf.read_dump_frame(traj,natoms,1,atype=traj_atype)
    a_type_m=a_type_a[mol_a==1]
        
    if use_weights:
        masses_m=weights
    else:
        masses_m=rf.get_masses_from_data(data,a_type_m)
    natoms_m=len(a_type_m) #number of atoms in one molecule

    box, box_dim=rf.get_box_from_dump(traj, natoms, 1)

    boxs=8
    boxe=box_dim[0][0]-8
    COM=np.zeros([mols_n, 3])


    for ii in range(mols_n):
        jj=ii
        coord_m=Coord_a[jj*natoms_m:(jj+1)*natoms_m,:]
        coord_m=geo.f_unwrap(coord_m,box_dim,scaled=scaled_coords)
        COM[ii,:]=geo.f_com(coord_m,masses_m,return_coords=False)

    Vmol=np.array(box_dim[0][0]**3/mols_n)
    epsilon=2
    ii=1
    Vout_old=0
    while epsilon > 0.001 or ii < 10:

        grid=geo.get_box_grid_MC(boxs, boxe, 5000)
        Dist=np.ones([grid.shape[0],COM.shape[0]])
        
        for row in range(COM.shape[0]):
            Dist[:,row]=np.linalg.norm(grid-COM[row,:],axis=1)-Rout

        access=f(Dist)
        access=access.all(axis=1)#returns an array with True for outside of molecule and false for inside
    
        inM=grid.shape[0]-np.sum(access,axis=0) # Number of points inside of the molecule
        outM=grid.shape[0]-inM #Number of points outside of the molecule
        Vout_single=outM/grid.shape[0]*100 #outer Volume calculated in this single run
       
        Vout_new=(ii-1)/ii*Vout_old+1/ii*Vout_single #new average outer Volume 
            
        epsilon=np.abs(Vout_old-Vout_new)
        Vout_old=Vout_new
        ii+=1

        
    return Vmol*Vout_new/100, Vmol 

def calc_Vout_CG(traj, Rout, mols_n=1000, scaled_coords=True, use_weights=False, weights=False, traj_atype=True, frame=1):
    """
    Calculates the outer volume of molecules in a CG trajectory.

    Parameters:
    traj (str): Path to the CG trajectory file.
    Rout (float): Outer radius for shell volume calculation.
    mols_n (int): Number of molecules to process.
    scaled_coords (bool): Whether the coordinates in the trajectory are scaled.
    use_weights (bool): If True, uses custom weights.
    weights (array): Custom weights to be used if use_weights is True.
    traj_atype (bool): Whether to read atom types from the trajectory.
    frame (int): Frame number to read.

    Returns:
    tuple: Outer volume and molecular volume.
    """
    
    def f(a):
        if a >= 0:
            return 1
        else:
            return 0
    f = np.vectorize(f)


    natoms=rf.get_natoms_from_dump(traj)
    mol_a, a_type_a,coords=rf.read_dump_frame(traj,natoms,frame,atype=traj_atype)
    a_type_m=a_type_a[mol_a==1]
     
    box, box_dim=rf.get_box_from_dump(traj, natoms, frame)
    coords=geo.f_unwrap(coords, box_dim, scaled=scaled_coords)


    def wrap(x, boxdim):
        if x > boxdim:
            return x-boxdim
        elif x < 0:
            return boxdim+x
        else:
            return x
        
    wrap=np.vectorize(wrap)

    coords=wrap(coords, box_dim[0][0])
    
    
    boxs=8
    boxe=box_dim[0][0]-8

    Vmol=np.array(box_dim[0][0]**3/mols_n)
    epsilon=2
    ii=1
    Vout_old=0
    while epsilon > 0.001 or ii < 10:

        grid=geo.get_box_grid_MC(boxs, boxe, 5000)
        Dist=np.ones([grid.shape[0],coords.shape[0]])
        
        for row in range(coords.shape[0]):
            Dist[:,row]=np.linalg.norm(grid-coords[row,:],axis=1)-Rout

        access=f(Dist)
        access=access.all(axis=1)#returns an array with True for outside of molecule and false for inside
    
        inM=grid.shape[0]-np.sum(access,axis=0) # Number of points inside of the molecule
        outM=grid.shape[0]-inM #Number of points outside of the molecule
        Vout_single=outM/grid.shape[0]*100 #outer Volume calculated in this single run
        
        Vout_new=(ii-1)/ii*Vout_old+1/ii*Vout_single #new average outer Volume 
            
        epsilon=np.abs(Vout_old-Vout_new)
        Vout_old=Vout_new
        ii+=1

        
    return Vmol*Vout_new/100, Vmol 

def get_mix_Vout_CG(traj, nbeads1, nbeads2, Rout1, Rout2, frame=1):
    """
    Calculates the mixed outer volume of molecules in a binary mixture from mapped (with VOTCA) CG trajectory.

    Parameters:
    traj (str): Path to the CG trajectory file.
    nbeads1 (int): Number of beads in component 1.
    nbeads2 (int): Number of beads in component 2.
    Rout1 (float): Outer radius for component 1.
    Rout2 (float): Outer radius for component 2.
    frame (int): Frame number to read.

    Returns:
    tuple: Outer volume and molecular volume.
    """
 
    def f(a):
        if a >=0:
            return 1
        else:
            return 0
    f = np.vectorize(f)

    atoms_n=rf.get_natoms_from_dump(traj)
    _, coords=rf.read_dump_frame(traj, atoms_n,frame, mol=False)
    box, box_dim=rf.get_box_from_dump(traj, atoms_n, frame)
    coords=geo.f_unwrap(coords, box_dim, scaled=False)

    def wrap(x, boxdim):
        if x > boxdim:
            return x-boxdim
        elif x < 0:
            return boxdim+x
        else:
            return x
        
    wrap=np.vectorize(wrap)
    coords=wrap(coords, box_dim[0][0])
        
    ngrid=5000
    ii=1
    Vout_old=0
    epsilon=20

    Rout1=np.repeat(Rout1,nbeads1)
    Rout2=np.repeat(Rout2,nbeads2)

    Rout=np.hstack((Rout1, Rout2))
    Vmol=box_dim[0][0]**3/(nbeads1+nbeads2)

    while epsilon > 0.005 or ii < 10:
        grid=get_box_grid_MC(8, box_dim[0][0]-8, ngrid)
        Dist=np.ones([ngrid,coords.shape[0]])
        
        for row in range(coords.shape[0]):
            Dist[:,row]=np.linalg.norm(grid-coords[row,:], axis=1)-Rout[row]
        access=f(Dist)
        access=access.all(axis=1)#returns an array with True for outside of molecule and false for inside

        outM=np.sum(access) # Number of points outside molecule
        inM=grid.shape[0]-np.sum(access) # Number of points inside of molecule

        Vout_single=outM/grid.shape[0]*100
        Vout_new=(ii-1)/ii*Vout_old+1/ii*Vout_single
        epsilon=np.abs(Vout_old-Vout_new)
        
        Vout_old=Vout_new
        ii+=1
        
    return Vmol*Vout_old/100, Vout_old





