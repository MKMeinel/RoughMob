import numpy as np
import numpy.matlib
from scipy import stats

#-------------------------------------------------------------------------------------------------
# Unwrap and Scale Coordinates to Box Dimension 
#-------------------------------------------------------------------------------------------------
def f_unwrap(coords, dim, scaled=True):
    """
    Unwraps and scales coordinates relative to the box dimensions.

    Parameters:
    coords (ndarray): Coordinates to be unwrapped.
    dim (ndarray): Dimensions of the simulation box.
    scaled (bool): If True, assumes lammps output xs, ys, zs.

    Returns:
    ndarray: Unwrapped and scaled coordinates.
    """
    def adjust_dist(dists):
        if dists <= -0.5:
            return +1
        elif dists >= 0.5:
            return -1
        else:
            return 0
    adjust_dist_vect = np.vectorize(adjust_dist)

    if scaled:
        dist = coords - coords[0, :]
        returncoords = (coords + adjust_dist_vect(dist)) * dim
    else:
        dist = (coords - coords[0, :]) / dim
        returncoords = coords + adjust_dist_vect(dist) * dim
    return returncoords

#-------------------------------------------------------------------------------------------------
# Calculate Center of Mass
#-------------------------------------------------------------------------------------------------
def f_com(coords, masses, return_coords=False):
    """
    Calculates the center of mass of a molecule or system.

    Parameters:
    coords (ndarray): Coordinates of the molecule/system.
    masses (ndarray): Array containing the masses of atoms in the molecule/system.
    return_coords (bool): If True, also returns adjusted coordinates.

    Returns:
    ndarray: Center of mass, and optionally coordinates shifted to center of mass.
    """
    com = 1 / np.sum(masses) * np.dot(masses, coords)
    if return_coords:
        return com, coords - com
    else:
        return com

#-------------------------------------------------------------------------------------------------
# Calculate Evenly Distributed Spherical Gridpoints
#-------------------------------------------------------------------------------------------------
def f_sphere_gridpoints(n, r):
    """
    Calculates evenly distributed grid points on a sphere.

    Parameters:
    n (int): Number of grid points.
    r (float): Radius of the sphere.

    Returns:
    ndarray: Coordinates of grid points with COM at (0,0,0).
    """
    golden = (1 + np.sqrt(5.0)) / 2
    ind = np.arange(-n / 2, n / 2)
    phi = 2.0 * np.pi * ind / golden
    sintheta = 2 * ind / n
    costheta = np.sqrt(1 - sintheta**2)
    p_grid = np.array([costheta * np.cos(phi), costheta * np.sin(phi), sintheta])
    return (p_grid * r).T

#-------------------------------------------------------------------------------------------------
# Ray-Sphere Intersection
#-------------------------------------------------------------------------------------------------
def f_ray_sphere(C, Rad, P0, rays):
    """
    Calculates the intersection points of rays with spheres.

    Parameters:
    C (ndarray): Centers of spheres.
    Rad (ndarray): Radii of spheres.
    P0 (ndarray): Start points of the rays.
    rays (ndarray): End points of the rays (directions).

    Returns:
    ndarray: Intersection points.
    """

    s_num=C.shape[0] 
    r_num=rays.shape[0] 
    
    P0 = np.repeat(P0[np.newaxis,:],r_num,axis=0)
    P0 = np.tile(P0.T,s_num).T
    
    Rays=np.tile(rays.T,s_num).T
    
    Rdir=(Rays-P0)/np.matlib.repmat(np.linalg.norm((Rays-P0),axis=1),3,1).T
    
    C=np.repeat(C[:,:],r_num,axis=0)

    
    Rad=np.repeat(Rad,r_num,axis=0)
    Rad2 = Rad*Rad

    a = np.sum(Rdir*Rdir, axis=1)
    b = 2*np.sum(Rdir*(P0-C),axis=1)
    c = np.sum((P0-C)*(P0-C),axis=1)-Rad2
    d = b*b-4*a*c
    

    def OK(a,b,d):
        epsilon=0.000001
        if (d < epsilon):
            return 0.0, 0.0
        else:
            t1 = (-b+np.sqrt(d))/(2*a)
            t2 = (-b-np.sqrt(d))/(2*a)
            if (t1 < epsilon):
                t1 = 0.0
            if (t2 < epsilon):
                t2 = 0.0
            return t1, t2

    vOK = np.vectorize(OK)
    
    t1, t2 = vOK(a,b,d)
    P1, P2 = P0+np.matlib.repmat(t1,3,1).T*Rdir, P0+np.matlib.repmat(t2,3,1).T*Rdir
    
    
    D_P1 = np.linalg.norm(P1-P0,axis=1)
    indR=np.arange(0,r_num)
    D_P1 = np.reshape(D_P1,(s_num,r_num),order='A')
    indS_P1=np.argmax(D_P1,axis=0)
    D_P1=np.max(D_P1,axis=0)
    ind=indS_P1*r_num+indR
        
    P1=P1[ind]
    
    return P1 

    
    
def get_box_grid_MC(s, e, n):
    """
    Generates a grid of points within a box using Monte Carlo method.

    Parameters:
    s (ndarray): Start coordinates of the box.
    e (ndarray): End coordinates of the box.
    n (int): Number of points.

    Returns:
    ndarray: Grid points within the box.
    """
    grid = (e - s) * np.random.random((n, 3)) + s
    return grid

#-------------------------------------------------------------------------------------------------
# Calculate the Roughness Difference
#-------------------------------------------------------------------------------------------------
def f_calc_Ra(surf, ref, raystart):
    """
    Calculates the molecular roughness difference (arithmetical mean deviation Ra) between two surfaces (all-atom and coarse-grained).

    Parameters:
    surf (ndarray): Surface points of the object.
    ref (ndarray): Reference surface points.
    raystart (ndarray): Origin of the rays used in measurements.

    Returns:
    tuple: Arithmetical mean deviation (Ra), Skewness (Rsk), Kurtosis (Rku), and average radius (rAA).
    """
    D_surf = np.linalg.norm(surf - raystart, axis=1)
    D_ref = np.linalg.norm(ref - raystart, axis=1)
    
    D = D_surf - D_ref

    Ra = np.mean(np.abs(D))  # Arithmetical mean deviation of the assessed profile
    Rq = np.sqrt(np.mean(D**2))  # Root mean squared
    Rsk = (1 / Rq**3) * np.mean(D**3)  # Skewness
    Rku = (1 / Rq**4) * np.mean(D**4)  # Kurtosis
    rAA = np.mean(D_surf)
    
    return Ra, Rsk, Rku, rAA

