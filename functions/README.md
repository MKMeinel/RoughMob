
# RoughMob Analysis Tools

This repository contains a set of Python tools designed for the usage of the RoughMob method with molecular dynamics simulations (carried out with LAMMPS), and coarse-grained models (generated with VOTCA). 

## Contents

- `RMgeometry.py`: Contains functions for geometric calculations, such as unwrapping coordinates, calculating the center of mass, generating spherical grid points, and ray-sphere intersection.

- `RMreadfiles.py`: Provides utility functions to read and process data from LAMMPS dump files and data files and VOTCA input files. It includes methods to extract atom types, masses, and coordinates from simulation files.

- `RManalysis.py`: A script that brings together functions from `geometry.py` and `readfiles.py` to perform the specific RoughMob tasks, such as calculating the roughness of single CG bead molecules in a trajectory and determining the outer volume of molecules.

## Getting Started

To use these tools, ensure you have Python installed along with the necessary libraries: NumPy and SciPy. Clone this repository to your local machine, and you can import these scripts into your Python environment to use their functions.

### Prerequisites

- Python 3.x
- NumPy
- SciPy

### Installation

Clone the repository using:

```bash
git clone https://github.com/MKMeinel/RoughMob
```

### Usage

Add the repository to your PYTHONPATH. 

## Contributing

Contributions to this project are welcome. Please feel free to fork the repository, make changes, and submit pull requests. For major changes, please open an issue first to discuss what you would like to change. 

## Contact

For any queries or further assistance, please reach out to us via GitHub issues or email.
