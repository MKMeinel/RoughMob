## File Descriptions

### `<component1>-<component2>-DAA.txt` and `<component1>-<component2>-DCG.txt`
These files contain diffusion coefficients for the specified component mixtures. 

- **DAA**: All-atom diffusion coefficients
- **DCG**: Coarse-grained diffusion coefficients

The data structure for both files is as follows:
- **First Column**: Mole fraction of component 1
- **Second Column**: Self-diffusion coefficient of component 1
- **Third Column**: Self-diffusion coefficient of component 2
- **Fourth Column**: Binary/mutual diffusion coefficient

### Error Files: `DAAError` and `DCGError`
These files represent the standard deviation between diffusion coefficients calculated from independent Cartesian coordinates. The column structure is identical to the corresponding DAA and DCG files.

### `<component1>-<component2>-alpha.txt`
This file details the acceleration factors (DCG/DAA) for the specified component mixtures.

The data is organized as follows:
- **First Column**: Mole fraction of component 1
- **Second Column**: Number averaged self-diffusion acceleration
- **Third Column**: Self-diffusion acceleration factor of component 1
- **Fourth Column**: Self-diffusion acceleration factor of component 2
- **Fifth Column**: Binary diffusion acceleration factor


