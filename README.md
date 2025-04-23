# Residue Displacement Visualization Tool

This Python script computes per-residue Cα (alpha carbon) displacement vectors between two aligned protein structures in PDB format. It generates:

- A `.bild` file for visualizing displacement vectors as arrows in ChimeraX.
- A `.cxc` ChimeraX command script to automate loading and styling of the structures and vectors.

## Features

- Calculates displacement vectors and magnitudes between corresponding Cα atoms.
- Color-codes arrows from blue (minimal displacement) to red (maximal displacement).
- Provides summary statistics: minimum, maximum, and average displacement magnitudes.
- Generates ChimeraX scripts for streamlined visualization.

## Requirements

- Python 3.6 or higher
- [Biopython](https://biopython.org/)
- [NumPy](https://numpy.org/)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) (for visualization)

## Installation

Install the required Python packages using pip:

```bash
pip install biopython numpy
```

Ensure ChimeraX is installed for visualization purposes.

## Usage

1. Prepare two aligned PDB files:
   - `model1-9fuq.pdb`: Reference structure.
   - `model2-9fur.pdb`: Structure aligned to the reference.

2. Place the script in the same directory as the PDB files.

3. Run the script:

   ```bash
   python displacement_visualization.py
   ```

4. The script will output:
   - `displacement_vectors.bild`: Arrow representations of displacements.
   - `visualize_displacements.cxc`: ChimeraX command script for visualization.

5. Open the `visualize_displacements.cxc` with ChimeraX.


## Output Description

- **displacement_vectors.bild**: Contains arrow graphics representing displacement vectors between corresponding Cα atoms.
- **visualize_displacements.cxc**: Automates the loading of PDB files, application of the displacement vectors, and sets up the visualization environment in ChimeraX.

## Notes

- Ensure that both PDB files are properly aligned before running the script. Misaligned structures will result in inaccurate displacement calculations.
- The script matches residues based on chain ID and residue number. Ensure consistency between the two PDB files.

## License

This project is licensed under the MIT License.
