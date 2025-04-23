from Bio.PDB import PDBParser, is_aa
import numpy as np

pdb1 = "model1-9fuq.pdb" #Provide model 1 PDB file
pdb2 = "model2-9fur.pdb" #Provide model 2 PDB file aligned to model 1

def get_ca_coords_by_residue(structure):
    ca_dict = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if is_aa(res, standard=True) and res.id[0] == ' ' and 'CA' in res:
                    res_id = (chain.id, res.id[1])  # (chain ID, residue number)
                    ca_dict[res_id] = res['CA'].coord
    return ca_dict

parser = PDBParser(QUIET=True)
structure1 = parser.get_structure("s1", pdb1)
structure2 = parser.get_structure("s2", pdb2)

coords1_dict = get_ca_coords_by_residue(structure1)
coords2_dict = get_ca_coords_by_residue(structure2)

# Find common residues
common_residues = set(coords1_dict.keys()) & set(coords2_dict.keys())

if not common_residues:
    raise ValueError("No matching residues found between the two structures.")

# Sort residues for consistent ordering
common_residues = sorted(common_residues)

coords1 = np.array([coords1_dict[res_id] for res_id in common_residues])
coords2 = np.array([coords2_dict[res_id] for res_id in common_residues])

displacements = coords2 - coords1
magnitudes = np.linalg.norm(displacements, axis=1)

min_mag = np.min(magnitudes)
max_mag = np.max(magnitudes)
avg_mag = np.mean(magnitudes)

print(f"Minimum displacement magnitude: {min_mag:.2f} Å")
print(f"Maximum displacement magnitude: {max_mag:.2f} Å")
print(f"Average displacement magnitude: {avg_mag:.2f} Å")

def color_from_magnitude(mag, max_mag):
    scale = min(1.0, mag / max_mag)
    r = scale
    g = 0.0
    b = 1.0 - scale
    return r, g, b

lines = []

for start, vec, mag in zip(coords1, displacements, magnitudes):
    end = start + vec
    shaft_length = 0.8 * mag
    shaft_end = start + (vec / mag) * shaft_length
    r, g, b = color_from_magnitude(mag, max_mag)

    # Set color
    lines.append(f".color {r:.3f} {g:.3f} {b:.3f}")

    # Cylinder (shaft)
    lines.append(f".cylinder {start[0]:.3f} {start[1]:.3f} {start[2]:.3f} "
                 f"{shaft_end[0]:.3f} {shaft_end[1]:.3f} {shaft_end[2]:.3f} 0.25")

    # Cone (arrowhead)
    lines.append(f".cone {shaft_end[0]:.3f} {shaft_end[1]:.3f} {shaft_end[2]:.3f} "
                 f"{end[0]:.3f} {end[1]:.3f} {end[2]:.3f} 0.5")

with open("displacement_vectors.bild", "w") as f:
    f.write("\n".join(lines))

# Generate ChimeraX command script
with open("visualize_displacements.cxc", "w") as f:
    f.write(f"open {pdb1}\n")
    f.write(f"open {pdb2}\n")
    f.write("matchmaker #1 to #2\n")
    f.write("open displacement_vectors.bild\n")
    f.write(f"key blue:{min_mag:.2f} red:{max_mag:.2f} size 0.018,0.27 pos 0.78,0.55 numericLabelSpacing proportional ticks true tickLength 10 tickThickness 2\n")
    f.write(f'2dlabels text "Displacement Magnitude (Å)" xpos 0.2 ypos 0.9 size 24\n')
    f.write(f'2dlabels text "Min: {min_mag:.2f} Å" xpos 0.7 ypos 0.46 size 20\n')
    f.write(f'2dlabels text "Max: {max_mag:.2f} Å" xpos 0.7 ypos 0.43 size 20\n')
    f.write(f'2dlabels text "Avg: {avg_mag:.2f} Å" xpos 0.7 ypos 0.40 size 20\n')
    f.write(f"hide #!1 models\n")
    f.write(f'windowsize 600 600\n')
    f.write(f'set bgColor white\n')
    f.write(f'hide atoms\n')
    f.write(f'preset cartoons/nucleotides licorice/ovals\n')
    f.write(f'color #2 darkgrey\n')
    f.write(f'view\n')
print("displacement_vectors.bild and visualize_displacements.cxc generated.")
