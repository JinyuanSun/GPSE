#!/usr/bin/env python

import subprocess
import numpy as np
from Bio.PDB import PDBIO, PDBParser, Atom, Residue, Structure, Chain, Model


def translate_and_rotate_structure(structure, t, u):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.coord
                    coord = np.dot(u, coord) + t
                    atom.coord = coord
    return structure


def parse_matrix(matrix_file):
    t = []
    u = []
    x = open(matrix_file, "r").readlines()
    for y in x[2:5]:
        t_i, u_x, u_y, u_z = [float(z) for z in y.strip().split()[1:]]
        t.append(t_i)
        u.append([u_x, u_y, u_z])
    return np.array(t), np.array(u)



def copy_lig_residue_to_structure(
    from_struct, to_struct, from_chain_id, from_residue_id, to_chain_id, to_residue_id=1
):
    # get the first model in the structure
    found_lig = 0
    from_chain = from_struct[0][from_chain_id]
    try:
        to_chain = to_struct[0][to_chain_id]
    except:
        to_chain = Chain.Chain(to_chain_id)
        to_struct[0].add(to_chain)
    for residue in from_chain:
        if residue.id[1] == from_residue_id:
            residue_to_copy = residue
            found_lig = 1
    if found_lig == 0:
        raise Exception("Sorry, Ligand not Found!")

    # create a residue object
    residue = Residue.Residue(
        (residue_to_copy.id[0], to_residue_id, residue_to_copy.id[2]),
        residue_to_copy.get_resname(),
        residue_to_copy.get_segid(),
    )
    # create atom objects
    for atom in residue_to_copy.get_atoms():
        atom = Atom.Atom(
            atom.name,
            atom.coord,
            atom.bfactor,
            atom.occupancy,
            atom.altloc,
            atom.fullname,
            atom.serial_number,
            element=atom.element,
        )
        residue.add(atom)

    # add the residue to the chain
    to_chain.add(residue)
    return to_struct, residue


def superimpose_structures(structure1, structure2):
    """
    structure2 is the reference (fixed)
    """
    outname = f"aligned_{structure1}"
    cmd = f"TMalign {structure1} {structure2} -m {structure1}_{structure2}_matrix.txt"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    output = result.stdout.decode()
    lines = output.split("\n")
    for line in lines:
        if "RMSD" in line:
            rmsd = float(line.split(",")[1].split("=")[1])
    t, u = parse_matrix(f"{structure1}_{structure2}_matrix.txt")
    # print(t, u)
    parser = PDBParser(QUIET=True)

    structure1 = parser.get_structure("1", structure1)
    structure2 = parser.get_structure("2", structure2)
    structure1 = translate_and_rotate_structure(structure1, t, u)

    return rmsd, structure1


class LigAligner:
    def __init__(self, s1_name, s2_name) -> None:
        self.parser = PDBParser(QUIET=True)
        self.io = PDBIO()
        self.s1_name = s1_name
        self.s2_name = s2_name

    def superimpose_structures(self):
        """
        structure2 is the reference (fixed)
        """
        self.outname = f"aligned_{self.s1_name}"
        cmd = f"TMalign {self.s1_name} {self.s2_name} -m {self.s1_name}_{self.s2_name}_matrix.txt"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        output = result.stdout.decode()
        lines = output.split("\n")
        for line in lines:
            if "RMSD" in line:
                rmsd = float(line.split(",")[1].split("=")[1])
        t, u = parse_matrix(f"{self.s1_name}_{self.s2_name}_matrix.txt")
        parser = PDBParser(QUIET=True)
        self.structure1 = parser.get_structure("1", self.s1_name)
        self.structure2 = parser.get_structure("2", self.s2_name)
        self.structure1 = translate_and_rotate_structure(self.structure1, t, u)
        return rmsd

    def copy_lig(self, from_chain_id, from_residue_id, to_chain_id, to_residue_id=1):
        self.structure1, self.new_residue = copy_lig_residue_to_structure(
            from_struct=self.structure2,
            from_chain_id=from_chain_id,
            from_residue_id=from_residue_id,
            to_struct=self.structure1,
            to_chain_id=to_chain_id,
            to_residue_id=to_residue_id,
        )

    def write_ligand(self, chain_id="B"):
        lig_structure = Structure.Structure("Lig")
        lig_model = Model.Model(0)
        lig_chain = Chain.Chain(chain_id)
        lig_chain.add(self.new_residue)
        lig_model.add(lig_chain)
        lig_structure.add(lig_model)
        self.io.set_structure(lig_structure)
        self.io.save("ligand_from_template.pdb")

    def write_out(self, name=None):
        self.io.set_structure(self.structure1)
        if name == None:
            self.io.save(file=self.outname)
        else:
            if name.endswith(".pdb"):
                self.io.save(file=name)
            else:
                self.io.save(file=f"{name}.pdb")


def get_center_of_coordinates(residue):
    """
    Get the center of coordinates of a residue.

    Args:
    residue (Residue): a Bio.PDB Residue object

    Returns:
    tuple: the center of coordinates of the residue as (x, y, z)
    """
    x_coords = [atom.coord[0] for atom in residue]
    y_coords = [atom.coord[1] for atom in residue]
    z_coords = [atom.coord[2] for atom in residue]

    center_x = sum(x_coords) / len(x_coords)
    center_y = sum(y_coords) / len(y_coords)
    center_z = sum(z_coords) / len(z_coords)

    size_x = max(x_coords) - min(x_coords) + 2
    size_y = max(y_coords) - min(y_coords) + 2
    size_z = max(z_coords) - min(z_coords) + 2

    return (center_x, center_y, center_z), (size_x, size_y, size_z)


def run_autodock_vina(
    ligand_file,
    receptor_file,
    output_file,
    center,
    size,
    num_modes=9,
    exhaustiveness=16,
    localonly=True,
):
    # Define the command to run Autodock Vina
    cmd = [
        "vina",
        "--ligand",
        ligand_file,
        "--receptor",
        receptor_file,
        "--out",
        output_file,
        "--exhaustiveness",
        str(exhaustiveness),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
    ]

    # Add the localonly option if requested
    if localonly:
        cmd.append("--local_only")
    # print(" ".join(cmd))

    # Run the command
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    for line in result.stdout.decode().split("\n"):
        if line.startswith("Affinity"):
            affinity = float(line.split(" ")[1])
    return affinity

def prepare_ligand(lig_file):
    cmd = ["prepare_ligand", "-l", lig_file]
    res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    return lig_file.replace(".pdb", ".pdbqt")


def prepare_protein(prot_file):
    cmd = ["prepare_receptor", "-r", prot_file]
    res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    return prot_file.replace(".pdb", ".pdbqt")

def read_fasta_headers(fasta_file):
    headers = []
    with open(fasta_file, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                headers.append(line.replace("\n", ""))
    return headers

if __name__ == "__main__":
    pass
