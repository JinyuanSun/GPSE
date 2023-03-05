#!/usr/bin/env python

import argparse
import glob, os
import pandas as pd
from utils import *
from joblib import Parallel, delayed
from modeller_mutation import *

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="Path to the fasta file")
    parser.add_argument("--active_chains", help="chains to build mutations")
    parser.add_argument(
        "--structure_dir", help="Path to the directory containing the structures"
    )
    parser.add_argument("--template_pdb", help="Path to the template PDB file")
    parser.add_argument("--ligand_chain", help="The chain ID of the target ligand")
    parser.add_argument("--ligand_resid", help="The residue ID of the target ligand", type=int)
    parser.add_argument(
        "--candidate_ligand_chain",
        help="The chain ID to assign of the ligand superimposed to the candidate structure",
        default="B",
    )
    parser.add_argument(
        "--candidate_ligand_resid",
        help="The residue ID to assign of the ligand superimposed to the candidate structure",
        default=1,
        type=int
    )

    args = parser.parse_args()
    return args


class Screener:
    def __init__(self, structure_dir, template_pdb, ligand_chain, ligand_resid, candidate_ligand_chain="X", candidate_ligand_resid=1) -> None:
        self.candidates = glob.glob(os.path.join(structure_dir, "*.pdb"))
        self.reference_structure = template_pdb
        self.lig_chain = ligand_chain
        self.lig_resid = ligand_resid

        self.candidate_ligand_chain = candidate_ligand_chain
        self.candidate_ligand_resid = candidate_ligand_resid

    def screen(self, id, candidate_structure):
        current_path = os.getcwd()
        id = str(id)
        try:
            os.mkdir('screening')
        except:
            pass
        child_dir = os.path.join("screening", id)
        os.makedirs(child_dir, exist_ok=True)
        subprocess.run(["cp", candidate_structure, child_dir])
        subprocess.run(["cp", self.reference_structure, child_dir])
        os.chdir(child_dir)
        candidate_structure = candidate_structure.split("/")[-1]
        ligaligner = LigAligner(candidate_structure, self.reference_structure)
        rmsd = ligaligner.superimpose_structures()
        ligaligner.write_out(f"apo-{candidate_structure}")
        rec_file = prepare_protein(f"apo-{candidate_structure}")
        ligaligner.copy_lig(self.lig_chain, self.lig_resid, self.candidate_ligand_chain, self.candidate_ligand_resid)
        ligaligner.write_ligand(self.candidate_ligand_chain)
        lig_file = prepare_ligand("ligand_from_template.pdb")
        ligaligner.write_out()
        center, size = get_center_of_coordinates(ligaligner.new_residue)
        affinity = run_autodock_vina(lig_file, rec_file, f"opt-{candidate_structure}qt", center, size)
        # print(f"Affinity: {affinity}")
        os.chdir(current_path)
        return {
            "id": id,
            "candidate-name": candidate_structure,
            "prot-rmsd": rmsd,
            "affinity": affinity
        }

def test():
    s = Screener(
            'candidates',
            '6jtt.pdb',
            'A',
            701
        )
    results = Parallel(n_jobs=-1)(delayed(s.screen)(*x) for x in enumerate(s.candidates))
    all_res = {
        "id": [],
        "candidate-name": [],
        "prot-rmsd": [],
        "affinity": []
    }
    for res_dict in results:
        for k, v in res_dict.items():
            all_res[k].append(v)
    df = pd.DataFrame(all_res)
    df.to_csv('screening_resluts.csv', index=None, sep=',')

def write_output(results:dict):
    all_res = {
        "id": [],
        "candidate-name": [],
        "prot-rmsd": [],
        "affinity": []
    }
    for res_dict in results:
        for k, v in res_dict.items():
            all_res[k].append(v)
    df = pd.DataFrame(all_res)
    df.to_csv('screening_resluts.csv', index=None, sep=',')

if __name__ == '__main__':
    args = get_args()

    ligand_chain = args.ligand_chain
    ligand_resid = args.ligand_resid
    template_pdb = args.template_pdb

    structure_dir = args.structure_dir
    if structure_dir is not None:
        s = Screener(structure_dir, template_pdb, ligand_chain, ligand_resid, args.candidate_ligand_chain, args.candidate_ligand_resid)
        results = Parallel(n_jobs=-1)(delayed(s.screen)(*x) for x in enumerate(s.candidates))
        write_output(results)
    else:
        fasta_file = args.fasta_file
        active_chains = args.active_chains
        headers = read_fasta_headers(fasta_file)
        try:
            os.mkdir('modeller_out')
        except:
            print("the dir ./modeller_out/ exists.")
            # pass
        all_mutations = []
        all_seq_names = []
        for header in headers:
            seq_name, mutations = convert_name_to_mutations(header_line=header, active_chains=active_chains.split(","))
            all_mutations.append(mutations)
            all_seq_names.append(seq_name)

        mutated_structures = Parallel(n_jobs=-1)(delayed(build_mutation)(template_pdb, mutations, seq_name) for seq_name, mutations in zip(all_seq_names, all_mutations))
        s = Screener('modeller_out', template_pdb, ligand_chain, ligand_resid, args.candidate_ligand_chain, args.candidate_ligand_resid)
        results = Parallel(n_jobs=-1)(delayed(s.screen)(*x) for x in enumerate(s.candidates))
        write_output(results)