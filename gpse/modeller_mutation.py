#!/usr/bin/env python
import sys
import os
from Bio.PDB.Polypeptide import one_to_three
from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched

# from https://bougui505.github.io/2016/09/06/mutate_a_single_residue_in_a_pdb_with_python_and_modeller.html

def optimize(atmsel, sched):
    # conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    # md
    refine(atmsel)
    cg = conjugate_gradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


# molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return="FINAL")
    init_vel = True
    for (its, equil, temps) in (
        (200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
        (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0)),
    ):
        for temp in temps:
            md.optimize(
                atmsel,
                init_velocities=init_vel,
                temperature=temp,
                max_iterations=its,
                equilibrate=equil,
            )
            init_vel = False


# use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    s = selection(mdl1)
    for typ in ("stereo", "phi-psi_binormal"):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ("omega", "chi1", "chi2", "chi3", "chi4"):
        rsr.make(
            s,
            restraint_type=typ + "_dihedral",
            spline_range=4.0,
            spline_dx=0.3,
            spline_min_points=5,
            aln=aln,
            spline_on_site=True,
        )


def convert_name_to_mutations(header_line: str, active_chains=["A", "B"]):
    # >mut_5_340A_341G
    converted_mutations = []
    seq_name = "_".join(header_line.split("_")[:2]).replace(">", "")
    mutations = header_line.split("_")[2:]
    for chain in active_chains:
        for mutation in mutations:
            pos = mutation[:-1]
            one_letter = mutation[-1]
            converted_mutations.append(f"{pos} {one_to_three(one_letter)} {chain}")
    return seq_name, converted_mutations


def build_mutation(modelname: str, mutations: list, seq_name: str, output_dir:str='modeller_out'):
    log.verbose()

    env = environ(rand_seed=-49837)

    env.io.hetatm = True

    env.edat.dynamic_sphere = False

    env.edat.dynamic_lennard = True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    env.libs.topology.read(file="$(LIB)/top_heav.lib")

    env.libs.parameters.read(file="$(LIB)/par.lib")

    mdl1 = model(env, file=modelname)
    ali = alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    for mutation in mutations:
        respos, restyp, chain = mutation.split()
        s = selection(mdl1.chains[chain].residues[respos])
        s.mutate(residue_type=restyp)

    ali.append_model(mdl1, align_codes=modelname)

    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    mdl1.transfer_xyz(ali)

    mdl1.build(initialize_xyz=False, build_method="INTERNAL_COORDINATES")

    mdl2 = model(env, file=modelname)

    mdl1.res_num_from(mdl2, ali)

    mdl1.write(file=f"{seq_name}.tmp")
    # time.sleep(0.1)
    # try:
    mdl1.read(file=f"{seq_name}.tmp")
    # except:
    #     time.sleep(1)
    #     mdl1.read(file=modelname + restyp + respos + ".tmp")

    make_restraints(mdl1, ali)

    mdl1.env.edat.nonbonded_sel_atoms = 1

    sched = autosched.loop.make_for_model(mdl1)

    s = selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()

    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms = 2
    optimize(s, sched)

    mdl1.env.edat.nonbonded_sel_atoms = 1
    optimize(s, sched)

    s.energy()
    
    outname = os.path.join(output_dir, modelname.split("/")[-1].replace(".pdb", f"_{seq_name}.pdb"))
    mdl1.write(file=outname)

    os.remove(f"{seq_name}.tmp")
    return outname


if __name__ == "__main__":
    model_name, header_line, out_dir = sys.argv[1:]
    seq_name, mutations = convert_name_to_mutations(header_line)
    build_mutation(modelname=model_name, mutations=mutations, seq_name=seq_name, output_dir=out_dir)
