#GASSER for proteins 20200203 build
# srtucture_prediction
print("""
        #######################################################
        #                      GASSER                         #
        #                                                     #
        #  Genome scAle Substrate Specific Enzyme pRediction  #
        #                                                     #
        #  Author: Sun Jinyuan, Ming Dengming*                #
        #                                                     #
        #######################################################
        
        """)
import os
import argparse

parser = argparse.ArgumentParser(
    description="""GASSER: Genome scAle Substrate Specific Enzyme pRediction""",
    usage='python3 GASSER.py proteins.fasta ligand.pdb')
parser.add_argument(
    "input",
    help=
    "input a seqence in fasta format and we highly suggest that you make the head linesimple",
    type=str)
parser.add_argument("ligand", help="ligand or substrate in pdb", type=str)
#parser.add_argument("filetype", help="the input file type, DNA or AA"
#,type=str)
parser.add_argument(
    "-nos",
    help=
    "Max Number Of Structures you wish to be predicted for each sequence, defult is 1",
    type=int,
    default=1)
parser.add_argument("-nt",
                    help="Number of Threads used to run BLAST, default is 1",
                    type=int,
                    default=32)
parser.add_argument(
    "-ic",
    help="the Identity Cutoff for structure prediction, default is 30",
    type=int,
    default=30)

args = parser.parse_args()

#inf = open(args.i)
infile = args.input
ligand = args.ligand
nos = args.nos
nt = args.nt
ic = args.ic
"""
print(str(nos),str(nt),(ic))
try:
    nos = int(args.nos)
except ValueError:
    nos = 1
    
try:
    nt = int(args.nt)
except ValueError:
    nos = 1
"""


def maindic(infile):
    in_file = open(infile)
    dic = {}
    for line in in_file:
        if line.startswith(">"):
            head = line.strip()
            dic[head] = ''
        else:
            dic[head] += line.strip().replace("*", "")
    return dic


def blast_pdb(query_file_name, nos, nt, ic):
    #a list for template
    target_lst = []
    #blastsearch
    blastout = os.popen(
        "psiblast -query " + query_file_name +
        " -db ../database/pdbaa/pdbaa -outfmt 6 -num_threads " + str(nt) +
        " -evalue 1e-5")
    #split blastout into list
    try:
        blastout = blastout.read()[:-2].split("\n")
        if len(blastout) < nos + 1:
            for line in blastout:
                identity = float(line.split("\t")[2])
                pdb = line.split("\t")[1]
                if identity > ic:
                    #download pdb templates
                    os.system("wget https://files.rcsb.org/download/" +
                              pdb.split("_")[0] + ".pdb")
                    target_lst.append(pdb)
        if len(blastout) > nos:
            blastout = blastout[0:nos]
            for line in blastout:
                identity = float(line.split("\t")[2])
                pdb = line.split("\t")[1]
                if identity > ic:
                    os.system("wget https://files.rcsb.org/download/" +
                              pdb.split("_")[0] + ".pdb")
                    target_lst.append(pdb)
        return target_lst
    except IndexError:
        return "No match"


def generate_seq(
    query_file_name
):  #prepare a sequence file of target meet the requirments of modeller
    query_fasta = open(query_file_name)
    seq = ''
    for line in query_fasta:
        if line.startswith(">"):
            head = line.replace(" ", "_").replace("\n", "").replace(
                ">", "").replace("-", "_").replace("|", "_")
        else:
            seq += line.replace("\n", "")
    sec_line = "sequence:" + head + ":::::::0.00: 0.00"
    with open(head + ".ali", "w+") as f:
        s = ">P1;" + head + "\n" + sec_line + "\n" + seq + "*"
        f.write(s)
    f.close()
    return head


def generate_align2d(query_file_name, target_code):
    target_name = target_code.split("_")[0]
    chain = target_code.split("_")[1]
    with open("align2d.py", "w+") as f:
        f.write("from modeller import *\n")
        f.write("env = environ()")
        f.write("\n")
        f.write("aln = alignment(env)")
        f.write("\n")
        f.write("mdl = model(env, file='" + target_name +
                "', model_segment=('FIRST:" + chain + "','LAST:" + chain +
                "'))")
        f.write("\n")
        f.write("aln.append_model(mdl, align_codes='" + target_code +
                "', atom_files='" + target_name + ".pdb')")
        f.write("\n")
        f.write("aln.append(file='" + head + ".ali', align_codes='" + head +
                "')")
        f.write("\n")
        f.write("aln.align2d()")
        f.write("\n")
        f.write("aln.write(file='" + head + "_" + target_code +
                ".ali', alignment_format='PIR')")
        f.write("\n")
        f.write("aln.write(file='" + head + "_" + target_code +
                ".pap', alignment_format='PAP')")
    f.close()
    ali_code = head + "_" + target_code
    return ali_code


def generatr_model(ali_code, target_code, head):
    with open("model.py", "w+") as f:
        f.write("from modeller import *")
        f.write("\n")
        f.write("from modeller.automodel import *")
        f.write("\n")
        f.write("env = environ()")
        f.write("\n")
        f.write("a = automodel(env, alnfile='" + ali_code + ".ali',knowns='" +
                target_code + "', sequence='" + head +
                "',assess_methods=(assess.DOPE,assess.GA341))")
        f.write("\n")
        f.write("a.starting_model = 1")
        f.write("\n")
        f.write("a.ending_model = 1")
        f.write("\n")
        #f.write("a.get_model_filename(sequence='"+head+"',)")
        f.write("a.make()")
    f.close()


def docking(receptor, ligand):
    os.system(
        "pythonsh ~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "
        + str(ligand))
    os.system(
        "pythonsh ~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r "
        + str(receptor))
    os.system(
        "pythonsh ~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -r "
        + receptor + "qt -l " + ligand + "qt")
    os.system("autogrid4 -p " + receptor.replace(".pdb", "") + ".gpf -l " +
              receptor.replace(".pdb", "") + ".grd.log")
    os.system("prankp predict -f " + receptor + " -o ./")
    pred = open(receptor + "_predictions.csv")
    for line in pred:
        lst = line.split(",")
        try:
            if float(lst[2]) > 10:
                rank = lst[1]
                centre_x = lst[5]
                centre_y = lst[6]
                centre_z = lst[7]
                with open(
                        receptor.replace(".pdb", "") + "_pocket" + rank +
                        "vina.cfg", "w+") as f:
                    f.write("receptor = " + receptor + "qt")
                    f.write("\n")
                    f.write("ligand = " + ligand + "qt")
                    f.write("\n")
                    f.write("center_x = " + str(centre_x))
                    f.write("\n")
                    f.write("center_y = " + str(centre_y))
                    f.write("\n")
                    f.write("center_z = " + str(centre_z))
                    f.write("\n")
                    f.write("""size_x = 30 
                               size_y = 30
                               size_z = 30""")
                    f.write("\n")
                    f.write("energy_range = 4")
                    f.write("\n")
                    f.write("out = " + receptor.replace(".pdb", "") +
                            "_pocket" + rank + "_" +
                            ligand.replace(".pdb", ".pdbqt"))
                f.close()
                vina = os.popen("vina --config " +
                                receptor.replace(".pdb", "") + "_pocket" +
                                rank + "vina.cfg")
                vina_out = vina.read().split("-+-")[-1].split("\n")[1:][:-1]
                with open(
                        receptor.replace(".pdb", "") + "_pocket" + rank + "_" +
                        ligand.replace(".pdb", ".out"), "w+") as dockout:
                    dockout.write(receptor + "_pocket" + " " + rank +
                                  " docking:\nmode\taffinity (kcal/mol)\trmsd")
                    for x in vina_out:
                        xlst = x.split()
                        if float(xlst[1]) < -4:
                            dockout.write("\t".join(str(a) for a in xlst))
                dockout.close()
        except ValueError:
            continue


#receptor = head + ".B99990001.pdb"
seqdic = maindic(infile)
for query_name in seqdic:
    query_file_name = query_name.replace(">", "") + ".fasta"
    #write the sequence in fasta format
    with open(query_file_name, "w+") as seqfile:
        seqfile.write(query_name)
        seqfile.write("\n")
        seqfile.write(seqdic[query_name])
    seqfile.close()
    #print(query_name+"\n"+seqdic[query_name],file=seqfile)
    target_lst = blast_pdb(query_file_name, nos, nt, ic)
    if target_lst == "No match":
        print(str(query_name) + "No match")
    else:
        head = generate_seq(query_file_name)
        for x in target_lst:
            outlst = []
            ali_code = generate_align2d(query_file_name, x)
            generatr_model(ali_code, x, head)
            os.system("/opt/modeller923/bin/mod9.23 align2d.py")
            os.system("/opt/modeller923/bin/mod9.23 model.py")
            os.system("mv " + head + ".B99990001.pdb " + head + "_" + x +
                      ".pdb")
            outlst.append(head + "_" + x + ".pdb")

        # docking
        py_wd = "~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/"

        ligand = "ZEN.pdb"
        for receptor in outlst:
            docking(receptor, ligand)
