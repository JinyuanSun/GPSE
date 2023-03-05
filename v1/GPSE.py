#!/usr/bin/python3
import os
import argparse
from GPSEcfg import ADT_PATH, PDBAA_PATH, MODELLER_PATH, BRENDA_PATH

print(
    """
        #########################################################
        #                         GPSE                          #
        #                                                       #
        # Genome-scale Prediction of Substrate-specific Enzymes #
        #                                                       #
        #  Author: Sun Jinyuan, Xia Yan, Ming Dengming*         #
        #                                                       #
        #########################################################
        
        """
)
def getparser():
    parser = argparse.ArgumentParser(
        description="""GPSE: Genome-scale Prediction of Substrate-specific Enzymes""",
        usage="python3 GPSE.py proteins.fasta ligand.pdb EC_number",
    )
    parser.add_argument(
        "input",
        help="input a seqence in fasta format and we highly suggest that you make the head linesimple",
        type=str,
    )
    parser.add_argument("ligand", help="ligand or substrate in pdb", type=str)
    # parser.add_argument("filetype", help="the input file type, DNA or AA"
    # ,type=str)
    parser.add_argument(
        "EC_number",
        help="the EC_number in form of x_x_x (e.g. 3_1_1 for 3.1.1. carboxylic-ester hydrolase)",
        type=str,
    )

    parser.add_argument(
        "-nos",
        help="Max Number Of Structures you wish to be predicted for each sequence, defult is 1",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-nt",
        help="Number of Threads used to run BLAST, default is 8",
        type=int,
        default=8,
    )
    parser.add_argument(
        "-evalue",
        help="E-Value cutoff for structure template, default is 1e-5",
        type=str,
        default=str(1e-5),
    )
    parser.add_argument(
        "-ic",
        help="the Identity Cutoff for structure prediction, default is 30",
        type=int,
        default=30,
    )

    args = parser.parse_args()
    return args


def maindic(infile):
    in_file = open(infile)
    dic = {}
    for line in in_file:
        if line.startswith(">"):
            head = line.strip()
            dic[head] = ""
        else:
            dic[head] += line.strip().replace("*", "")
    return dic


def blast_pdb(query_file_name, nos, nt, ic):
    # a list for template
    target_lst = []
    # blastsearch
    blastout = os.popen(
        "psiblast -query "
        + query_file_name
        + " -db "+PDBAA_PATH+" -outfmt 6 -num_threads "
        + str(nt)
        + " -evalue "
        + str(evalue)
    )
    # split blastout into list
    try:
        blastout = blastout.read()[:-2].split("\n")
        if len(blastout) < nos + 1:
            for line in blastout:
                identity = float(line.split("\t")[2])
                pdb = line.split("\t")[1]
                if identity > ic:
                    # download pdb templates
                    # os.system("wget https://files.rcsb.org/download/" +
                    # pdb.split("_")[0] + ".pdb")
                    target_lst.append(pdb)
        if len(blastout) > nos:
            blastout = blastout[0:nos]
            for line in blastout:
                identity = float(line.split("\t")[2])
                pdb = line.split("\t")[1]
                if identity > ic:
                    # os.system("wget https://files.rcsb.org/download/" +
                    # pdb.split("_")[0] + ".pdb")
                    target_lst.append(pdb)
        return target_lst
    except IndexError:
        return "No match"


def generate_seq(query_file_name):  
    # prepare a sequence file of target meet the requirments of modeller
    query_fasta = open(query_file_name)
    seq = ""
    for line in query_fasta:
        if line.startswith(">"):
            head = (
                line.replace(" ", "_")
                .replace("\n", "")
                .replace(">", "")
                .replace("-", "_")
                .replace("|", "_")
            )
        else:
            seq += line.replace("\n", "")
    sec_line = "sequence:" + head + ":::::::0.00: 0.00"
    with open(head + ".ali", "w+") as f:
        s = ">P1;" + head + "\n" + sec_line + "\n" + seq + "*"
        f.write(s)
    f.close()
    return head


def generate_align2d(query_file_name, target_code, head):
    target_name = target_code.split("_")[0]
    chain = target_code.split("_")[1]
    with open("align2d.py", "w+") as f:
        f.write("from modeller import *\n")
        f.write("env = environ()")
        f.write("\n")
        f.write("aln = alignment(env)")
        f.write("\n")
        f.write(
            "mdl = model(env, file='"
            + target_name
            + "', model_segment=('FIRST:"
            + chain
            + "','LAST:"
            + chain
            + "'))"
        )
        f.write("\n")
        f.write(
            "aln.append_model(mdl, align_codes='"
            + target_code
            + "', atom_files='"
            + target_name
            + ".pdb')"
        )
        f.write("\n")
        f.write("aln.append(file='" + head + ".ali', align_codes='" + head + "')")
        f.write("\n")
        f.write("aln.align2d()")
        f.write("\n")
        f.write(
            "aln.write(file='"
            + head
            + "_"
            + target_code
            + ".ali', alignment_format='PIR')"
        )
        f.write("\n")
        f.write(
            "aln.write(file='"
            + head
            + "_"
            + target_code
            + ".pap', alignment_format='PAP')"
        )
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
        f.write(
            "a = automodel(env, alnfile='"
            + ali_code
            + ".ali',knowns='"
            + target_code
            + "', sequence='"
            + head
            + "',assess_methods=(assess.DOPE,assess.GA341))"
        )
        f.write("\n")
        f.write("a.starting_model = 1")
        f.write("\n")
        f.write("a.ending_model = 1")
        f.write("\n")
        # f.write("a.get_model_filename(sequence='"+head+"',)")
        f.write("a.make()")
    f.close()


def docking(receptor, ligand):
    os.system(
        "pythonsh "+ADT_PATH+"prepare_ligand4.py -l "
        + str(ligand)
    )
    os.system(
        "pythonsh "+ADT_PATH+"prepare_receptor4.py -A checkhydrogens -r "
        + str(receptor)
    )
    #os.system(
    #    "pythonsh ~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -r "
    #    + receptor
    #    + "qt -l "
    #    + ligand
    #    + "qt"
    #)
    #os.system(
    #    "autogrid4 -p "
    #    + receptor.replace(".pdb", "")
    #    + ".gpf -l "
    #    + receptor.replace(".pdb", "")
    #    + ".grd.log"
    #)
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
                    receptor.replace(".pdb", "") + "_pocket" + rank + "vina.cfg", "w+"
                ) as f:
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
                    f.write(
                        """size_x = 15 
                           size_y = 15
                           size_z = 15"""
                    )
                    f.write("\n")
                    f.write("energy_range = 4")
                    f.write("\n")
                    f.write(
                        "out = "
                        + receptor.replace(".pdb", "")
                        + "_pocket"
                        + rank
                        + "_"
                        + ligand.replace(".pdb", ".pdbqt")
                    )
                f.close()
                vina = os.popen(
                    "vina --config "
                    + receptor.replace(".pdb", "")
                    + "_pocket"
                    + rank
                    + "vina.cfg"
                )
                vina_out = vina.read().split("-+-")[-1].split("\n")[1:][:-1]
                with open(
                    receptor.replace(".pdb", "")
                    + "_pocket"
                    + rank
                    + "_"
                    + ligand.replace(".pdb", ".out"),
                    "w+",
                ) as dockout:
                    dockout.write(
                        receptor
                        + "_pocket"
                        + " "
                        + rank
                        + " docking:\nmode\taffinity (kcal/mol)\trmsd\n"
                    )
                    for x in vina_out:
                        xlst = x.split()
                        print(xlst)
                        if float(xlst[1]) < -4:
                            dockout.write("\t".join(str(a) for a in xlst))
                            dockout.write("\n")
                dockout.close()
        except ValueError:
            continue


def run_(query_name):
    start_path = os.getcwd()
    query_file_name = query_name.replace(">", "") + ".fasta"
    dir_name = query_name.replace(">", "")
    # make a new dir for files
    os.system("mkdir " + dir_name)
    # go to the new dir
    os.chdir(dir_name)
    # write the sequence in fasta format
    with open(query_file_name, "w+") as seqfile:
        seqfile.write(query_name)
        seqfile.write("\n")
        seqfile.write(seqdic[query_name])
    seqfile.close()
    # run BLAST for template target_lst contains
    template_lst = blast_pdb(query_file_name, nos, nt, ic)
    if template_lst == "No match":
        print(str(query_name) + " No match")
    else:
        head = generate_seq(query_file_name)
        outlst = []
        wgetcode = ""
        for x in template_lst:
            try:
                os.mkdir("structure")
                os.chdir("structure")
                # os.system("cp ../" + head + ".ali ./")
                # os.system("mv ../"+x.split("_")[0]+".pdb ./")
                wgetcode = os.system(
                    "wget -q https://files.rcsb.org/download/"
                    + x.split("_")[0]
                    + ".pdb"
                )
                # check if the pdb file exists
                if str(wgetcode) == "2048":
                    print(x + " maybe in cif format, please check it manually")
                    os.chdir("..")
                    break
                else:
                    os.system("cp ../" + head + ".ali ./")

                    ali_code = generate_align2d(query_file_name, x, head)
                    generatr_model(ali_code, x, head)
                    os.system(MODELLER_PATH+" align2d.py")
                    os.system(MODELLER_PATH+" model.py")
                    os.system(
                        "mv " + head + ".B99990001.pdb " + head + "_" + x + ".pdb"
                    )
                    outlst.append(head + "_" + x + ".pdb")
                    # os.system("rm *.py")
                    # os.system("rm *.log")
                    os.chdir("..")
            except FileExistsError:
                # os.system("mv " + x + " " + x + "_bak")
                # os.mkdir(x)
                os.chdir("structure")
                # os.system("cp ../" + head + ".ali ./")
                # os.system("cp ../"+x.split("_")[0]+".pdb ./")
                wgetcode = os.system(
                    "wget -q https://files.rcsb.org/download/"
                    + x.split("_")[0]
                    + ".pdb"
                )
                if str(wgetcode) == "2048":
                    print(x + " maybe in cif format, please check it manually")
                    os.chdir("..")
                    break
                else:
                    os.system("cp ../" + head + ".ali ./")

                    ali_code = generate_align2d(query_file_name, x, head)
                    generatr_model(ali_code, x, head)
                    os.system(MODELLER_PATH+" align2d.py")
                    os.system(MODELLER_PATH+" model.py")
                    os.system(
                        "mv " + head + ".B99990001.pdb " + head + "_" + x + ".pdb"
                    )
                    outlst.append(head + "_" + x + ".pdb")
                    os.chdir("..")

            # docking
        # py_wd = "~/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/"
        # ligand = "ZEN.pdb"
        for receptor in outlst:
            if wgetcode != "2048":
                os.mkdir(receptor)
                os.chdir(receptor)
                os.system("cp " + start_path + "/" + ligand + " ./")
                os.system("cp ../structure/" + receptor + " ./")
                docking(receptor, ligand)
                os.chdir("..")
            else:
                print("No docking preformed!")
    os.chdir("..")
    return str(query_name) + " finished!"


def ec_anno_filter(dic):
    out_dic = {}
    with open("EC_assignmet.out", "w+") as ecf:
        for query in dic:
            idic = {
                k: v
                for k, v in sorted(
                    dic[query].items(), key=lambda item: item[1], reverse=True
                )
            }
            # print(idic)
            eclist = list(idic.keys())
            # print(eclist)
            outlist = []
            i = 0
            # print(idic)
            # print(len(eclist))

            for i in range(len(eclist)):
                # print("yes")
                # print(i)
                try:
                    if idic[eclist[i]] > idic[eclist[i + 1]] > 1:
                        outlist.append(eclist[i])
                        break
                        i += 1
                    if idic[eclist[i]] == 1:
                        i = i + 1
                        break
                    if idic[eclist[i]] == idic[eclist[i + 1]] > 1:
                        outlist.append(eclist[i])
                        i = i + 1
                except IndexError:
                    if idic[eclist[i]] > 1:
                        outlist.append(eclist[i])
                        break
                # else:
                # i = i + 1
            # print(query+ ",".join(outlist))
            out_dic[">" + str(query)] = ",".join(outlist)
            ecf.write(query + "\t" + ",".join(outlist) + "\n")
            # print(query+"\t"+",".join(outlist))
    ecf.close()
    return out_dic


def ec_assignment(querys):
    os.system(
        "blastp -query "
        + querys
        + " -db "+BRENDA_PATH+" -out "
        + querys
        + ".out -evalue 1e-5 -outfmt 6 -num_threads "
        + str(nt)
    )
    dic = {}
    query_list = []
    outfile = open(querys + ".out")
    for line in outfile:
        lst = line.split("\t")
        qname = lst[0]
        hname = "_".join(lst[1].split(".")[0:3])
        if qname in query_list:
            try:
                dic[qname][hname] += 1
            except KeyError:
                dic[qname][hname] = 1
        else:
            dic[qname] = {hname: 1}
            query_list.append(qname)
    return dic
    # return sorted(dic['tr|Q8NKB0|Q8NKB0_BIOOC'], key=dic['tr|Q8NKB0|Q8NKB0_BIOOC'].get, reverse=True)


# inf = open(args.i)
if __name__ == "__main__":
    args = getparser()
    infile = args.input
    ligand = args.ligand
    nos = args.nos
    nt = args.nt
    ic = args.ic
    EC_number = args.EC_number
    evalue = args.evalue
    ec_dic = ec_assignment(infile)
    seqdic = maindic(infile)

    print("file read!")
    out_dic = ec_anno_filter(ec_dic)  # {">query":"ec_number"}
    for query_name in out_dic:
        ec_list = out_dic[query_name].split(",")
        if EC_number in ec_list:
            run_(query_name)
