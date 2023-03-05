#!/usr/bin/env python

import os

PATH_TO_BRENDA = ''


def run_blast(querys, evalue_cutoff, identity_cutoff, threads, db_path=PATH_TO_BRENDA):
    """
    """
    blast_cmd  = f"blastp -query {querys} -db {db_path} -out {querys}_vs_brenda.txt -evalue {evalue_cutoff} -outfmt 6 -num_threads {threads}"
    os.system(blast_cmd)

    query_list = []
    hit_list = []
    blast_outfile = open(f'{querys}_vs_brenda.txt', 'r')
    for line in blast_outfile:
        query, hit = line.split("\t")[0,1]
        query_list.append(query)
        hit_list.append(hit)
    return query_list, hit_list


def ec_filter(query_list, hit_list, ec_number):
    """
    """
    passed_query = query_list.copy()
    for i, query in enumerate(query_list):
        anno = hit_list[i]
        anno_ec_number = anno.split("|")[0].replace(">EC:", "").split(".")
        if len(anno_ec_number) < 3:
            passed_query.pop(i)
        elif ec_number == "_".join(anno_ec_number):
            continue
        else:
            passed_query.pop(i)
    return set(passed_query)


def read_query(query_file):
    """
    """
    query_file = open(query_file, 'r')
    name_map = open('name_map.txt', 'w+')
    protein_dict = {}
    i = 1
    for line in query_file:
        if line.startswith(">"):
            i += 1
            name_map.write(f"protein_{i}\t{line}")
            protein_dict[f"protein_{i}"] = ""
        else:
            protein_dict[f"protein_{i}"] += line.strip().replace("*", "")
    name_map.close()
    query_file.close()
    return protein_dict

