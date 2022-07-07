#!/usr/bin/env python3
# output should follow the convention detailed here:
# https://docs.aws.amazon.com/neptune/latest/userguide/bulk-load-tutorial-format-gremlin.html

# import sys
import re

# import json
import urllib.parse

# import os.path
import logging
import argparse

# from numpy import append
import pandas as pd


# def read_external_fasta(file_name):
#     fasta = {}
#     with open(file_name) as fp:
#         lines = fp.read().splitlines()
#         current_fasta = ""
#         for line in lines:
#             if line[0] == ">":
#                 current_fasta = line[1:]
#             elif len(line) > 5:
#                 fasta[current_fasta] = line
#     return fasta


def parse_gff_record(line):
    segment = re.split("\t", line)
    attributesString = re.split("\n|;", segment[8])
    attributes = {}
    for s in attributesString:
        e = s.split("=")
        if len(e) >= 2:
            attributes.update({e[0]: urllib.parse.unquote(e[1])})
    if len(attributes) == 0:
        attributes = {
            "ID": f"{segment[0]}-{segment[2]}-{str(int(segment[3])- 1)}-{segment[4]}",
            "product": segment[2],
        }
    record = {
        "seqName": segment[0],
        "source": segment[1],
        "feature": segment[2],
        "start": int(segment[3]) - 1,
        "end": int(segment[4]),
        "score": segment[5],
        "strand": segment[6],
        "frame": segment[7],
        "attribute": attributes,
    }
    return record


def get_gene_node(line: str) -> dict:
    record = parse_gff_record(line)

    if "ID" not in record["attribute"]:
        logging.info(
            f"Could not extract relevant GENE information from this line: {line}"
        )
        logging.info(f"Here is what was parsed: {record}")
        return None

    label = "GENE"
    if record["feature"]:
        label = f"{label};{record['feature'].upper()}"

    node = {
        "~id": record["attribute"]["ID"],
        "~label": label,
        "type:String": record["feature"],
        "contig:String": record["seqName"],
        "source:String": record["source"],
        "location_start:Int": record["start"],
        "location_end:Int": record["end"],
        "location_strand:String": record["strand"],
    }

    if "product" in record["attribute"]:
        node["product:String"] = record["attribute"]["product"]

    if "locus_tag" in record["attribute"]:
        node["locus_tag:String"] = record["attribute"]["locus_tag"]

    return node


def parse_tabular_pfam_record(line):
    columns = [
        "gene_oid",
        "gene_length",
        "query_start",
        "query_end",
        "subj_start",
        "subj_end",
        "evalue",
        "bit_score",
        "pfam_id",
        "pfam_name",
        "pfam_length",
    ]
    segment = re.split("\t", line)
    record = {}
    for idx, s in enumerate(segment):
        column = columns[idx]
        record.update({column: urllib.parse.unquote(s)})

    return record


def get_pfam_node(line: str) -> tuple:
    record = parse_tabular_pfam_record(line)

    if "pfam_id" not in record:
        logging.info(
            f"Could not extract relevant PFAM information from this line: {line}"
        )
        logging.info(f"Here is what was parsed: {record}")
        return None, None, None

    node_id = f"{record['gene_oid']}__{record['pfam_id']}__{record['query_start']}"

    node = {
        "~id": node_id,
        "~label": "PFAM;HMM",
        "location_start:Int": record["query_start"],
        "location_end:Int": record["query_end"],
        "evalue:Double": record["evalue"],
        "bit_score:Double": record["bit_score"],
        "pfam_id:String": record["pfam_id"],
        "name:String": record["pfam_name"],
        "pfam_length:Int": record["pfam_length"],
    }

    gene_has_hmm = {
        "~id": f"{record['gene_oid']}_r_{node_id}",
        "~from": f"{record['gene_oid']}",
        "~to": node_id,
        "~label": "HAS_HMM",
    }

    return record["gene_oid"], node, gene_has_hmm


# def parse_tabular_tigrfam_record(line):
#     columns = [
#         "gene_oid",
#         "gene_length",
#         "query_start",
#         "query_end",
#         "evalue",
#         "bit_score",
#         "tigrfam_id",
#         "tigrfam_name",
#     ]
#     segment = re.split("\t", line)
#     record = {}
#     for idx, s in enumerate(segment):
#         column = columns[idx]
#         record.update({column: urllib.parse.unquote(s)})

#     return record


# def get_tigrfam_node(line):
#     content = ""
#     record = parse_tabular_tigrfam_record(line)

#     if "tigrfam_id" not in record:
#         logging.info(
#             f"Could not extract relevant TIGRFAM information from this line: {line}"
#         )
#         logging.info(f"Here is what was parsed: {record}")
#         return None, None, None

#     node_id = f"{record['gene_oid']}__{record['tigrfam_id']}__{record['query_start']}"

#     node = {
#         "type": "node",
#         "id": node_id,
#         "labels": ["TIGRFAM", "HMM"],
#         "properties": {
#             "location_start": record["query_start"],
#             "location_end": record["query_end"],
#             "evalue": record["evalue"],
#             "bit_score": record["bit_score"],
#             "tigrfam_id": record["tigrfam_id"],
#             "name": record["tigrfam_name"],
#             # "gene_length": record.gene_length,
#             # "subj_start": record.subj_start,
#             # "subj_end": record.subj_end,
#         },
#     }
#     content = json.dumps(node) + "\n"

#     relation = {
#         "id": f"{record['gene_oid']}_r_{record['tigrfam_id']}",
#         "type": "relationship",
#         "label": "HAS_HMM",
#         "start": {"id": f"{record['gene_oid']}", "labels": ["GENE"]},
#         "end": {"id": node_id, "labels": ["TIGRFAM"],},
#     }

#     content += json.dumps(relation) + "\n"

#     return record["gene_oid"], node["id"], content


def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def ingest_gff_file(gff_file: str) -> tuple:
    vertices = list()
    edges = list()
    gene_neighbors = dict()
    with open(gff_file) as gff:
        for line in nonblank_lines(gff):
            if re.match("##FASTA", line):
                break
            if line[0] == "#":
                continue
            gene_node = get_gene_node(line)
            if gene_node is None:
                continue

            vertices.append(gene_node)

            contig_id = gene_node["contig:String"]
            # If genes are on the same contig, connect them with a "NEXT" relationship
            if contig_id not in gene_neighbors:
                gene_neighbors.update({contig_id: list()})

            gene_neighbors[contig_id].append(gene_node["~id"])

    edges += generate_next_relationships(gene_neighbors, "GENE")

    return vertices, edges


def ingest_pfam_file(pfam_file):
    vertices = list()
    edges = list()
    neighbors = dict()
    with open(pfam_file) as pfam:
        for idx, line in enumerate(pfam.readlines()):
            # skip the header line
            if idx == 0:
                continue

            line = line.strip()

            gene_id, pfam_node, gene_pfam_edge = get_pfam_node(line)
            if gene_id is not None:
                vertices.append(pfam_node)
                edges.append(gene_pfam_edge)

                if gene_id not in neighbors:
                    neighbors.update({gene_id: list()})

                neighbors[gene_id].append(pfam_node["~id"])

    edges += generate_next_relationships(neighbors, "PFAM")
    return vertices, edges


# def ingest_tigrfam_file(tigrfam_file):
#     content = ""
#     neighbors = dict()
#     with open(tigrfam_file) as tigrfam:
#         for idx, line in enumerate(tigrfam.readlines()):
#             # skip the header line
#             if idx == 0:
#                 continue

#             line = line.strip()

#             gene_id, tigrfam_node_id, tigrfam_content = get_tigrfam_node(line)
#             if gene_id is not None:
#                 content += tigrfam_content

#                 if gene_id not in neighbors:
#                     neighbors.update({gene_id: list()})

#                 neighbors[gene_id].append(tigrfam_node_id)

#     content += generate_next_relationships(neighbors, "TIGRFAM")
#     return content


def generate_next_relationships(neighbors: dict, label: str) -> list:
    relationships = list()
    for _, child_list in neighbors.items():
        total_children = len(child_list)
        idx = 0
        while idx <= (total_children - 2):
            this_child = child_list[idx]
            next_child = child_list[idx + 1]
            relationships.append(
                {
                    "~id": f"{this_child}_r_{next_child}",
                    "~from": f"{this_child}",
                    "~to": f"{next_child}",
                    "~label": "NEXT",
                }
            )
            idx += 1

    return relationships


def usage():
    parser = argparse.ArgumentParser(description="convert gff file format to json")

    parser.add_argument("--gff", help="gff file")
    parser.add_argument("--pfam", help="JGI PFam tabular output")
    # parser.add_argument("--tigrfam", help="JGI TIGRFam tabular output")
    parser.add_argument("--prefix", help="specify the output prefix")

    args = parser.parse_args()

    return args


# def save_to_json(content, json_file, append=False):
#     if append:
#         output = open(json_file, "w+")
#     else:
#         output = open(json_file, "w")

#     output.write(content)
#     output.close()

#     return


def main():
    args = usage()

    gene_nodes_out = f"{args.prefix}.genes_nodes.csv.gz"
    gene_edges_out = f"{args.prefix}.genes_edges.csv.gz"

    pfam_nodes_out = f"{args.prefix}.pfam_nodes.csv.gz"
    pfam_edges_out = f"{args.prefix}.pfam_edges.csv.gz"

    # tigrfam_json = f"{args.prefix}.tigrfam.json"
    logging.info("Ingesting Gene data ...")
    gene_nodes, gene_edges = ingest_gff_file(args.gff)
    pd.DataFrame(gene_nodes).to_csv(gene_nodes_out, index=False, compression="gzip")
    pd.DataFrame(gene_edges).to_csv(gene_edges_out, index=False, compression="gzip")
    logging.info(f"\tSaving {gene_nodes_out} and {gene_edges_out}")

    logging.info("Ingesting PFam data ...")
    pfam_nodes, pfam_edges = ingest_pfam_file(args.pfam)
    pd.DataFrame(pfam_nodes).to_csv(pfam_nodes_out, index=False, compression="gzip")
    pd.DataFrame(pfam_edges).to_csv(pfam_edges_out, index=False, compression="gzip")
    logging.info(f"\tSaving {pfam_nodes_out} and {pfam_edges_out}")

    # tigrfam_content = ingest_tigrfam_file(args.tigrfam)
    # save_to_json(tigrfam_content, tigrfam_json)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    main()
