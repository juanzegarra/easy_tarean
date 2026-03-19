from Bio import SeqIO
import re
import pandas as pd
import argparse
from collections import defaultdict
from bs4 import BeautifulSoup
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        description="Easy post-processing of TAREAN html report")
    parser.add_argument("html_file", help="Path to the input html file")
    parser.add_argument("-o", type=str, default="output", help="Base name for output directory and files")
    parser.add_argument("--g_size", type=int, default=180000000, help="Genome size for estimated copy number calculation (default for average Drosophila genome size)")
    parser.add_argument("--plot", action="store_true", help="Plot clusters by size (number of reads) and annotation")
    parser.add_argument("--repeatM_file", type=str, help="path to custom annotation file")
    parser.add_argument("--custom_anno", action="store_true", help="Enable overwriting of TAREAN annotation for a external annotation file,\n more information about annotation file format on read.me")
    parser.add_argument("--custom_sat", action="store_true", help="Enable use of custom satDNA fasta file that annotates (note that it must be a exact sequence match)")
    parser.add_argument("--fasta_sats", type=str, help="path to custom satDNA file")
    return parser.parse_args()

class TAREAN_post_processing:
    def __init__(self, html_file, output_base, g_size):
        self.html_file = html_file
        self.output_base = output_base
        self.g_size = g_size
        #set for excluding other types of headers
        self.types = {" Putative satellites (high confidence)", " Putative satellites (low confidence)", " Other"}
        self.clusters, self.metrics = self.parse_html(html_file)

    def parse_html(self, html_file):
        metrics = defaultdict(int)
        with open(html_file, "r") as f:
            html_content = f.read()
            soup = BeautifulSoup(html_content, "html.parser")
            number_of_reads = int(re.search(r'Number of input reads: (\d+)', html_content).group(1))
            metrics["Number of input reads"] = number_of_reads
            analyzed_reads = int(re.search(r'Number of analyzed reads: (\d+)', html_content).group(1))
            metrics["Number of analyzed reads"] = analyzed_reads
            clusters = []
            #iterates over each in h2 header that delimitates tables
            for h2 in soup.find_all("h2"):
                if h2.text in self.types:
                    table = h2.find_next("tbody")
                    #excludes header line of the table
                    for row in table.find_all("tr", class_=lambda x: x != "firstline"):
                        cluster = defaultdict(int)
                        #divides row and extracts relevant cells
                        c = row.find("td", class_="firstcolumn")
                        cells = row.find_all("td", class_="cellinside")
                        cluster["Cluster"] = c.text.strip()
                        proportion = cells[2].text.strip()
                        number_reads = cells[3].text.strip()
                        sat_prob = cells[4].text.strip()
                        consensus = cells[6].text.strip()
                        C_index = cells[9].text.strip()
                        P_index = cells[10].text.strip()
                        annotation_cell = cells[15]

                        #extracts highest ranking annotation
                        clean_annotation = annotation_cell.find_all("b")
                        max_annotation = defaultdict(int)
                        annotation = ""
                        for i in clean_annotation:
                            if len(clean_annotation) <= 1:
                                annotation = clean_annotation[0].text.split("%")[1]
                            elif not i.text:
                              annotation = "NA"
                            else:
                                annotation_value = float(i.text.split("%")[0])
                                annotation_name = i.text.split("%")[1]
                                max_annotation[annotation_name] = annotation_value
                        if max_annotation:
                            max_annotation = max(max_annotation.items(), key=lambda x: x[1])
                            if max_annotation[1] >= 1:
                                annotation = max_annotation[0]

                        cluster["Proportion"] = proportion
                        cluster["Number of reads"] = number_reads
                        cluster["SAT prob"] = sat_prob
                        cluster["consensus"] = consensus
                        cluster["C-index"] = C_index
                        cluster["P-index"] = P_index
                        cluster["Annotation"] = annotation
                        cluster["Classification"] = h2.text.strip()
                        clusters.append(cluster)

        return clusters, metrics
    #corrects number of analyzed reads by excluding contamination
    def correct_metrics(self, metrics, clusters):
        contamination = 0
        for cluster in clusters:
            if cluster["Annotation"] == "Contamination":
                contamination += int(cluster["Number of reads"])
        self.metrics["Number of analyzed reads"] = int(metrics["Number of analyzed reads"]) - contamination

    def estimate_copy_number(self, clusters, metrics, g_size):
        for cluster in clusters:
            cluster["Proportion"] = (float(cluster["Proportion"]) / self.metrics["Number of analyzed reads"]) * 100
            consensus_size = len(cluster["consensus"])
            if consensus_size > 0:
                estimated_copy_number = ((g_size * float(cluster["Number of reads"])) / 100) / consensus_size
            else:
                estimated_copy_number = 0
            cluster["Estimated copy number"] = estimated_copy_number

    def write_files(self, clusters, metrics):
        with open(f"{self.output_base}.tsv", "w") as output_file, open(f"{self.output_base}.fasta", "w") as fasta_out:
            output_file.write("Cluster\tProportion\tNumber of reads\tSAT prob\tconsensus\tC-index\tP-index\tAnnotation\tClassification\tEstimated copy number\n")
            for cluster in clusters:
                output_file.write(f"{cluster['Cluster']}\t{cluster['Proportion']}\t{cluster['Number of reads']}\t{cluster['SAT prob']}\t{cluster['consensus']}\t{cluster['C-index']}\t{cluster['P-index']}\t{cluster['Annotation']}\t{cluster['Classification']}\t{cluster['Estimated copy number']}\n")
                if cluster["consensus"] != "":
                  fasta_out.write(f">{cluster['Cluster']}\t{cluster['Classification']}\n{cluster['consensus']}\n")

    def overwrited_annotation(self, repeatM_file, clusters):
        with open(repeatM_file, "r") as infile:
            custom_annotation = defaultdict(str)
            for line in infile:
                parts = line.strip().split("\t")
                if parts[0] == "Cluster":
                    continue
                else:
                    cluster = parts[0]
                    cluster = re.sub(r"^CL0*", r"", cluster)
                    annotation = parts[1]
                    custom_annotation[cluster] = annotation

        for cluster in clusters:
            cluster_tarean = cluster["Cluster"]
            annotation_tarean = cluster["Annotation"]
            cluster_tarean_clean = re.sub(r"^CL0*", r"", cluster_tarean)
            if cluster_tarean_clean in custom_annotation:
                cluster["Annotation"] = custom_annotation[cluster_tarean_clean]

    def sat_annotation(self, clusters, fasta_sats):
        fasta_seqs = SeqIO.parse(fasta_sats, "fasta")

        for cluster in clusters:
            cluster_seq = cluster["consensus"]
            for seq in fasta_seqs:
                if cluster_seq == seq.seq:
                    cluster["Annotation"] = f"{seq.id}"


    def plot_clusters(self):
        subprocess.run(["plot_clusters.R", f"{self.output_base}.tsv", f"plot{self.output_base}.png"])


def main():
    args = parse_args()
    html_file = args.html_file
    output_base = args.o
    g_size = args.g_size
    #sanity checking custom annotation files


    tarean = TAREAN_post_processing(html_file, output_base, g_size)
    tarean.correct_metrics(tarean.metrics, tarean.clusters)
    tarean.estimate_copy_number(tarean.clusters, tarean.metrics, g_size)
    #if custom library is chosen and provided
    if args.custom_anno and args.repeatM_file:
        repeatM_file = args.repeatM_file
        tarean.overwrited_annotation(repeatM_file, tarean.clusters)
    elif args.custom_anno and not args.repeatM_file:
        print("Provide a valid annotation file")

    if args.custom_sat and args.fasta_sats:
        fasta_sats = args.fasta_sats
        tarean.sat_annotation(tarean.clusters, fasta_sats)
    elif args.custom_sat and not args.fasta_sats:
        print("Provide a valid fasta file")

    tarean.write_files(tarean.clusters, tarean.metrics)

    if args.plot:
        tarean.plot_clusters()

main()
