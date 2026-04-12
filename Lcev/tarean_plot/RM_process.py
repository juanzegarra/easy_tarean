from collections import defaultdict
import re

def parse_repeatmasker_output(file_path):
    file_index = file_path.replace(".tabular", "")
  #Record every cluster and its TE score and annotation
    cluster_annotation = defaultdict(list)
    print(f"executing file {file_path}")
    with open(file_path, "r") as f:
        line_counter = 0
      #skips header
        for line in f:
          line_counter += 1
          if line_counter >= 3:
            #handles mal-formatted tsv files
            parts = line.strip().replace("C\t", "").split("\t")
            score, query, te = parts[0], parts[4], parts[9]

            #keeps only the cluster numer, contig number is not needed
            cluster = re.sub(r"^(CL\d+).*", r"\1", query)

            #excludes simple repeats and low complexity regions
            if te not in {"Simple_repeat", "Low_complexity"}:
              cluster_annotation[cluster].append((score, te))

    with open(f"{file_index}.tsv", "w") as output_file:

      output_file.write("Cluster\tAnnotation\tScore\n")

      for cluster in cluster_annotation:
          cluster_sum = defaultdict(float)

          for score, te in cluster_annotation[cluster]:
              cluster_sum[te] += float(score)

          if cluster_sum:  # avoid error if empty
              max_key, max_value = max(cluster_sum.items(), key=lambda x: x[1])
              output_file.write(f"{cluster}\t{max_key}\t{max_value}\n")

