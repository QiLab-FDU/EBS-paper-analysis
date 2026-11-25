import csv
from pathlib import Path

script_dir = Path(__file__).resolve().parent
cut_tag_dir = script_dir.parents[1]
gene_set_dir = cut_tag_dir / "peak_annotation" / "gene_set"
bed_ref = cut_tag_dir.parent / "tair12_genome" / "tair12_final.bed"

# Profile gene sets use geneId. The three-way motif set is a separate input.
csv_files = sorted([
    *gene_set_dir.glob("*_intersect_*.csv"),
    *gene_set_dir.glob("*_setdiff_*.csv"),
])

for csv_path in csv_files:
    output_bed = script_dir / f"{csv_path.stem}.bed"

    with open(csv_path, newline="") as f:
        specific_genes = {row["geneId"] for row in csv.DictReader(f) if row["geneId"]}

    output_lines = []
    with open(bed_ref, 'r') as f_bed:
        for line in f_bed:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            # The fourth BED column contains the TAIR gene ID.
            gene_id = fields[3]
            if gene_id in specific_genes:
                output_lines.append(line)

    with open(output_bed, 'w') as f_out:
        f_out.writelines(output_lines)
    print(f"Wrote {output_bed}")

