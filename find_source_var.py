import argparse
import re
import splice_lib as sl




def import_protein_dict(full_proteins):

	protein_dict = {}

	with open(full_proteins, 'r') as file:

		next(file)

		for line in file:

			entry = line.strip().split()

			gene = entry[0].strip()

			transcript = entry[1].strip()

			cds = entry[2].strip()

			seq = entry[3].strip()

			if cds not in protein_dict:

				protein_dict[cds] = {"gene": gene, "transcript": transcript, "seq": seq}

	return protein_dict


def check_peptides(peptides, protein_dict, outdir):

	peptide_pos_out = open(outdir + "/peptide_positions.bed", 'w')

	with open(peptides, 'r') as file:

		for line in file:

			entry = line.strip().split()

			peptide = entry[0].strip()
			gene = entry[1].strip()
			transcript = entry[2].strip()
			cds = entry[3].strip()

			cds_expanded = cds.split("_")
			chrom = cds_expanded[0]
			strand = cds_expanded[-1]
			cds_exons_flat = map(int, cds_expanded[1:-1])
			cds_exons = []

			for i in range(0, len(cds_exons_flat) - 1, 2):

				cds_exons.append([cds_exons_flat[i], cds_exons_flat[i+1]])


			if cds in protein_dict:

				if peptide in protein_dict[cds]["seq"]:

					for match in re.finditer(peptide, protein_dict[cds]["seq"]):

						position = match.span()

						start = sl.genome_to_transcript_coords(position[0]*3, strand, cds_exons, direction = "TG")
						end = sl.genome_to_transcript_coords(position[1]*3-1, strand, cds_exons, direction = "TG")

						if start < end:

							peptide_pos_out.write("\t".join([chrom, str(start), str(end), peptide + "_" + protein_dict[cds]["gene"], "1000", strand]) + "\n")

						else:

							peptide_pos_out.write("\t".join([chrom, str(end), str(start), peptide + "_" + protein_dict[cds]["gene"], "1000", strand]) + "\n")



	peptide_pos_out.close()








def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--peptides", type = str, help = "TSV containing peptide sequences of interest", required = True)
	parser.add_argument("--full_proteins", type = str, help = "TSV containing full protein sequences", required = True)
	parser.add_argument("--outdir", type = str, help = "Path to output directory", required = True)


	args = parser.parse_args()

	check_peptides(args.peptides, import_protein_dict(args.full_proteins), args.outdir)


if __name__ == "__main__":

	main()