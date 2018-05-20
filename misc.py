import argparse

def get_kmer(cid,seq,k):
    mer_add=[]
    mer=[]
    ln=len(seq)
    for n in range(ln-k):
        mer_add = seq[n:n+k]
        mer.append( mer_add )
    mer=list(set(mer))
    result = {"cid":cid, "mer":mer}
    return result

def parse_kmer(output_path, gene_id,transcript_id, result,mode = "a+"):
    lines=[]
    with open(output_path,mode) as ofn:
        cid = result["cid"]
        for mer in result["mer"]:
            line="\t".join([gene_id,transcript_id,cid,mer])
        lines.append(line+"\n")
        ofn.writelines(lines)

def parse_kmer2(outdir, result):

    
    with open(outdir + "/kmers_by_gene.tsv", "a+") as ofn:
        ofn.writelines(result)


def get_kmer_2(input_path, mode = "r" ):
    final=[]
    with open(input_path, "r") as f:
        first_line = f.readline()
        #print first_line
        for line in f:
            field= line.split("\t")
            gene_id,transcript_id,cds_id,seq = field
            #print(cds_id,seq)
            mers=get_kmer(cds_id,seq,8)
            #print(mers)
            for mer in  mers["mer"]:
                #print(mer)
                add_on= '\t'.join([mer, gene_id, transcript_id, cds_id ])
                final.append(add_on+"\n")
        return final

def main():


    parser = argparse.ArgumentParser()
    parser.add_argument("--input_tsv", type = str, required = True, help = "4-field tsv with peptide sequence in 4th field")
    parser.add_argument("--outdir", type = str, required = True, help = "Path to output directory")
    args = parser.parse_args()

    input_tsv = args.input_tsv
    outdir = args.outdir

    kmer_list = get_kmer_2(input_tsv)
    parse_kmer2(kmer_list, outdir)
	
	
if __name__ == "__main__":

    main()
