import argparse
import os
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

def parse_kmer2(outdir, file_name ,result):

    #outpath= os.path.join(outdir,"bla_snp_transcript_aa_kmer_1.tsv")
    outpath= os.path.join(outdir,file_name)
    with open(outpath , "a+") as ofn:
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
def map_kmer_gen(input_path, mer_list):

    """
    run example:
    result = misc.map_kmer_gen("/home/yup/hackerthon/src/DamTheRiver/dummy_kmer.txt",mer_list)

    :param input_path:
    :param mer_list:
    :return:

    """
    mer_list.sort()
    result=[]
    with open(input_path,"r") as f:
        for line in f:
            field= line.split("\t")
            mer, gene_id, transcript_id, cds_id = field
            if mer in mer_list:
                result.append('\t'.join([gene_id, mer])+'\n')
    return result

def main():


    parser = argparse.ArgumentParser()
    parser.add_argument("--input_tsv", type = str, required = True, help = "4-field tsv with peptide sequence in 4th field")
    parser.add_argument("--outdir", type = str, required = True, help = "Path to output directory")
    args = parser.parse_args()

    input_tsv = args.input_tsv
    outdir = args.outdir

    kmer_list = get_kmer_2(input_tsv)
    parse_kmer2(outdir, kmer_list)
	
	
if __name__ == "__main__":

    main()
