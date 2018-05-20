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

def parse_kmer(output_path, gene_id,transcript_id, result,mode="a+"):
    lines=[]
    with open(output_path,mode) as ofn:
        cid = result["cid"]
        for mer in result["mer"]:
            line="\t".join([gene_id,transcript_id,cid,mer])
        lines.append(line+"\n")
        ofn.writelines(lines)

