
import copy
import re
from operator import itemgetter
import gzip
import sys
import translate_lib as tl

def get_junctions(exon_list):


    '''
        Takes a list of exons (list of lists with each element of the inner list containing the start, stop position of an exon) and returns a list of the junction positions (essentially the inner coordinates of any pair of exons).  This facilitates comparison of objects that may differ in start/stop positions but contain identical splice junctions.
    '''
    
    junction_list = []
    for i in range(0,len(exon_list) - 1):
        junction_list.append([exon_list[i][1],exon_list[i+1][0]])
    return junction_list


def junction_overlap(subset_jl, superset_jl):
    ##TODO replace with subset operation
    #print "I am the junction_overlap function!"

    '''
        Takes as input two lists of junctions - expected subset and expected superset.  Checks whether subset is truly a subset of superset.  Returns boolean.
    '''

    subset_subset = [i for i in subset_jl if i in superset_jl]

    if subset_subset == subset_jl and len(subset_subset) > 0:
        #print "Subset subset contains all of the junctions - checking for connectivity . . . "

        first_subset_junction_index = superset_jl.index(subset_subset[0])
        last_subset_junction_index = superset_jl.index(subset_subset[-1])
        superset_junction_subset = superset_jl[first_subset_junction_index:last_subset_junction_index+1]
    
        if superset_junction_subset == subset_subset:
            #print "Exact junction match!"
            return True
        else:
            #print "Connectivity interrupted - no match"
            return False
    else:
        #print "subset_subset missing junctions - no match"
        return False


def position_contained(exon_list, position):

    '''
        Takes as input list of exon coordinates (list of lists, each element is a list containing start, stop position of exon), position (single genomic coordinate)

        Checks whether position is contained within any of the exons, returns list containing boolean (True if the position was contained in one of the exons), the containing exon entry (or None of the position was not contain)
    '''

    for exon in exon_list:


        if ((int(position) >= int(exon[0])) and (int(position) <= int(exon[1]))):

            containing_exon = exon
            return [True, containing_exon]

    else:

        return [False, None]


def translate_ORF(transcript_seq, stop_codons, valid_adj_cds_start):
    '''
        Takes as input transcript sequence, list of stop codons (e.g. ["TAA", "TAG", "TGA"]), start codon.  Attempts to find in-frame stop-codon by translating the sequence from the start codon.

        Returns a list containing boolean (True if in-frame stop codon is found), stop codon coordinate (None if nothing is found) in transcript-centric coords

    '''

    CDS_seq_left_bound = transcript_seq[valid_adj_cds_start:].upper()

    if CDS_seq_left_bound[0:3] != "ATG":

        return [False, "Not a real start", CDS_seq_left_bound[0:5]]

    codon_list = []

    for i in range(0,len(CDS_seq_left_bound) - 2,3):

        codon = CDS_seq_left_bound[i:i+3]
        if codon not in stop_codons:
            codon_list += [codon]
        else:
            #print "In-frame stop codon found!"
            #codon_list += [codon]
            #print "".join(codon_list)
            adj_cds_end = valid_adj_cds_start + 3*len(codon_list) - 1 ###Again check this
            codon_list=("").join(codon_list)
            result= [True,adj_cds_end,codon_list]

    else:
        #print "No in-frame stop codon found"
        result=[False, None,codon_list]

    AA= tl.translate_seq(codon_list)
    result.append(AA)
    return result

def genome_to_transcript_coords(position, strand, transcript_exons, direction = "TG"): ##Where "TG" is transcript -> genome and "GT" is genome -> transcript
    '''
        Interchanges transcript coords (distance along the spliced transcript with reference to the 5'-end of the transcript) with genomic coordinates
    '''

    position = int(position)
    
    if strand == "+":
        gtc_transcript_exons = copy.deepcopy(transcript_exons)
    elif strand == "-":
        gtc_transcript_exons = copy.deepcopy(transcript_exons)[::-1]

    exon_lengths = [0] ## zero is just to prevent list index out of range on first iteration of loop below (cumulative lengths)

    individual_exon_lengths = []

    for i in gtc_transcript_exons:

        prev_length = exon_lengths[-1]
        exon_lengths += [prev_length + (int(i[1]) - int(i[0]) + 1)]
        individual_exon_lengths += [(int(i[1]) - int(i[0]) + 1)]

    #print individual_exon_lengths


    if direction == "TG":    
        
        for i in range(0,len(exon_lengths)-1):
            if position < exon_lengths[i+1]:
                if strand == "+":
                    return int(gtc_transcript_exons[i][0]) + position - exon_lengths[i]  ###double check this!!!  ####MAKE IT TRIPLE!!!!!!!!!!!!!!!!!!!
                elif strand == "-":
                    return int(gtc_transcript_exons[i][1]) -  position + exon_lengths[i] ###double check this!!!!

        else:
            print "GARBAGE - position not found. genome_to_transcript_coords"    
            print position
            print strand
            print transcript_exons
            print direction
            print exon_lengths
            sys.exit()
                
    elif direction == "GT":

        for i in gtc_transcript_exons:
            if position in range(int(i[0]),int(i[1]) + 1):
                if strand == "+":
                    return abs(position - int(i[0])) + int(exon_lengths[gtc_transcript_exons.index(i)])
                elif strand == "-":
                    return abs(position - int(i[1])) + int(exon_lengths[gtc_transcript_exons.index(i)])

        else:
            print "GARBAGE - position not found. genome_to_transcript_coords"    
            print position
            print strand
            print transcript_exons
            print direction
            print exon_lengths
            sys.exit()

#####QUICK EXAMPLE
#
#    genomic position: 236
#    start of exon 1: 231 (meaning that the 1st base in the exon is 231)
#    stop of exon 1: 270
#    start of CDS: 236
#    exon_lengths[transcript_exons.index([231,270])] is equal to exon_lengths[0], which is always going to be 0 the way it's written above
#    therefore, we return 236 - 231 + 0 = 5
#    However, the correct 1-base answer would be 6
#    We can fix this by supplanting 0 with 1 in the exon_lengths list after all the lenglth calculations are finished### NEVERMIND - just add 1 going GT and subtract 1 going TG
#
#    Going backards:
#    transcript coordinate = 6 + 231 - 1 = 236
#
#    More complicated: exon1: [231,270]; exon 2: [321,370]
#        exon_lengths = [0, 40, 90]
#        exon_lengths = [1, 40, 90] # after swapping out the 0 #####AGAIN DONT BOTHER WITH THIS - it's fixed by adding/subtracting 1 as below
#        genomic position = 326 ---> transcript position = 46
#        GT calculation: 326 - 321 + 40 + 1 = 5 + 40 + 1 = 46
#        TG calculation: 321 + 46 - 40 - 1 = 321 + 6 - 1 = 326
#
#    OK now what about example on the minus strand:
#        same exon positions (but the order is now swapped)    
#        exon_lengths = [0, 50, 90]
#        transcript position = 46
#        genomic position = 325
#        GT calculation: 370 - 325 + (length term which is 0 since it's in the first exon) + 1
#        TG calculation: 370 - 46 + (length term - again 0 here)
#
#
#

def get_feature_start_end(strand, feature_exons):
    '''
        Extracts the start, end positions of features based on the exons of the feature.  Start is the left-most position of the leftmost exon for features on the "+" strand, and 
        the right-most position of the right-most exon for features on the "-" strand.  End is the reverse for both cases.
    '''

    if strand == "+":
        return [feature_exons[0][0],feature_exons[-1][1]]

    elif strand == "-":
        return [feature_exons[-1][1], feature_exons[0][0]]


def get_cds_utr_exons(strand, feature_exons, cds_start, cds_end, start_container, end_container):
    '''

    '''


    feature_exons = copy.deepcopy(feature_exons)

    if strand == "+":

        five_utr_subset = copy.deepcopy(feature_exons[0:feature_exons.index(start_container) + 1])
        #try:
        cds_subset = copy.deepcopy(feature_exons[feature_exons.index(start_container):feature_exons.index(end_container) + 1])
        #except ValueError:
        #    print "Disaster afoot in get_cds_utr_exons"
        #    print strand
        #    print feature_exons
        #    print cds_start
        #    print cds_end
        #    print start_container
        #    print end_container
        #    sys.exit()
        three_utr_subset = copy.deepcopy(feature_exons[feature_exons.index(end_container):])

        five_utr_subset[-1][-1] = int(cds_start) - 1

        cds_subset[0][0] = cds_start
        cds_subset[-1][-1] = cds_end

        three_utr_subset[0][0] = int(cds_end) + 1

    if strand == "-":

        three_utr_subset = copy.deepcopy(feature_exons[0:feature_exons.index(end_container) + 1])

        cds_subset = copy.deepcopy(feature_exons[feature_exons.index(end_container):feature_exons.index(start_container) + 1])
        five_utr_subset = copy.deepcopy(feature_exons[feature_exons.index(start_container):])


        five_utr_subset[0][0] = int(cds_start) + 1

        cds_subset[0][0] = cds_end
        cds_subset[-1][-1] = cds_start

        three_utr_subset[-1][-1] = int(cds_end) - 1

    return [five_utr_subset, cds_subset, three_utr_subset]

def get_exon_subset(feature_exons, genomic_start, genomic_end):

    if not genomic_start < genomic_end:
        raise ValueError("genomic end not larger than genomic start in call to get_exon_subset")

    start_found = False
    end_found = False

    for exon in feature_exons:

        if not start_found and genomic_start in range(exon[0], exon[1] + 1):

            start_exon = exon
            start_found = True

        if not end_found and genomic_end in range(exon[0], exon[1] + 1):

            end_exon = exon
            end_found = True

        if start_found and end_found:

            break

    else:
        print "Either start or end not found in provided exon list"
        print feature_exons
        print genomic_start
        print genomic_end
        sys.exit()




    exon_subset = copy.deepcopy(feature_exons[feature_exons.index(start_exon):feature_exons.index(end_exon) + 1])

    exon_subset[0][0] = genomic_start
    exon_subset[-1][-1] = genomic_end

    return exon_subset


def calc_length_exon_list(exon_list):
    '''Takes a given list of exons (e.g. [[1000,1100],[2000,2100]]) where each sublist contains the start/stop position of each exon (1-based) and finds the total exonic length of the feature'''

    length = 0

    exon_list = [map(int, i) for i in exon_list]

    for i in exon_list:

        length += i[1] - i[0] + 1

    return length



def generate_standard_event_dict(event_gtf_filename):
    '''
        Generates a dictionary of events indexed by event ID
    '''

    event_types = ["A3", "A5", "AF", "AL", "RI", "MX", "SE", "MS", "MF", "ML", "CF", "CL", "CO", "AT", "AP"]

    standard_event_dict = {}

    if event_gtf_filename.endswith(".gz"):

        gtf_file = gzip.open(event_gtf_filename, 'rb')

    else:

        gtf_file = open(event_gtf_filename, 'r')

    for line in gtf_file:

        entry = line.split()

        if entry[2] == "exon":

            start = int(entry[3])
            end = int(entry[4])

            chrom = re.sub("_","&",entry[0].strip())
            strand = entry[6].strip()

            transcript_id_entry = re.findall('transcript_id\s\"[^;\"]+\";', line)

            if len(transcript_id_entry) == 0:

                sys.exit("Bad gtf file - 'transcript_id' not found in exon entry")

            elif len(transcript_id_entry) > 1:

                print "Bad transcript_id - setting transcript_id to 'bad_transcript_id'"

                transcript_id = "bad_transcript_id"

            else:

                    transcript_id = re.sub('[";]', '', transcript_id_entry[0].strip().split()[1])    


            event_id = "_".join(transcript_id.split("_")[0:-1])
            form = transcript_id.split("_")[-1]

            for i in event_types:

                if i in event_id:

                    event_type = i
                    break

            if event_id not in standard_event_dict:

                standard_event_dict[event_id] = {
                    "event_type": event_type,
                    "included_exons": [],
                    "excluded_exons": [],
                    "chrom": chrom,
                    "strand": strand,
                    "included_count": 0,
                    "excluded_count": 0,
                    "included_form_transcripts": [],
                    "excluded_form_transcripts": [],
                    "sources": []
                }

            if form == "alternative1" or form == "iso2" or form == "included":

                standard_event_dict[event_id]["included_exons"].append([start, end])

            if form == "alternative2" or form == "iso1" or form == "excluded":

                standard_event_dict[event_id]["excluded_exons"].append([start, end])

    gtf_file.close()

    return standard_event_dict


def find_exons_unique_to_form(standard_event_dict):
    '''
        post-hoc method to find exons unique to included, excluded forms respectively and add them to new fields in event_dict
    '''

    for event in standard_event_dict:


        included_set = set(["_".join(map(str, i)) for i in standard_event_dict[event]["included_exons"]])
        excluded_set = set(["_".join(map(str, i)) for i in standard_event_dict[event]["excluded_exons"]])


        included_unique = included_set - excluded_set
        excluded_unique = excluded_set - included_set

        standard_event_dict[event]["included_unique_exons"] = [map(int, i.split("_")) for i in included_unique]
        standard_event_dict[event]["excluded_unique_exons"] = [map(int, i.split("_")) for i in excluded_unique]


def assess_event_types(standard_event_dict):

    event_type_counts = {
    "SE": 0,
    "A3": 0,
    "A5": 0,
    "MS": 0,
    "MX": 0,
    "RI": 0,
    "AF": 0,
    "AL": 0,
    "CO": 0,
    "MF": 0,
    "ML": 0,
    "CF": 0,
    "CL": 0,
    "UF": 0,
    "UL": 0,
    "AT": 0,
    "AP": 0,
    "AB": 0,
    "total": 0
    }
    

    for i in standard_event_dict:

        event_type_counts[standard_event_dict[i]["event_type"]] += 1
        event_type_counts["total"] += 1

    return event_type_counts


def complete_event_dict(standard_event_dict):

    '''
        Sorts event exons, infers junctions - necessary to complete event dict

        Other functions also have this effect - e.g. collapse_redundant_junctions.
    '''

    for event in standard_event_dict:

        included_exons = [map(str, i) for i in sorted(standard_event_dict[event]["included_exons"], key=itemgetter(0))]
        standard_event_dict[event]["included_exons"] = copy.deepcopy(included_exons)

        excluded_exons = [map(str,i) for i in sorted(standard_event_dict[event]["excluded_exons"], key=itemgetter(0))]
        standard_event_dict[event]["excluded_exons"] = copy.deepcopy(excluded_exons)

        included_jl = get_junctions(included_exons)
        excluded_jl = get_junctions(excluded_exons)

        standard_event_dict[event]["included_junctions"] = copy.deepcopy(included_jl)
        standard_event_dict[event]["excluded_junctions"] = copy.deepcopy(excluded_jl)

        standard_event_dict[event]["included_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in included_jl}
        standard_event_dict[event]["excluded_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in excluded_jl}

        ###handle RI included "junctions"

        if standard_event_dict[event]["event_type"] == "RI":  

            chrom = standard_event_dict[event]["chrom"]
            strand = standard_event_dict[event]["strand"]
            start = str(int(standard_event_dict[event]["excluded_junctions"][0][0]) + 1)
            end = str(int(standard_event_dict[event]["excluded_junctions"][0][1]) - 1)

            standard_event_dict[event]["included_junctions"] = included_jl = [[int(start), int(start) + 1], [int(end) - 1, int(end)]]
            standard_event_dict[event]["included_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in included_jl}

        standard_event_dict[event]["included_jn_count"] = str(len(included_jl)) ##confusing terminology! this is just the number of junctions
        standard_event_dict[event]["excluded_jn_count"] = str(len(excluded_jl))



def check_integrity(event_dict):

    incomplete_events = []
    
    for event in event_dict:

        if event_dict[event]["event_type"] == "SE":

            if len(event_dict[event]["included_exons"]) != 3 or len(event_dict[event]["excluded_exons"]) != 2:

                incomplete_events.append(event)

                print event_dict[event]

        elif event_dict[event]["event_type"] == "MS":

            if len(event_dict[event]["included_exons"]) <= 3 or len(event_dict[event]["excluded_exons"]) != 2:

                incomplete_events.append(event)            

        elif event_dict[event]["event_type"] in ["A3", "A5", "AL", "AF"]:

            if len(event_dict[event]["included_exons"]) != 2 or len(event_dict[event]["excluded_exons"]) != 2:

                incomplete_events.append(event)


        elif event_dict[event]["event_type"] == "RI":

            if len(event_dict[event]["included_exons"]) != 1 or len(event_dict[event]["excluded_exons"]) != 2:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] == "MX":

            if len(event_dict[event]["included_exons"]) != 3 or len(event_dict[event]["excluded_exons"]) != 3:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] in ["AT", "AP"]:

            if len(event_dict[event]["included_exons"]) != 1 or len(event_dict[event]["excluded_exons"]) != 1:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] in ["MF", "ML"]:

            if (len(event_dict[event]["included_exons"]) < 3 and len(event_dict[event]["excluded_exons"]) < 3) or len(event_dict[event]["included_exons"]) < 2 or len(event_dict[event]["excluded_exons"]) < 2:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] in ["UF", "UL"]:

            if len(event_dict[event]["excluded_exons"]) != 1 or len(event_dict[event]["included_exons"]) < 2:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] in ["CF", "CL"]:

            if (len(event_dict[event]["included_exons"]) < 2 and len(event_dict[event]["excluded_exons"]) < 2) or len(event_dict[event]["included_exons"]) < 1 or len(event_dict[event]["excluded_exons"]) < 1:

                incomplete_events.append(event)

        elif event_dict[event]["event_type"] == "CO":

            if (len(event_dict[event]["included_exons"]) < 2 and len(event_dict[event]["excluded_exons"]) < 2) or len(event_dict[event]["included_exons"]) < 1 or len(event_dict[event]["excluded_exons"]) < 1:

                incomplete_events.append(event)


    for event in incomplete_events:

        del event_dict[event]



def generate_standard_transcript_dict(transcript_gtf_filename, ccds = False, feature = "exon"):
    '''
        Generates a dictionary of transcripts indexed by transcript ID
        Presently limited to stringtie formatting
    '''

    standard_transcript_dict = {}

    if transcript_gtf_filename.endswith(".gz"):

        gtf_file = gzip.open(transcript_gtf_filename, 'rb')

    else:

        gtf_file = open(transcript_gtf_filename, 'r')


    for line in gtf_file:

        if not line.startswith("#"):

            if ccds:

                if "CCDS" not in line:

                    continue

            entry = line.split()

            if entry[2] == feature:

                start = int(entry[3])
                end = int(entry[4])

                chrom = re.sub("_","&",entry[0].strip())
                strand = entry[6].strip()

                transcript_id_entry = re.findall('transcript_id\s\"[^;\"]+\";', line)

                if len(transcript_id_entry) == 0:

                    sys.exit("Bad gtf file - 'transcript_id' not found in exon entry")

                elif len(transcript_id_entry) > 1:

                    print "Bad transcript_id - setting transcript_id to 'bad_transcript_id'"

                    transcript_id = "bad_transcript_id"

                else:

                    transcript_id = re.sub('[";]', '', transcript_id_entry[0].strip().split()[1])

                gene_name_entry = re.findall('gene_name\s\"[^;\"]+\";', line)

                if len(gene_name_entry) == 0:

                    gene_id_entry = re.findall('gene_id\s\"[^;\"]+\";', line)

                    if len(gene_id_entry) == 0:

                        gene = ""

                    elif len(gene_id_entry) > 1:

                        print "Bad gene name - setting gene name to 'bad_gene_name'"

                        gene = "bad_gene_name"

                    else:

                        gene = re.sub('[";]', '', gene_id_entry[0].strip().split()[1])

                elif len(gene_name_entry) > 1:

                    print "Bad gene name - setting gene name to 'bad_gene_name'"

                    gene = "bad_gene_name"

                else:

                    gene = re.sub('[";]', '', gene_name_entry[0].strip().split()[1])


                if transcript_id not in standard_transcript_dict:

                    standard_transcript_dict[transcript_id] = {
                        "exons": [],
                        "chrom": chrom,
                        "strand": strand,
                        "CDS": {},
                        "nonstop": {},
                        "included_form_events": [],
                        "excluded_form_events": [],
                        "gene": gene
                    }

                standard_transcript_dict[transcript_id]["exons"].append([start, end])

    gtf_file.close()

    return standard_transcript_dict



def sort_transcript_dict_exons(transcript_dict):

    for transcript in transcript_dict:

        transcript_dict[transcript]["exons"] = sorted(transcript_dict[transcript]["exons"], key=itemgetter(0))



def add_junctions_to_transcript_dict(transcript_dict):

    for transcript in transcript_dict:

        transcript_dict[transcript]["junctions"] = get_junctions(transcript_dict[transcript]["exons"])

        transcript_dict[transcript]["junction_counts"] = dict(zip([transcript_dict[transcript]["chrom"] + "_" + "_".join(map(str, junction)) + "_" + transcript_dict[transcript]["strand"] for junction in transcript_dict[transcript]["junctions"]], [0]*len(transcript_dict[transcript]["junctions"])))

        transcript_dict[transcript]["total_junction_count"] = 0



def index_transcripts_by_junctions(full_transcript_dict):
    '''
        Intended to take transcript dict supplied by generate_standard_transcript_dict method.

        Returns dictionary where splice junctions (specified by chrom, strand, and inner coordinates of flanking exons) are keys which point to lists of matching transcripts

        Takes the transcript dictionary as input.

        This facilitates the matching of CDS without needing prior gene information.

    '''

    junction_dict = {}

    for transcript in full_transcript_dict:

        transcript_jl = copy.deepcopy(full_transcript_dict[transcript]["junctions"])

        for junction in transcript_jl:

            junction_key = full_transcript_dict[transcript]["chrom"] + "_" + "_".join(map(str, junction)) + "_" + full_transcript_dict[transcript]["strand"]

            if junction_key not in junction_dict:

                junction_dict[junction_key] = [transcript]

            else:
                junction_dict[junction_key].append(transcript)

    return junction_dict



def index_transcripts_by_donor_acceptor(full_transcript_dict):
    '''
        Returns list of three dictionaries: first indexes transcripts by each individual donor/acceptor coordinate (without identifying specifically as a donor or acceptor)
        second two index transcripts by specifically either the first acceptor, or the last donor.  The latter two dictionaries allows inference of pairwise AF/AL events. 
        The first of the latter two indexes by the second from the leftmost half junction, while the second indexes with the second from the last half junction

    '''

    donor_acceptor_dict = {}
    first_acceptor_last_donor_dict_left = {} ##AF if strand +, AL if strand -
    first_acceptor_last_donor_dict_right = {} ##AL if strand +, AF if strand -

    for transcript in full_transcript_dict:

        transcript_exons = sorted(copy.deepcopy(full_transcript_dict[transcript]["exons"]), key=itemgetter(0))
        full_transcript_dict[transcript]["exons"] = transcript_exons

        transcript_jl = get_junctions(transcript_exons)
        full_transcript_dict[transcript]["junctions"] = copy.deepcopy(transcript_jl)
        flat_jl = [i for j in transcript_jl for i in j]
        full_transcript_dict[transcript]["flat_junctions"] = copy.deepcopy(flat_jl)

        #print transcript_jl

        for index, half_junction in enumerate(flat_jl):

            half_junction_key = full_transcript_dict[transcript]["chrom"] + "_" + str(half_junction) + "_" + full_transcript_dict[transcript]["strand"]

            if half_junction_key not in donor_acceptor_dict:

                donor_acceptor_dict[half_junction_key] = []

            if transcript not in donor_acceptor_dict[half_junction_key]:

                donor_acceptor_dict[half_junction_key].append(transcript)

            if index == 1:

                if half_junction_key not in first_acceptor_last_donor_dict_left:

                    first_acceptor_last_donor_dict_left[half_junction_key] = []

                if transcript not in first_acceptor_last_donor_dict_left[half_junction_key]:

                    first_acceptor_last_donor_dict_left[half_junction_key].append(transcript)

            elif index == (len(flat_jl) - 2):

                if half_junction_key not in first_acceptor_last_donor_dict_right:

                    first_acceptor_last_donor_dict_right[half_junction_key] = []

                if transcript not in first_acceptor_last_donor_dict_right[half_junction_key]:

                    first_acceptor_last_donor_dict_right[half_junction_key].append(transcript)

    return [donor_acceptor_dict, first_acceptor_last_donor_dict_left, first_acceptor_last_donor_dict_right]


def collapse_redundant_junction_events(standard_event_dict, outdir):


    '''
        Searches for events with identical included+excluded form junction sets (happens a lot with AF events (probably AL too)), collapses to a single event
    '''

    event_type_clash = open(outdir + "/event_type_clash.tsv", 'w')

    junction_set_dict = {}

    for event in standard_event_dict:

        included_exons = [map(str, i) for i in sorted(standard_event_dict[event]["included_exons"], key=itemgetter(0))]
        standard_event_dict[event]["included_exons"] = copy.deepcopy(included_exons)

        excluded_exons = [map(str,i) for i in sorted(standard_event_dict[event]["excluded_exons"], key=itemgetter(0))]
        standard_event_dict[event]["excluded_exons"] = copy.deepcopy(excluded_exons)

        included_jl = get_junctions(included_exons)
        excluded_jl = get_junctions(excluded_exons)

        standard_event_dict[event]["included_junctions"] = copy.deepcopy(included_jl)
        standard_event_dict[event]["excluded_junctions"] = copy.deepcopy(excluded_jl)

        standard_event_dict[event]["included_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in included_jl} #this is the number of reads aligning to these junctions
        standard_event_dict[event]["excluded_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in excluded_jl}

        if standard_event_dict[event]["event_type"] == "RI":  

            chrom = standard_event_dict[event]["chrom"]
            strand = standard_event_dict[event]["strand"]
            start = str(int(standard_event_dict[event]["excluded_junctions"][0][0]) + 1)
            end = str(int(standard_event_dict[event]["excluded_junctions"][0][1]) - 1)

            standard_event_dict[event]["included_junctions"] = included_jl = [[int(start), int(start) + 1], [int(end) - 1, int(end)]]
            standard_event_dict[event]["included_junction_counts"] = {(standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,i)) + "_" + standard_event_dict[event]["strand"]):0 for i in included_jl}

        standard_event_dict[event]["included_jn_count"] = str(len(included_jl)) ##confusing terminology! this is just the number of junctions
        standard_event_dict[event]["excluded_jn_count"] = str(len(excluded_jl))


        flat_joined_included_jl = "_".join([str(i) for j in included_jl for i in j])
        flat_joined_excluded_jl = "_".join([str(i) for j in excluded_jl for i in j])

        chrom = standard_event_dict[event]["chrom"]
        strand = standard_event_dict[event]["strand"]

        if flat_joined_included_jl != "" and flat_joined_excluded_jl != "":

            key = chrom + "_" + flat_joined_included_jl + "_" +  flat_joined_excluded_jl + "_" + strand

            if key not in junction_set_dict:

                junction_set_dict[key] = []

            junction_set_dict[key].append(event)

    for key in junction_set_dict:

        if len(junction_set_dict[key]) > 1:

            new_key = ",".join(junction_set_dict[key])

            event_types = []

            for event in junction_set_dict[key]:

                if standard_event_dict[event]["event_type"] not in event_types:

                    event_types.append(standard_event_dict[event]["event_type"])

            if not all([i == event_types[0] for i in event_types]):

                event_type_clash.write("_".join(junction_set_dict[key]) + "\t" + key + "\n")
                standard_event_dict[new_key] = copy.deepcopy(standard_event_dict[junction_set_dict[key][1]])
                standard_event_dict[new_key]["event_type"] = "AB"

            else:

                standard_event_dict[new_key] = copy.deepcopy(standard_event_dict[junction_set_dict[key][1]])

            for old_key in junction_set_dict[key]:

                del standard_event_dict[old_key]

    event_type_clash.close()


def generate_junction_indexed_event_dict(standard_event_dict):

    '''
        Returns dictionary where splice junctions (specified by chrom, strand, and inner coordinates of flanking exons) are keys which point to lists of matching transcripts

        Takes the transcript dictionary as input.

        This facilitates the matching of CDS without needing prior gene information.

    '''

    junction_dict = {}

    for event in standard_event_dict:

        included_jl = standard_event_dict[event]["included_junctions"] 
        excluded_jl = standard_event_dict[event]["excluded_junctions"] 

        for junction in included_jl:

            junction_key = standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,junction)) + "_" + standard_event_dict[event]["strand"]

            if junction_key not in junction_dict:

                junction_dict[junction_key] = [event + "_" + "included"]


            else:
                junction_dict[junction_key].append(event + "_" + "included")

        for junction in excluded_jl:

            junction_key = standard_event_dict[event]["chrom"] + "_" + "_".join(map(str,junction)) + "_" + standard_event_dict[event]["strand"]

            if junction_key not in junction_dict:

                junction_dict[junction_key] = [event + "_" + "excluded"]

            else:
                junction_dict[junction_key].append(event + "_" + "excluded")

    return junction_dict



def filter_overlapping_se_alt_donacc_events(standard_event_dict, junction_dict):

    events_to_filter = []
    for junction in junction_dict:

        if len(junction_dict[junction]) > 1:

            events = []
            event_types = []
            eventsoi = [] #events of interest - filtered for events that are SE, A5, or A3\
            event_typesoi = []

            for isoform in junction_dict[junction]:

                event = isoform.split("_")[0]

                events.append(event)

                event_types.append(standard_event_dict[event]["event_type"])

            if "SE" in event_types and ("A3" in event_types or "A5" in event_types):

                for i, item in enumerate(event_types):

                    if item in ["SE","A5","A3"]:

                        eventsoi.append(events[i])
                        event_typesoi.append(item)

                for i in range(0, len(eventsoi)):

                    for j in range(i+1, len(eventsoi)):

                        if event_typesoi[i] == "SE" and (event_typesoi[j] == "A3" or event_typesoi[j] == "A5"):

                            if standard_event_dict[eventsoi[i]]["excluded_junctions"] == standard_event_dict[eventsoi[j]]["excluded_junctions"]:

                                if standard_event_dict[eventsoi[j]]["included_junctions"][0] in standard_event_dict[eventsoi[i]]["included_junctions"]:


                                    if eventsoi[j] not in events_to_filter:
                                        #print "Filtering", eventsoi[j]
                                        events_to_filter.append(eventsoi[j])


                        elif event_typesoi[j] == "SE" and (event_typesoi[i] == "A3" or event_typesoi[i] == "A5"):

                            if standard_event_dict[eventsoi[i]]["excluded_junctions"] == standard_event_dict[eventsoi[j]]["excluded_junctions"]:

                                if standard_event_dict[eventsoi[i]]["included_junctions"][0] in standard_event_dict[eventsoi[j]]["included_junctions"]:

                                    if eventsoi[i] not in events_to_filter:

                                        #print "Filtering", eventsoi[i]
                                        events_to_filter.append(eventsoi[i])

    for event in events_to_filter:

        del standard_event_dict[event]


def rename_events(standard_event_dict):

    counting_dict = {
    "SE": 1,
    "A3": 1,
    "A5": 1,
    "MS": 1,
    "MX": 1,
    "RI": 1,
    "AF": 1,
    "AL": 1,
    "CO": 1,
    "MF": 1,
    "ML": 1,
    "CF": 1,
    "CL": 1,
    "UF": 1,
    "UL": 1,
    "AT": 1,
    "AP": 1,
    "AB": 1
    }
    
    new_event_dict = {}

    for event in standard_event_dict.keys():

        if standard_event_dict[event]["event_type"] in counting_dict:

            new_event_name = standard_event_dict[event]["event_type"] + "." + str(counting_dict[standard_event_dict[event]["event_type"]]).zfill(7)

            counting_dict[standard_event_dict[event]["event_type"]] += 1

            new_event_dict[new_event_name] = copy.deepcopy(standard_event_dict[event])
        else: 
            print "WHOOPS"

    return new_event_dict


def output_event_gtf(standard_event_dict, outdir, name="splice_lib_events"):

    with open(outdir + "/" + name + ".gtf", 'w') as file:

        for event in standard_event_dict:

            chrom = re.sub("&","_",standard_event_dict[event]["chrom"])
            strand = standard_event_dict[event]["strand"]

            for exon in standard_event_dict[event]["included_exons"]:

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_included"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + isoform_id + '";']) + "\n")

            for exon in standard_event_dict[event]["excluded_exons"]:

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_excluded"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + isoform_id + '";']) + "\n")

def output_event_bedfile(standard_event_dict, outdir, name = "splice_lib_events"):

    with open(outdir + "/" + name + ".bed", 'w') as file:

        for event in standard_event_dict:

            chrom = re.sub("&","_",standard_event_dict[event]["chrom"])
            strand = standard_event_dict[event]["strand"]    

            start = str(standard_event_dict[event]["included_exons"][0][0])
            end = str(standard_event_dict[event]["included_exons"][-1][1])

            file.write("\t".join([chrom, start, end, event, "1000", strand]) + "\n")    


def output_event_gtf_with_transcripts(standard_event_dict, outdir, name="splice_lib_events_with_transcripts"):

    with open(outdir + "/" + name + ".gtf", 'w') as file:

        for event in standard_event_dict:

            chrom = re.sub("&","_",standard_event_dict[event]["chrom"])
            strand = standard_event_dict[event]["strand"]

            gene_transcript_start = str(standard_event_dict[event]["included_exons"][0][0])
            gene_transcript_end = str(standard_event_dict[event]["included_exons"][-1][1])


            file.write("\t".join([chrom, "splice_lib_event", "transcript", gene_transcript_start, gene_transcript_end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + event + "_included" + '";']) + "\n")

            for exon in standard_event_dict[event]["included_exons"]:

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_included"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + isoform_id + '";']) + "\n")

            file.write("\t".join([chrom, "splice_lib_event", "transcript", gene_transcript_start, gene_transcript_end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + event + "_excluded" + '";']) + "\n")

            for exon in standard_event_dict[event]["excluded_exons"]:

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_excluded"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "gene_id " + '"' + event + '"; transcript_id ' + '"' + isoform_id + '";']) + "\n")


def output_miso_event_gff3(standard_event_dict, outdir, name="splice_lib_events"):

    with open(outdir + "/" + name + ".gff3", 'w') as file:

        file.write("##gff-version 3" + "\n")

        for event in standard_event_dict:


            chrom = re.sub("&","_",standard_event_dict[event]["chrom"])
            strand = standard_event_dict[event]["strand"]


            gene_transcript_start = str(standard_event_dict[event]["included_exons"][0][0])
            gene_transcript_end = str(standard_event_dict[event]["included_exons"][-1][1])

            gene = parent = event

            file.write("\t".join([chrom, "splice_lib_event", "gene", gene_transcript_start, gene_transcript_end, ".", strand, ".", "ID=" + event + ";Name=" + event]) + "\n")
            file.write("\t".join([chrom, "splice_lib_event", "transcript", gene_transcript_start, gene_transcript_end, ".", strand, ".", "ID=" + event + "_included" + ";Parent=" + event]) + "\n")


            for index, exon in enumerate(standard_event_dict[event]["included_exons"]):

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_included"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "ID=exon:" + event + "_included:" + str(index + 1) + ";Parent=" + event + "_included"]) + "\n")


            file.write("\t".join([chrom, "splice_lib_event", "transcript", gene_transcript_start, gene_transcript_end, ".", strand, ".", "ID=" + event + "_excluded" + ";Parent=" + event]) + "\n")


            for index, exon in enumerate(standard_event_dict[event]["excluded_exons"]):

                start = str(exon[0])
                end = str(exon[1])
                isoform_id = event + "_excluded"

                file.write("\t".join([chrom, "splice_lib_event", "exon", start, end, ".", strand, ".", "ID=exon:" + event + "_excluded:" + str(index + 1) + ";Parent=" + event + "_excluded"]) + "\n")



def output_transcript_gtf(standard_transcript_dict, outdir, name = "splice_lib_transcripts"):

    with open(outdir + "/" + name + ".gtf", 'w') as file:

        for transcript in standard_transcript_dict:

            chrom = re.sub("&", "_", standard_transcript_dict[transcript]["chrom"])
            strand = standard_transcript_dict[transcript]["strand"]

            start = str(standard_transcript_dict[transcript]["exons"][0][0])
            end = str(standard_transcript_dict[transcript]["exons"][-1][1])

            try:
                gene = standard_transcript_dict[transcript]["gene"]
            except KeyError:
                sys.exit("Tried to output a transcript gtf with no gene name. Exiting . . . ")

            file.write("\t".join([chrom, "splice_lib_transcript", "transcript", start, end, ".", strand, ".", "gene_name " + '"' + gene + '"; transcript_id ' + '"' + transcript + '";']) + "\n")

            for exon in standard_transcript_dict[transcript]["exons"]:

                start = str(exon[0])
                end = str(exon[1])

                file.write("\t".join([chrom, "splice_lib_transcript", "exon", start, end, ".", strand, ".", "gene_name " + '"' + gene + '"; transcript_id ' + '"' + transcript + '";']) + "\n")