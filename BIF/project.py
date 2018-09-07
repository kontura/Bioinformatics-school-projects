import sys

def chromozomes_legth():
    chrom_lens = [ 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415 ]
    total = 0
    for chrom_len in chrom_lens:
        total += chrom_len
    return total

class GeneBiotypes:
    rotein_coding                      =  0
    G_C_gene                           =  1
    G_D_gene                           =  2
    G_J_gene                           =  3
    G_V_gene                           =  4
    R_C_gene                           =  5
    R_D_gene                           =  6
    R_J_gene                           =  7
    R_V_gene                           =  8
    nRNA                               =  9
    RNA                                =  10
    noRNA                              =  11
    iRNA                               =  12
    isc_RNA                            =  13
    incRNA                             =  14
    on_coding                          =  15
    rocessed_transcript                =  16
    ntisense                           =  17
    prime_overlapping_ncrna            =  18
    ense_intronic                      =  19
    ense_overlapping                   =  20
    nown_ncrna                         =  21
    nitary_pseudogene                  =  22
    G_C_pseudogene                     =  23
    ranslated_processed_pseudogene     =  24
    olymorphic_pseudogene              =  25
    R_J_pseudogene                     =  26
    G_J_pseudogene                     =  27
    R_V_pseudogene                     =  28
    G_V_pseudogene                     =  29
    seudogene                          =  30
    nprocessed_pseudogene              =  31
    ranscribed_unprocessed_pseudogene  =  32
    ranslated_unprocessed_pseudogene   =  33
    ranscribed_processed_pseudogene    =  34
    rocessed_pseudogene                =  35
    ranscribed_unitary_pseudogene      =  36

    coding_genes_names = ["protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
    small_non_coding_genes_names = ["snRNA", "rRNA", "snoRNA", "miRNA", "misc_RNA"]
    long_non_coding_genes_names = ["lincRNA", "non_coding", "processed_transcript", "antisense", "3prime_overlapping_ncrna", "sense_intronic", "sense_overlapping", "known_ncrna"]
    pseudo_genes_names = ["unitary_pseudogene", "IG_C_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "TR_J_pseudogene", "IG_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "transcribed_processed_pseudogene", "processed_pseudogene", "transcribed_unitary_pseudogene"]


class GeneFeatures:
    transcript = 0
    CDS = 1

class GRange:
    def __init__(self, chromozom, start, stop, feature):
        self.chrom = chromozom
        self.start = start
        self.stop = stop
        self.feature = feature

class GRangeList:
    def __init__(self):
        self.list = []

    def amount_of_bp(self):
        amount = 0
        for rang in self.list:
            amount = amount + (rang.stop - rang.start)
        return amount

    def amount_of_bp_of_feature(self, feature):
        amount = 0
        for rang in self.list:
            if rang.feature == feature:
                amount = amount + (rang.stop - rang.start)
        return amount

    def reduce(self):
        c_chromozom = 0
        c_start = 0
        c_stop = 0
        for rang in self.list:
            c_chromozom, c_start, c_stop = self.update_intervals(c_chromozom, c_start, c_stop, rang)

    def update_intervals(self, current_chromozom, current_start, current_stop, rang):
        if rang.chrom == current_chromozom:
            if current_stop < rang.start:
                current_start = rang.start
                return current_chromozom, current_start, rang.stop
            else:
                if current_stop < rang.stop:
                    rang.start = current_stop+1 #potential +1 should be here, cause now they overlap by one, the transition point
                    return current_chromozom, current_start, rang.stop
                else:
                    rang.start = 0
                    rang.stop = 0
                    return current_chromozom, current_start, current_stop
                    # or possible delete?
        else:
            return rang.chrom, 0, 0

    def reduce_feature(self, feature):
        c_chromozom = 0
        c_start = 0
        c_stop = 0
        for rang in self.list:
            if rang.feature == feature:
                c_chromozom, c_start, c_stop = self.update_intervals(c_chromozom, c_start, c_stop, rang)


    def add_range(self, gr):
        self.list.append(gr)

    def count(self):
        return len(self.list)

    def count_feature(self, feature):
        count = 0
        for rang in self.list:
            if rang.feature == feature:
                count += 1
        return count


def argument():
    if len(sys.argv) != 2:
        print("Usage: project.py <file.gtf>")
        quit()
    else:
        return sys.argv[1]

def find_gene_biotype(attributes):
    #possible heuristic, first try most often fields, if not the search all
    for attribute in attributes:
        if attribute.find('gene_biotype') != -1:
            text = attribute.split('"')
            #matching witout first letter, cause one starts with "3", and that can be no attr name
            return biotype_id(text[1])

def parse_line(line, lists):
   values = line.split('\t') 
   if (len(values) == 9):
       attributes = values[8].split(';')
       gene_biotype = find_gene_biotype(attributes)
       if gene_biotype != None:
           r = GRange(values[0], int(values[3]), int(values[4]), getattr(GeneFeatures, values[2], None))
           lists[gene_biotype].add_range(r)

def parse_by_line(file_name, dict_of_GRangeLists):
    with open(file_name) as file_gtf:
        for line in file_gtf:
            parsed_line = parse_line(line, dict_of_GRangeLists)

def initialize_lists():
    dict_list = {}
    for num in range(0, 37):
        dict_list[num] = GRangeList()
    return dict_list

def ssv_print(list):
    for value in list:
        sys.stdout.write(str(value))
        sys.stdout.write(";")
    sys.stdout.write("\n")
    sys.stdout.flush()

def biotype_id(name):
    return getattr(GeneBiotypes, name[1:], None)

def print_ssv_array(array, total_chrom_len):
    total_size_bp = 0
    total_count = 0
    total_coverage = 0
    for name in array:
        gene_biotype = biotype_id(name)
        dict_of_GRangeLists[gene_biotype].reduce()
        size_bp = dict_of_GRangeLists[gene_biotype].amount_of_bp()
        count = dict_of_GRangeLists[gene_biotype].count()
        coverage = (size_bp / total_chrom_len) * 100
        ssv_print([name, count, size_bp, coverage])
        total_size_bp += size_bp
        total_count += count
        total_coverage += coverage
    ssv_print(["<b>Summary</b>", total_count, total_size_bp, total_coverage])


def csv_out(list_dict):
    dict_of_GRangeLists[0].reduce_feature(getattr(GeneFeatures, "transcript"))
    transcript_size_bp = dict_of_GRangeLists[0].amount_of_bp_of_feature(getattr(GeneFeatures, "transcript"))

    dict_of_GRangeLists[0].reduce_feature(getattr(GeneFeatures, "CDS"))
    CDS_size_bp = dict_of_GRangeLists[0].amount_of_bp_of_feature(getattr(GeneFeatures, "CDS"))

    total_chrom_len = chromozomes_legth()

    for key in dict_of_GRangeLists:
        dict_of_GRangeLists[key].reduce()

    ssv_print(["Coding genes", "Count", "Size [bp]", "G.cov. [%]"])
    print_ssv_array(GeneBiotypes.coding_genes_names, total_chrom_len)
    ssv_print([])
    ssv_print(["Small non-coding genes", "Count", "Size [bp]", "G.cov. [%]"])
    print_ssv_array(GeneBiotypes.small_non_coding_genes_names, total_chrom_len)
    ssv_print([])
    ssv_print(["Long non-coding genes", "Count", "Size [bp]", "G.cov. [%]"])
    print_ssv_array(GeneBiotypes.long_non_coding_genes_names, total_chrom_len)
    ssv_print([])
    ssv_print(["Pseudo genes", "Count", "Size [bp]", "G.cov. [%]"])
    print_ssv_array(GeneBiotypes.pseudo_genes_names, total_chrom_len)
    ssv_print([])

    ssv_print(["Proteing coding genes", "Count", "Size [bp]", "G.cov. [%]"])
    coverage = (transcript_size_bp / total_chrom_len) * 100
    ssv_print(["Coding transcripts", dict_of_GRangeLists[0].count_feature(getattr(GeneFeatures, "transcript")), transcript_size_bp, coverage])
    coverage = (CDS_size_bp / total_chrom_len) * 100
    ssv_print(["CDS", dict_of_GRangeLists[0].count_feature(getattr(GeneFeatures, "CDS")), transcript_size_bp, coverage])

file_name = argument()
dict_of_GRangeLists = initialize_lists()
parse_by_line(file_name, dict_of_GRangeLists)
csv_out(dict_of_GRangeLists)
