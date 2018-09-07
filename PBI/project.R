library(AnnotationHub)
library(rlist)
library(ggplot2)
library(biomaRt)
library(Biostrings)
library(pqsfinder)

pqs_size_filer = function(pqs_data, size){
  scores = score(pqs_data)
  widths = width(pqs_data)
  c = scores[widths==12]
  result = sum(c)/length(c)
  # if Pqsfinder didnt find any matches set score to 0, another alternative would be to remove these sequences from result. Pqsfinder with default scoring system, will not find sequences, where there are multiple zero len inner loops, where as our regular expression will. We could also change scoring system of Pqsfinder.
  if (is.nan(result)){
    return(0)
  }
  else{
    return(result)
  }
}

sort_into_data_frame = function(ratios, seqs){
  ratios.withoutNULLs = list.clean(ratios, fun=is.null, recursive=FALSE)
  seqs.withoutNULLs = list.clean(seqs, fun=is.null, recursive=FALSE)

  df = data.frame(id=names(ratios.withoutNULLs), ratio=unlist(ratios.withoutNULLs))
  df$seq = seqs.withoutNULLs
  df = df[ order(-df$ratio), ]

  return(df)
}

printAllKLen = function(alphabet, prefix, n, k){
  if (k == 0){
    return(prefix)
  }
  l = c()
  for(i in 1:n){
    new_prefix = paste(prefix, alphabet[i], sep="")
    l = c(printAllKLen(alphabet, new_prefix, n, k-1), l)
  }
  return(l)
}

expected_number_of_occurances = function(substring_len, string_len, probabilities){
  #https://math.stackexchange.com/questions/220521/expected-number-of-times-random-substring-occurs-inside-of-larger-random-string
  sum = sum(`^`(probabilities, 2))
  sum = `^`(sum, substring_len)
  return(sum*(string_len - substring_len +1))
}

filter_matches_by_len = function(x, len){
  if (x[[1]][1] != -1){
    df = data.frame(start=unlist(x), length=unlist(attr(x[[1]], "match.length")))

    #longer = lapply(x, function(y) )
    subdf = subset(df, length == len)

    return(if(nrow(subdf)){subdf}else{NULL})
  }else{
    return(NULL)
  }
}

generate_uniq_AAs = function(){
  print("Generating amino acids complient to specified regex")
  nucleotides = printAllKLen(c('A', 'T', 'G', 'C'), "", 4, 8) #we only need 8, pairs of GG will be always pre/appended to start/end
  complient_ncds = lapply(nucleotides, function(x) gregexpr(".{0,4}GG.{0,4}GG.{0,4}", x))
  filtered_com_ncds = lapply(complient_ncds, function(x) filter_matches_by_len(x, 8))
  ncds_with_match = mapply(function(x, y) if(!is.null(y)){x}, nucleotides, filtered_com_ncds)
  ncds_with_match.clean = list.clean(ncds_with_match, fun=is.null, recursive=FALSE)
  ncds_with_match.clean.added_GGs = lapply(ncds_with_match.clean, function(x) paste("GG","GG",sep=x))
  translated = lapply(ncds_with_match.clean.added_GGs, function(x) translate(DNAString(x)))
  c_m_str = lapply(translated, function(x) toString(x))
  u_q = unique(unlist(c_m_str, use.names=FALSE))
  cat(u_q2, sep="|")
  return(u_q)
}

generate_uniq_AAs_old_alternative = function(human_cdna_strings_with_match, human_cdna_matches_of_12){
  # OLD NOT USED, but gave the same (at least same amount) of uniq 4AA seqs for regex
  matched_aminos = mapply(function(x,y) if(!is.null(x)){translate(substr(x,head(y$start,1), head(y$start,1)+11))}, human_cdna_strings_with_match, human_cdna_matches_of_12)
  c_m = list.clean(matched_aminos, fun=is.null, recursive=FALSE)
  c_m_str = lapply(c_m, function(x) toString(x))
  u_q = unique(unlist(c_m_str, use.names=FALSE))
}

##############################################
                 #  MAIN  #
##############################################

ah = AnnotationHub()
print("Quering AnnotationHub")
ah2 = query(ah, c("fasta", "homo sapiens", "Ensembl"))
ah3 = query(ah, c("Saccharomyces cerevisiae", "fasta", "Ensembl"))
ah4 = query(ah, c("Arabidopsis thaliana"))
print("Fetching human fasta file")
human_cdna_fasta = ah2[["AH18522"]]
human_pep_fasta = ah2[["AH21570"]]

print("Fetching yeast fasta file")
yeast_cdna_fasta = ah3[["AH21745"]]
yeast_pep_fasta = ah3[["AH21750"]]

print("Fetching arabidopsis fasta file")
if (!file.exists("arabidopsis_cdna.fa")){
  download.file("ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/TAIR10_cdna_20101214_updated", "arabidopsis_cdna.fa", "wget")
}
if (!file.exists("arabidopsis_pep.fa")){
  download.file("ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/TAIR10_pep_20101214_updated", "arabidopsis_pep.fa", "wget")
}

print("Loading sequences from fasta")
arabidopsis_cdna = readDNAStringSet("arabidopsis_cdna.fa")
arabidopsis_pep = readAAStringSet("arabidopsis_pep.fa")
human_cdna = readDNAStringSet(path(human_cdna_fasta))
human_pep = readAAStringSet(path(human_pep_fasta))
yeast_cdna = readDNAStringSet(path(yeast_cdna_fasta))
yeast_pep = readAAStringSet(path(yeast_pep_fasta))

print("Matching G-Quadruplexes in cdna")
regex_cdna = "(GG..GG.GG.GG)|(GG.GG..GG.GG)|(GG.GG.GG..GG)|(GGGGGG....GG)|(GGGG....GGGG)|(GG....GGGGGG)|(GG..GG..GGGG)|(GGGG..GG..GG)|(GG..GGGG..GG)|(GG...GG.GGGG)|(GG...GGGG.GG)|(GG.GGGG...GG)|(GG.GG...GGGG)|(GGGG...GG.GG)|(GGGG.GG...GG)"
print("for human..")
human_cdna_matches = lapply(human_cdna, function(x) gregexpr(regex_cdna, x))
print("for yeast..")
yeast_cdna_matches = lapply(yeast_cdna, function(x) gregexpr(regex_cdna, x))
print("for arabidopsis..")
arabidopsis_cdna_matches = lapply(arabidopsis_cdna, function(x) gregexpr(regex_cdna, x))
print("Matching G-Quadruplexes in peptitedes")
#AA sequences obtained by generate_uniq_AAs function, AA sequences with * (stop codon) are removed
regex = "G(PGG|PAG|PGR|PGW|PVG|PEG|RGG|RPG|RRR|RRG|RRW|RLG|RQG|RAG|RGR|RGW|RVG|REG|RSG|RWR|RWG|RWW|RTG|RMG|RKG|LGG|LAG|LGR|LGW|LVG|LEG|HGG|QAG|QGR|QGG|QGW|QVG|QEG|AGG|AAG|AGR|AGW|AVG|AEG|GPG|GRR|GRG|GRW|GLG|GQG|GAG|GGR|GGG|GGW|GVG|GEG|GSG|GWR|GWG|GWW|GTG|GMG|GKG|GAR|GAW|GVR|GVW|GDR|GDG|GDW|GER|GEW|VGG|VAG|VGR|VGW|VVG|VEG|DGG|EAG|EGR|EGG|EGW|EVG|EEG|SGG|SAG|SGR|SGW|SVG|SEG|CGG|WPG|WRR|WRG|WRW|WLG|WQG|WAG|WGR|WGG|WGW|WVG|WEG|WSG|WWR|WWG|WWW|WTG|WMG|WKG|FGG|YGG|TGG|TAG|TGR|TGW|TVG|TEG|IGG|MAG|MGR|MGG|MGW|MVG|MEG|NGG|KAG|KGR|KGG|KGW|KVG|KEG|APG|ARR|ARG|ARW|ALG|AQG|ASG|AWR|AWG|AWW|ATG|AMG|AKG|AAR|AAW|AVR|AVW|ADR|ADG|ADW|AER|AEW|GPR|GPW|GLR|GLW|GHR|GHG|GHW|GQR|GQW|GSR|GSW|GCR|GCG|GCW|GFR|GFG|GFW|GYR|GYG|GYW|GTR|GTW|GIR|GIG|GIW|GMR|GMW|GNR|GNG|GNW|GKR|GKW|VPG|VRR|VRG|VRW|VLG|VQG|VSG|VWR|VWG|VWW|VTG|VMG|VKG|VAR|VAW|VVR|VVW|VDR|VDG|VDW|VER|VEW|DPG|DRR|DRG|DRW|DLG|DQG|DAG|DGR|DGW|DVG|DEG|DSG|DWR|DWG|DWW|DTG|DMG|DKG|EPG|ERR|ERG|ERW|ELG|EQG|EAR|EAW|EVR|EVW|EDR|EDG|EDW|EER|EEW|ESG|EWR|EWG|EWW|ETG|EMG|EKG)"
print("for human..")
human_pep_matches = lapply(human_pep, function(x) gregexpr(regex, x))
print("for yeas..")
yeast_pep_matches = lapply(yeast_pep, function(x) gregexpr(regex, x))
print("for arabidopsis..")
arabidopsis_pep_matches = lapply(arabidopsis_pep, function(x) gregexpr(regex, x))

print("Filtering pep seqs with G-Quadruplexes")
#coverts to data.frame, because of legacy reasons check length
human_pep_matches_of_4 = lapply(human_pep_matches, function(x) filter_matches_by_len(x, 4))
yeast_pep_matches_of_4 = lapply(yeast_pep_matches, function(x) filter_matches_by_len(x, 4))
arabidopsis_pep_matches_of_4 = lapply(arabidopsis_pep_matches, function(x) filter_matches_by_len(x, 4))

human_pep_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, human_pep, human_pep_matches_of_4)
yeast_pep_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, yeast_pep, yeast_pep_matches_of_4)
arabidopsis_pep_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, arabidopsis_pep, arabidopsis_pep_matches_of_4)

print("Filtering seqs with G-Quadruplexes to length 12")
#coverts to data.frame, because of legacy reasons check length
human_cdna_matches_of_12 = lapply(human_cdna_matches, function(x) filter_matches_by_len(x, 12))
yeast_cdna_matches_of_12 = lapply(yeast_cdna_matches, function(x) filter_matches_by_len(x, 12))
arabidopsis_cdna_matches_of_12 = lapply(arabidopsis_cdna_matches, function(x) filter_matches_by_len(x, 12))

human_cdna_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, human_cdna, human_cdna_matches_of_12)
yeast_cdna_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, yeast_cdna, yeast_cdna_matches_of_12)
arabidopsis_cdna_strings_with_match = mapply(function(x, y) if(!is.null(y)){x}, arabidopsis_cdna, arabidopsis_cdna_matches_of_12)

print("Counting expected occurance rates")
human_expected_occurs_cdna = lapply(human_cdna_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(12, length(x), alphabetFrequency(x, as.prob=TRUE, baseOnly=TRUE))}else{NULL})
yeast_expected_occurs_cdna = lapply(yeast_cdna_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(12, length(x), alphabetFrequency(x, as.prob=TRUE, baseOnly=TRUE))}else{NULL})
arabidopsis_expected_occurs_cdna = lapply(arabidopsis_cdna_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(12, length(x), alphabetFrequency(x, as.prob=TRUE, baseOnly=TRUE))}else{NULL})

human_expected_occurs_pep = lapply(human_pep_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(4, length(x), alphabetFrequency(x, as.prob=TRUE))}else{NULL})
yeast_expected_occurs_pep = lapply(yeast_pep_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(4, length(x), alphabetFrequency(x, as.prob=TRUE))}else{NULL})
arabidopsis_expected_occurs_pep = lapply(arabidopsis_pep_strings_with_match, function(x) if(!is.null(x)){expected_number_of_occurances(4, length(x), alphabetFrequency(x, as.prob=TRUE))}else{NULL})

print("Counting occurances")
human_cdna_occurs = lapply(human_cdna_matches_of_12, function(x) if(!is.null(x)){nrow(x)})
yeast_cdna_occurs = lapply(yeast_cdna_matches_of_12, function(x) if(!is.null(x)){nrow(x)})
arabidopsis_cdna_occurs = lapply(arabidopsis_cdna_matches_of_12, function(x) if(!is.null(x)){nrow(x)})

human_pep_occurs = lapply(human_pep_matches_of_4, function(x) if(!is.null(x)){nrow(x)})
yeast_pep_occurs = lapply(yeast_pep_matches_of_4, function(x) if(!is.null(x)){nrow(x)})
arabidopsis_pep_occurs = lapply(arabidopsis_pep_matches_of_4, function(x) if(!is.null(x)){nrow(x)})

print("Computing rations")
human_cdna_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, human_cdna_occurs, human_expected_occurs_cdna)
yeast_cdna_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, yeast_cdna_occurs, yeast_expected_occurs_cdna)
arabidopsis_cdna_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, arabidopsis_cdna_occurs, arabidopsis_expected_occurs_cdna)

human_pep_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, human_pep_occurs, human_expected_occurs_pep)
yeast_pep_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, yeast_pep_occurs, yeast_expected_occurs_pep)
arabidopsis_pep_ratio = mapply(function(x,y) if(!is.null(x)){x/y}, arabidopsis_pep_occurs, arabidopsis_expected_occurs_pep)

print("Sorting by rations")
human_cdna_df = sort_into_data_frame(human_cdna_ratio, human_cdna_strings_with_match)
yeast_cdna_df = sort_into_data_frame(yeast_cdna_ratio, yeast_cdna_strings_with_match)
arabidopsis_cdna_df = sort_into_data_frame(arabidopsis_cdna_ratio, arabidopsis_cdna_strings_with_match)

human_pep_df = sort_into_data_frame(human_pep_ratio, human_pep_strings_with_match)
yeast_pep_df = sort_into_data_frame(yeast_pep_ratio, yeast_pep_strings_with_match)
arabidopsis_pep_df = sort_into_data_frame(arabidopsis_pep_ratio, arabidopsis_pep_strings_with_match)

print("Printing top 50")
mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][4], ratio)), head(human_cdna_df$ratio,50), head(human_cdna_df$id,50))
mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][4], ratio)), head(yeast_cdna_df$ratio,50), head(yeast_cdna_df$id,50))
mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][1], ratio)), head(arabidopsis_cdna_df$ratio,50), head(arabidopsis_cdna_df$id,50))

mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][4], ratio)), head(human_pep_df$ratio,50), head(human_pep_df$id,50))
mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][4], ratio)), head(yeast_pep_df$ratio,50), head(yeast_pep_df$id,50))
mapply(function(ratio, name) print(sprintf("%s, ratio: %f", strsplit(toString(name), ' ')[[1]][1], ratio)), head(arabidopsis_pep_df$ratio,50), head(arabidopsis_pep_df$id,50))

print("Computing scores with pqsfinder")
print("for human..")
human_pqss = lapply(human_cdna_strings_with_match, function(x) if(!is.null(x)){invisible(pqsfinder(x, max_len=12, run_max_len=2, run_min_len=2, loop_min_len=0,loop_max_len=4, strand="+", overlapping=TRUE, run_re="GG", min_score=1))}else{NULL})
print("for yeas..")
yeast_pqss = lapply(yeast_cdna_strings_with_match, function(x) if(!is.null(x)){invisible(pqsfinder(x, max_len=12, run_max_len=2, run_min_len=2, loop_min_len=0,loop_max_len=4, strand="+", overlapping=TRUE, run_re="GG", min_score=1))}else{NULL})
print("for arabidopsis..")
arabidopsis_pqss = lapply(arabidopsis_cdna_strings_with_match, function(x) if(!is.null(x)){invisible(pqsfinder(x, max_len=12, run_max_len=2, run_min_len=2, loop_min_len=0,loop_max_len=4, strand="+", overlapping=TRUE, run_re="GG", min_score=1))}else{NULL})

human_pqss_cleaned = list.clean(human_pqss, fun=is.null, recursive=FALSE)
yeast_pqss_cleaned = list.clean(yeast_pqss, fun=is.null, recursive=FALSE)
arabidopsis_pqss_cleaned = list.clean(arabidopsis_pqss, fun=is.null, recursive=FALSE)

human_cdna_df$pqss = unlist(lapply(human_pqss_cleaned, function(x) pqs_size_filer(x, 12)), use.names=F)
yeast_cdna_df$pqss = unlist(lapply(yeast_pqss_cleaned, function(x) pqs_size_filer(x, 12)), use.names=F)
arabidopsis_cdna_df$pqss = unlist(lapply(arabidopsis_pqss_cleaned, function(x) pqs_size_filer(x, 12)), use.names=F)

print("All pqsFinder scores are 32. There is no point in showing graphs or computing Pearson, however the code is present, just commented out.")

#print("Displaying relationship between ratio and pqsFinder score")
#human_cdna_df.filteredZeros = human_cdna_df[human_cdna_df$pqss != 0, ]
#yeast_cdna_df.filteredZeros = yeast_cdna_df[yeast_cdna_df$pqss != 0, ]
#arabidopsis_cdna_df.filteredZeros = arabidopsis_cdna_df[arabidopsis_cdna_df$pqss != 0, ]
#
#human_graph = ggplot(data=human_cdna_df.filteredZeros, aes(x=ratio, y=pqss)) + geom_point()
#yeast_graph = ggplot(data=yeast_cdna_df.filteredZeros, aes(x=ratio, y=pqss)) + geom_point()
#arabidopsis_graph = ggplot(data=arabidopsis_cdna_df.filteredZeros, aes(x=ratio, y=pqss)) + geom_point()
#
#print("Computing Pearson")
#human_cor_sig_test = cor.test(human_cdna_df.filteredZeros$pqss, human_cdna_df.filteredZeros$ratio, method = "pearson")
#yeast_cor_sig_test = cor.test(yeast_cdna_df.filteredZeros$pqss, yeast_cdna_df.filteredZeros$ratio, method = "pearson")
#arabidopsis_cor_sig_test = cor.test(arabidopsis_cdna_df.filteredZeros$pqss, arabidopsis_cdna_df.filteredZeros$ratio, method = "pearson")
#
#print("For Human: Pearson correlation coefficient %f" human_cor_sig_test$estimate)
#if (human_cor_sig_test$p.value > 0.05){
#  print("Not statistically significant")
#}
#
#print("For Yeast: Pearson correlation coefficient %f" yeast_cor_sig_test$estimate)
#if (yeast_cor_sig_test$p.value > 0.05){
#  print("Not statistically significant")
#}
#
#print("For Arabidopsis: Pearson correlation coefficient %f" arabidopsis_cor_sig_test$estimate)
#if (arabidopsis_cor_sig_test$p.value > 0.05){
#  print("Not statistically significant")
#}
