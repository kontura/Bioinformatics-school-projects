############

#Incomplete version
#Rozpracovana verze

###########


library(GenomicRanges)
library(ggplot2)
library(biomaRt)
library(Gviz)
library(Biostrings)
library(gridGraphics)
library(gridExtra)

plot_track = function(x){
  plotTracks(x)
  gt = grid.grab()
  return(gt)
}

parse_exons = function(grange){
  exons = lapply(grange$unformatted_exons[[1]]$gene_exon, function(x) matchPattern(DNAString(x), DNAString(grange$sequence[[1]])))
  exons_vec = lapply(exons, function(x) start(x) )
  exons.df = as.data.frame(matrix(exons_vec))
  exons.df$end = lapply(exons, function(x) end(x) )
  i = IRanges(start=unlist(exons.df$V1), end=unlist(exons.df$end))
  gRa = GRanges(seqnames=rep(paste("chr",seqnames(grange), sep=""), length(i)), ranges=i)
  return(gRa)
}

cg_seq_graph = function(sequence, name){
  positions = gregexpr('CG', sequence)

  values = c(1:nchar(sequence))
  distance_to_closest_CG = lapply(values, function(x) min(abs(positions[[1]] - x)))
  distance = unlist(distance_to_closest_CG)
  myDF = data.frame(distance)
  myDF = cbind(position = rownames(myDF), myDF)
  graph = ggplot(data=myDF, aes(x=position, y=distance)) + geom_point(aes(stroke=0.01)) + scale_x_discrete(breaks=seq(from=0, to=length(myDF$position), by=round(length(myDF$position)/15)), limits=seq(0,length(myDF$position))) + ggtitle(name)
  return(ggplotGrob(graph))
#  ggsave("plot.png", plot=graph)
}

graph_to_grob = function(graph){
  g1= do.call(rbind, c(c(graph), size="first"))
  g1$widths = do.call(unit.pmax, lapply(c(graph), "[[", "widths"))

  return(g1)
}


ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
print("Fetching all genes")
data = getBM(attributes=c("chromosome_name", "transcript_start", "transcript_end", "hgnc_symbol", "strand", "ensembl_gene_id"), mart=ensembl)

print("Selecting only HOX genes")
data = subset(data, grepl("HOX", data$hgnc_symbol))

ir = IRanges(start=data$transcript_start, end=data$transcript_end)
gr = GRanges(seqnames=data$chromosome_name, ranges = (ir), name=data$hgnc_symbol, strand=data$strand, ensembl_gene_id=data$ensembl_gene_id)

gr$len = width(gr)

gr.ordered = gr[order(-gr$len), ]
gr.uniq = gr.ordered[!duplicated(gr.ordered$name), ]
gr = gr.uniq

#loaded = read.table("./hox_tableBrowser")
#ir.loaded = IRanges(start=c(loaded[3])$V3, end=c(loaded[4])$V4)
#gr.loaded = GRanges(seqnames=substring(loaded$V1, 4), ranges = ir.loaded, name=c(loaded[5])$V5, strand=loaded$V2)
#
#gr.loaded$len = width(gr.loaded)
#
#
#gr.loaded.ordered = gr.loaded[order(-gr.loaded$len), ]
#
#gr.loaded.uniq = gr.loaded.ordered[!duplicated(gr.loaded.ordered$name), ]
#
#histogram generation
png("histogram.png")
hist(gr$len, breaks=10, xlab="Length", main="DNA Lengths")
dev.off()

print("Fetching whole sequences for selected HOX genes")
seqs = lapply(gr, function(x) getSequence(id=x$ensembl_gene_id, mart=ensembl, seqType="gene_exon_intron", type="ensembl_gene_id"))

print("Creating CG distance graphs")
gr$sequence = lapply(seqs, function(x) x$gene_exon_intron[1])
gr$graph = lapply(gr, function(x) if (runValue(strand(x)) == "+") cg_seq_graph(x$sequence, x$name) else cg_seq_graph(chartr("ATCG", "TAGC", x$sequence), x$name))

print("Fetching exons only")
gr$unformatted_exons = lapply(gr, function(x) getSequence(id=x$ensembl_gene_id, mart=ensembl, seqType="gene_exon", type="ensembl_gene_id"))
print("Parsing exons")
gr$GRange_exons = lapply(gr, function(x) parse_exons(x))
gr$exon_count = unlist(lapply(gr, function(x) unlist(length(x$GRange_exons[[1]]))))

print("Creating tracks")
tracks = lapply(gr, function(x) AnnotationTrack(reduce(disjoin(x$GRange_exons[[1]], ignore.strand=TRUE, with.revmap=TRUE)), name=x$name))

print("Plotting tracks")
plotted_tracks = lapply(tracks, function(x) plot_track(x))

print("Printing to PDF")
graphs = lapply(gr, function(x) graph_to_grob(x$graph))

pdf()
for (i in 1:length(graphs)){
  grid.arrange(graphs[[i]], plotted_tracks[[i]], ncol=1)
}
dev.off()

print("All done!")

b = lapply(gr$unformatted_exons, function(x) x[[1]])
c = paste(b)
d = gsub("[\n]", "", c)
df = data.frame(name=gr$name, chromosome=seqnames(gr), start=start(gr)-1, end=end(gr), exon_count=gr$exon_count, exon_sequences=d)
write.table(df, file="hox_table", quote=F, sep="\t", row.names=F, col.names=F)
