rm(list=ls())
.libPaths(c("/oak/stanford/groups/ccurtis2/users/ydzhao/R_Library_4.2/R_Library_4.2",.libPaths()[2]))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Construct the function Granges
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

library(gUtils)
library(gTrack)
library(strawr)
library(igraph)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(ggrepel)
library(DCG)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene

tmpinf1 = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/ecDNA/output/SPCG-OS152_11M_classification_bed_files/SPCG-OS152_11M_amplicon3_ecDNA_1_intervals.bed"
tmpinf1 = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/ecDNA_refine/Data/Info/H000530_amp2_cyc3_hg38.bed"
res1 = read.table(tmpinf1,sep="\t",quote=NULL)
res = rbind(res1)
colnames(res) = c("chr1","pos1","pos2")
res[,"dir1"] = rep("-",nrow(res))
res[,"dir2"] = rep("+",nrow(res))
res[,"classification"] = rep("DUP",nrow(res))
res[,"chr2"] = res[,"chr1"]
res = res[,c("chr1","pos1","dir1","chr2","pos2","dir2","classification")]

info = res
row.names(info) = seq(1,nrow(info),1)

target = c("chr17")
tag = which(info$chr1 %in% target)
info = info[tag,]
row.names(info) = seq(1,nrow(info),1)
info = info[order(info$pos1),]
#info = info[c(1,2),]
#info = info[c(5,6,7),]

tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Genome/Hg38/Hg38_gene_position_regular_chr.bed"
genes = read.table(tmpinf,sep="\t",quote=NULL)
genes[,"V5"] = rep(20,nrow(genes))
tag = genes$V4 == "ERBB2"
genes = genes[tag,]
genes[,"Group"] = "Gene"
info_sup = data.frame(chr1="chr17",pos1=genes$V2,dir1="+",chr2="chr17",pos2=genes$V3,dir2="+",classification="DUP")
info = rbind(info,info_sup)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Reconstruct the SV 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SV_summary = info[,c("chr1","chr2")]
SV_summary = unique(SV_summary)
SV_summary$Label = gsub("chr","",SV_summary$chr1)
SV_summary$Label = as.numeric(SV_summary$Label)
SV_summary = SV_summary[order(SV_summary$Label),]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#chr17 to chr17 translocation
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
tag1 = info$chr1 == "chr17" & info$chr2 == "chr17"
tag2 = info$chr1 == "chr17" & info$chr2 == "chr17"

info_section_1 = info[tag1,]
info_section_2 = info[tag2,]
tmp = rbind(info_section_1,info_section_2)

tmp = unique(tmp)

distance_sv = 500000

for(i in 1 : nrow(tmp))
{
	cat("\r",i)
	
	chr1 = tmp[i,"chr1"]
	sta1 = tmp[i,"pos1"]
	
	chr2 = tmp[i,"chr2"]
	sta2 = tmp[i,"pos2"]
	
	sv_interval_1 = data.frame(chr=chr1, sta=sta1 - distance_sv, end= sta1 + distance_sv,classification= tmp$classification[i])
	sv_interval_2 = data.frame(chr=chr2, sta=sta2 - distance_sv, end=sta2 + distance_sv,classification= tmp$classification[i])
	sv_interval = rbind(sv_interval_1,sv_interval_2)
	
	if(i == 1)
	{
		chr10_chr10_trans_sv_interval = sv_interval
	
	}else{
	
		chr10_chr10_trans_sv_interval = rbind(chr10_chr10_trans_sv_interval ,sv_interval)
	
	}
}

sv_interval_final = rbind(chr10_chr10_trans_sv_interval)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#sv edge
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
target = seq(1,nrow(sv_interval_final),2)
sv_edge = matrix(0,length(target),4)
colnames(sv_edge) = c("from","to","col","type")
sv_edge = as.data.frame(sv_edge)

for(i in 1 : length(target))
{
	cat("\r",i)
	
	sv_edge[i,"from"] = target[i]
	sv_edge[i,"to"] = target[i]+1
	sv_edge[i,"type"] = sv_interval_final$classification[target[i]]
}
sv_edge$col = rep("red",nrow(sv_edge))
sv_edge$h = 10

target = seq(1,nrow(sv_interval_final),1)
graph = matrix(0, length(target),length(target))

for(i in 1 : nrow(sv_edge))
{
	cat("\r",i)
	row_pos = sv_edge$from[i]
	col_pos = sv_edge$to[i]
	
	graph[row_pos,col_pos] =1

}

sv_interval_final = peakDF2GRanges(sv_interval_final)

Temp_sv = as.data.frame(sv_interval_final)
target = c("chr17")

Region = rep(0,length(target))

for(i in 1 : length(target))
{
	tag = Temp_sv$seqnames == target[i]
	tmp = Temp_sv[tag,]
	min_start = max((min(tmp$start) - 100000),1)
	max_end = max(tmp$end) + 100000
	
	Region[i] = paste0(target[i],":",min_start,":",max_end)

}
Region = gsub("chr","",Region)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Construct the matrix One time
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
options(scipen=999)
tmpinf = "/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/isabl_configs/2022-10-24_mapping_donor_ids.txt"
map = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = grep("H000530",map$Individual.System.ID)
map = map[tag,]
#TCGA-A2-A0EY

tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/hichip_libraries_metadata.txt"
annotation = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = annotation$submitter_id == "TCGA-A2-A0EY"
annotation[tag,]
#BRCA-8D1E6006

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#chr10 : chr10 interaction
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/hic/TCGA_HiChIP_hic/BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac.allValidPairs.hic"
#data = straw("NONE", tmpinf, Region[1], Region[1], "BP", 10000, matrix = "observed")
data = straw("KR", tmpinf, Region[1], Region[1], "BP", 10000, matrix = "OE")

data$From = paste0("chr17","_",data$x,"_",data$x+10000)
data$To = paste0("chr17","_",data$y,"_",data$y+10000)
data = data[,c("From","To","counts")]
tag = is.na(data$counts)
data = data[tag==0,]
chr10_chr10_interaction = data

data = rbind(chr10_chr10_interaction)

#############################
#Post process matrix
#############################
# Step 1: Create canonical (sorted) pair ID for each row
data$pair_id <- apply(data[, c("From", "To")], 1, function(x) paste(sort(x), collapse = "_"))

# Step 2: Average counts across symmetric pairs
data_avg <- aggregate(counts ~ pair_id, data = data, FUN = mean)

# Step 3: Recover From and To from the pair_id
split_cols <- do.call(rbind, strsplit(data_avg$pair_id, "_"))
data_avg$From <- paste0(split_cols[, 1],"_",split_cols[, 2],"_",split_cols[, 3])
data_avg$To   <- paste0(split_cols[, 4],"_",split_cols[, 5],"_",split_cols[, 6])

# Step 4: Reorder columns for clarity
data_avg <- data_avg[, c("From", "To", "counts")]

edgelist = data_avg
mygraph <- graph.data.frame(edgelist,directed=F)
heatMap = as_adjacency_matrix(mygraph,sparse=T,attr='counts',type='both' )
#isSymmetric(as.matrix(heatMap))
data = log2(heatMap+1)

#myoutf = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/HiC/3D_contact/3D_contact_chr6_chr8_One_Time_Simple.Rda"
#save(data,file=myoutf)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Build the node annotation
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#load("/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/HiC/3D_contact/Node_collection_chr10_chr11_chr19.Rda")

Node_collection = colnames(data)
chr = sapply(Node_collection,function(x) strsplit(x,"_")[[1]][1])
sta = sapply(Node_collection,function(x) strsplit(x,"_")[[1]][2])
end = sapply(Node_collection,function(x) strsplit(x,"_")[[1]][3])

chr = as.vector(chr)
sta = as.numeric(sta)
end = as.numeric(end)

node = data.frame(chr=chr,sta=sta,end=end)
window = node
node = peakDF2GRanges(node)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Build total window
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
window = as.data.frame(sv_interval_final)
target_chr = unique(window$seqnames)

for(i in 1 : length(target_chr))
{
	tag = window$seqnames == target_chr[i]
	tmp = window[tag,]
	
	target_sta = min(tmp$start, tmp$end)
	target_end = max(tmp$start, tmp$end)
	
	target_sta = target_sta - 10000
	target_end = target_end + 10000
	
	if(i==1)
	{
		tmp_window = data.frame(chr=target_chr[i],sta=target_sta,end=target_end)
		target_window = tmp_window
	
	}else{
		tmp_window = data.frame(chr=target_chr[i],sta=target_sta,end=target_end)
		target_window = rbind(target_window,tmp_window)
	
	}

}

target_window = peakDF2GRanges(target_window)
target_window = sv_interval_final
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Create the RUNX2 track
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#exons <- exons(Homo.sapiens, columns = c("SYMBOL"))

load("/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/Data/Whole_exon_data.Rda")

genes <- exons
genes <- as.data.frame(genes)
tag = grep("ERBB2",genes$SYMBOL)
genes = genes[tag,]
tag = genes$SYMBOL == "ERBB2"
genes = genes[tag,]
genes$exon = seq(1,nrow(genes),1)
genes = genes[,c("seqnames", "start","end","exon")]
genes = peakDF2GRanges(genes)
RSPO2_genes = genes

grl = GRangesList(RSPO2_genes=RSPO2_genes)
gt.genes = gTrack(grl)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Create the ATAC track
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(rtracklayer)
#tmpinf = "/oak/stanford/groups/howchang/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/ATAC_celline/bigwig/085_merge.bigwig"
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/ATAC/Bigwig/BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X012_S03_L029_B1_T1_P024.insertions.bw"
ATAC = import.bw(tmpinf)
#ATAC_node = subsetByOverlaps(ATAC, target_window, ignore.strand = TRUE) 
ATAC_node = ATAC
ATAC_node = as.data.frame(ATAC_node)

colnames(ATAC_node) = c("V1","V2","V3","V4","V5","score")
ATAC_node = peakDF2GRanges(ATAC_node)


##############################
#create the H3K27Ac track
##############################
#tmpinf = "/oak/stanford/groups/howchang/users/ydzhao/Resources/Pub_Dat/ChIPAltas/OS/H3K27Ac/SRX6480643.bw"
#H3K27Ac = import.bw(tmpinf)
#ATAC_node = subsetByOverlaps(ATAC, target_window, ignore.strand = TRUE) 

tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/H3K27Ac_1D/hg38/BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac_1D-signal.norm.bw"
H3K27Ac = import.bw(tmpinf)

H3K27Ac_node = H3K27Ac
H3K27Ac_node = as.data.frame(H3K27Ac_node)

colnames(H3K27Ac_node) = c("V1","V2","V3","V4","V5","score")
H3K27Ac_node$score = as.numeric(H3K27Ac_node$score)
H3K27Ac_node = peakDF2GRanges(H3K27Ac_node)

#################################
#Create the H3K27Ac track
##################################
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/H3K27Ac_peaks/BRCA-8D1E6006-85CB-484A-8B5C-30766D90137B-X005-S02-B1-T1_H3K27ac_peaks.narrowPeak"
H3K27Ac_peaks = read.table(tmpinf,sep="\t",quote=NULL)
H3K27Ac_peaks = H3K27Ac_peaks[,c("V1","V2","V3","V5")]
H3K27Ac_peaks =  peakDF2GRanges(H3K27Ac_peaks)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
op = findOverlaps(H3K27Ac_peaks,promoter)
op = as.data.frame(op)
H3K27Ac_peaks$Group = rep("enhancer",length(H3K27Ac_peaks))
label = unique(op$queryHits)
H3K27Ac_peaks$Group[label] = "promoter"

#################################
#Create the ATAC track
##################################
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/ATAC/Peaks/Cancer_specific_peaks/BRCA_peakCalls.txt"
ATAC_peaks = read.table(tmpinf,sep="\t",quote=NULL,header=T)
ATAC_peaks = ATAC_peaks[,c("seqnames","start","end","score")]
colnames(ATAC_peaks)[4] = "V5"
ATAC_peaks =  peakDF2GRanges(ATAC_peaks)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
op = findOverlaps(ATAC_peaks,promoter)
op = as.data.frame(op)
ATAC_peaks$Group = rep("enhancer",length(ATAC_peaks))
label = unique(op$queryHits)
ATAC_peaks$Group[label] = "promoter"

###################################
#Create the Gene Position
###################################
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Genome/Hg38/Hg38_gene_position_regular_chr.bed"
genes = read.table(tmpinf,sep="\t",quote=NULL)
genes[,"V5"] = rep(20,nrow(genes))
tag = genes$V4 == "ERBB2"
genes = genes[tag,]
genes[,"Group"] = "Gene"

tmpinf1 = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/ecDNA_refine/Data/Info/H000530_amp2_cyc3_hg38.bed"
res1 = read.table(tmpinf1,sep="\t",quote=NULL)
res = res1
colnames(res) = c("V1","V2","V3","V4","V5")
res[,"V4"] = rep("ecDNA",nrow(res))
res[,"V5"] = rep(20,nrow(res))
res[,"Group"] = rep("ecDNA",nrow(res))

genes = rbind(genes,res)
genes = peakDF2GRanges(genes)

#################################
#Create the CNV track
##################################
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/Clin/somatic.cnv.all.tsv"
tmpinf = "/oak/stanford/groups/ccurtis2/users/ydzhao/Resources/Pub_Dat/TCGA/ASCAT_ploidy_normalized_HiChIP_Samples.txt"
info = read.table(tmpinf,sep="\t",quote=NULL,header=T)
tag = info$case_submitter_id == "TCGA-AO-A03L"
info = info[tag,]

info = info[,c("Chromosome","Start","End","Relative_Copy_Number")]
colnames(info)[4] = "copyNumber"
info$Label = rep("deletion",nrow(info))
info$col = rep("#377EB8",nrow(info))
tag = info$copyNumber > 1.5
info$Label[tag] = "amplification"
info$col[tag] = "#E41A1C"
tag = info$copyNumber < 1.5 & info$copyNumber > 0.9
info$Label[tag] = "Normal"
info$col[tag] = "grey50"
info = peakDF2GRanges(info)

################################
#ecDNA track 
################################
tmpinf1 = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/ecDNA_refine/Data/Info/H000530_amp2_cyc3_hg38.bed"
#tmpinf2 = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/Data/ecDNA/output/SPCG-OS186_89CL_classification_bed_files/SPCG-OS186_89CL_amplicon2_ecDNA_2_intervals.bed"
#tmpinf3 = "/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/Sarcoma_ERV/Data/Data/ecDNA/output/SPCG-OS186_89CL_classification_bed_files/SPCG-OS186_89CL_amplicon2_ecDNA_3_intervals.bed"
ecDNA1 = read.table(tmpinf1,sep="\t",quote=NULL)
#ecDNA2 = read.table(tmpinf2,sep="\t",quote=NULL)
#ecDNA3 = read.table(tmpinf3,sep="\t",quote=NULL)
#ecDNA = rbind(ecDNA1,ecDNA2,ecDNA3)
ecDNA = ecDNA1
ecDNA = ecDNA[,c("seqnames","start","end")]
colnames(ecDNA) = c("V1","V2","V3")
ecDNA$V5 = ecDNA$V3 - ecDNA$V2
ecDNA =  peakDF2GRanges(ecDNA)
ecDNA$Group = "ecDNA"

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Change color
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(RColorBrewer)
Color_panel = brewer.pal(9,"Set1")[1:4]
tag1 = sv_edge$type == "DUP"
tag2 = sv_edge$type == "INV"
tag3 = sv_edge$type == "INTRX"
tag4 = sv_edge$type == "DEL"

sv_edge[tag1,"col"] = Color_panel[1]
sv_edge[tag2,"col"] = Color_panel[2]
sv_edge[tag3,"col"] = Color_panel[3]
sv_edge[tag4,"col"] = Color_panel[4]


sv_interval = sv_interval_final     #chr8_chr10_trans_sv_interval
sv_edge = sv_edge
node = node

tag = is.na(data)
data[tag] = 0
data = as.matrix(data)

bottom_value = quantile(data,0.05)
up_value = quantile(data,0.95)
tag1 = data <= bottom_value
tag2 = data >= up_value

data[tag1] = bottom_value
data[tag2] = up_value
heatMap = data

H3K27Ac_peaks = gTrack(H3K27Ac_peaks, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,colormaps = list(Group = c(enhancer = "#984EA3",promoter="#C1E5C2")),col = NA,yaxis.pretty = 1)
ATAC_peaks = gTrack(ATAC_peaks, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,colormaps = list(Group = c(enhancer = "blue",promoter="#C1E5C2")), col = NA,yaxis.pretty = 1)

#ecDNA_track = gTrack(ecDNA, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,colormaps = list(Group = c(ecDNA = "darkred")), col = NA,yaxis.pretty = 1)
cnv_track = gTrack(info, y.field = 'copyNumber',stack.gap=5,bars=TRUE, colormaps = list(Label= c("amplification" = "#E41A1C","Normal"="grey50","deletion" = "#377EB8")),y1=5,yaxis.pretty = 1,col = NA)

ATAC_track = gTrack(ATAC_node, y.field = 'score',stack.gap=5,y1=200,bars=TRUE,col = 'blue',yaxis.pretty = 1)
H3K27Ac_track = gTrack(H3K27Ac_node, y.field = 'score',stack.gap=5,y1=15,bars=TRUE,col = '#984EA3')

gene_track =  gTrack(genes, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,colormaps = list(Group = c(Gene = "blue",ecDNA="#C1E5C2")), col = NA,yaxis.pretty = 1)
SV_track = gTrack(sv_interval , edges = sv_edge , stack.gap = 5)
contact_map = gTrack(node, mdata = heatMap, stack.gap = 5,colormaps = colorRampPalette(c("white", brewer.pal(n = 9, name = "Reds")))(100))

myoutf ="/oak/stanford/groups/ccurtis2/users/ydzhao/Co_Lab/ecDNA_refine/Figures/Integration_Example/Hic_heatmap_H000530_amp2_cyc3_25kb_Simple_H3K27Ac_Count.pdf"
pdf(myoutf,width=25,height=75)
#plot(c(gTrack(info, y.field = 'copyNumber',stack.gap=5,bars=TRUE,gr.colorfield = 'Label',colormaps = list('Label' = c(amplification = "#E41A1C",normal="grey50",deletion = "#377EB8")),y1=30,yaxis.pretty = 1),gTrack(ATAC_node, y.field = 'score',stack.gap=5,y1=50,bars=TRUE,col = 'blue',yaxis.pretty = 1),gTrack(H3K27Ac_node, y.field = 'score',stack.gap=5,y1=1,bars=TRUE,col = '#984EA3',yaxis.pretty = 1),gTrack(H3K27Ac_peaks, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,col = '#984EA3',yaxis.pretty = 1),gTrack(ATAC_peaks, y.field = 'V5',stack.gap=5,y1=10,bars=TRUE,col = 'blue',yaxis.pretty = 1),gt.genes,gTrack(sv_interval , edges = sv_edge , stack.gap = 5),gTrack(node, mdata = heatMap, stack.gap = 5,colormaps = colorRampPalette(c("white", brewer.pal(n = 9, name = "Reds")))(100))),window=target_window,height=c(0.025,0.025,0.025,0.0125,0.0125,0.0125,0.035,0.35))
plot(c(cnv_track,ATAC_peaks,ATAC_track,H3K27Ac_peaks,H3K27Ac_track,gene_track,SV_track,contact_map),window=target_window,height=c(0.025,0.0125,0.025,0.0125,0.025,0.0125,0.035,0.35))
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Expression Rank RUNX2 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls())
.libPaths(c("/oak/stanford/groups/howchang/users/ydzhao/R_Library_4.2",.libPaths()[2]))
    

