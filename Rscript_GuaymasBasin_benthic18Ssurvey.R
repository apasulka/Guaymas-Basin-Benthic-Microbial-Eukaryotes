# Set working directory to contents of repository
setwd("~/Guaymas-Basin-Benthic-Microbial-Eukaryotes/")
library(vegan); library(reshape2); library(ggplot2)
#
# Download materials from github repository and supplementary from the manuscript
## Table S1 - raw OTU table
## NameSchematic.txt
# Initial OTU QC ------------------------------------------------------
count<-read.delim("TableS1.txt",sep="\t",header=TRUE)
names(count)
#
# Downstream sample analysis
colnames(count)[25]<-c("PR2") #rename taxonomy column
colsum<-apply(count[2:24],2,sum);colsum #number of sequences per sample
#
# Removed OTUs with fewer than 2 sequences total (global singletons)
count_nosingleton<-count[! rowSums(count[2:24])<2,]
colsum_nosingle<-apply(count_nosingleton[2:24],2,sum)
colsum_nosingle #can compare number of sequences lost from colsum and colsum_nosingle
lost<-colsum - colsum_nosingle
lost;hist(lost) #total number of sequences lost per sample
#
# Format count file with proper labels
head(count_nosingleton[1:2,])
df1<-melt(count_nosingleton)
# Import name schematic text to incorporate a more informative sample name.
names<-read.delim("NameSchematic.txt"); head(names)
sam_order<-as.character(names$Num)
sam_order
unique(df1$variable)
rename<-as.character(names$SampleNameR)
df1$label<-factor(df1$variable, levels = sam_order, labels = rename)
head(df1[1:2,]); names(df1)
count_nosingleton_wname<-dcast(df1[c(1,2,5,6)], OTU.ID+PR2~label)
names(count_nosingleton_wname)
head(count_nosingleton_wname[1:3,])
#
#OTU counts and information:
names(count_nosingleton_wname)
counts_only<-count_nosingleton_wname[3:25] #using singleton-less data
seq_total<-apply(counts_only,2,sum) #number of sequences per sample
OTU_count<-colSums(counts_only>0);OTU_count
OTU_single<-colSums(counts_only==1); OTU_single
OTU_double<-colSums(counts_only==2);OTU_double
OTU_true<-colSums(counts_only>2);OTU_true
#
#compile sample information
sample_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)
sample_info
# write.csv(sample_info, file="OTUstats.csv")
#
# Supplementary figures to demonstrate OTU and seq distribution among samples
head(sample_info)
sample_info$samples<-factor(row.names(sample_info), levels = rename)
allM<-melt(sample_info[c(6,2:4)])
head(allM)
#
# Figure S1
bar_stats<- ggplot(allM, aes(x=samples, y=value, fill=variable))+
  geom_bar(stat="identity",position="stack",color="black")+
  labs(title="", x="",y="Total OTUs")+theme_bw()+
  theme(axis.text.x = element_text(angle = 0,hjust=1,vjust=0.5,size=12, color="black"),axis.text.y=element_text(size=12, color="black"))+
  scale_fill_manual("",values=c("#e41a1c","#fee08b","#4393c3"),labels = c("Singletons", "Doubletons","OTUs > 2 seqs"))+
  theme(legend.position = "top")+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))+
  geom_text(x = 4, y = 220, label = "*", size=8)+
  geom_text(x = 5, y = 220, label = "*", size=8)+
  geom_text(x = 19, y = 220, label = "*", size=8)+
  geom_text(x = 20, y = 220, label = "*", size=8)+
  geom_text(x = 21, y = 220, label = "*", size=8)+
  geom_text(x = 22, y = 540, label = "x", size=4)
bar_stats
#
#
bar_seqs<-ggplot(sample_info, aes(x=samples, y=seq_total))+
  geom_bar(stat="identity", fill="#80cdc1", color="black")+
  labs(title="Total number\nof sequences per sample", x="", y="Total sequences")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=7),axis.text.y=element_text(size=12))+
  theme(legend.position = "top")
bar_seqs
#
#Binning by various taxonomic groups:
head(count_nosingleton_wname[1:2,]) #make sure taxonomy column is labeled "PR2" - done in line 39 above
length(unique(count_nosingleton_wname$OTU.ID)) #Number of distinct OTUs
library(vegan)
#
load("R_obj_rawdata.RData",verbose=T)
head(data_binned)
#
pr2_rename_taxa<-function(df){
  library(reshape2)
  split<-colsplit(df$PR2, "; ", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8","Level9", "Level10", "Level11", "Level12"))
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "No blast hit"]="No blast hit"
  split$Taxa[split$Level1 == "Unassigned"]="Unassigned"
  split$Taxa[split$Level1 == "None"]="None"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_X"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_Mikro"]="Other/unknown"
  split$Taxa[split$Level2=="Stramenopiles"]="Other Stramenopiles"
  split$Taxa[split$Level2=="Alveolata"]="Other Alveolates"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonts"
  split$Taxa[split$Level2=="Archaeplastida"]="Other Archaeplastids"
  split$Taxa[split$Level2=="Excavata"]="Excavates"
  split$Taxa[split$Level2=="Rhizaria"]="Other Rhizaria"
  split$Taxa[split$Level2=="Hacrobia"]="Other/unknown"
  split$Taxa[split$Level3=="Haptophyta"]="Haptophytes"
  split$Taxa[split$Level3=="Fungi"]="Fungi"
  split$Taxa[split$Level3=="Foraminifera"]="Foraminifera"
  split$Taxa[split$Level3=="Dinophyta"]="Dinoflagellates"
  split$Taxa[split$Level3=="Cryptophyta"]="Crytophytes"
  split$Taxa[split$Level3=="Ciliophora"]="Ciliates"
  split$Taxa[split$Level3=="Apicomplexa"]="Apicomplexa"
  split$Taxa[split$Level3=="Chlorophyta"]="Chlorophytes"
  split$Taxa[split$Level3=="Cercozoa"]="Cercozoa"
  split$Taxa[split$Level3=="Centroheliozoa"]="Centrohelida"
  split$Taxa[split$Level4=="Acantharea"]="Acantharia"
  split$Taxa[split$Level4=="Chrysophyceae-Synurophyceae"]="Other Stramenopiles"
  split$Taxa[split$Level4=="Bacillariophyta"]="Diatoms"
  split$Taxa[split$Level4=="MAST"]="MAST"
  split$Taxa[split$Level4=="Polycystinea"]="Polycystines"
  split$Taxa[split$Level4=="RAD-C"]="RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-B"]="RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-A"]="RAD (A,B,C)"
  return(split)
} #function generates Taxa column with binned taxonomic group names for PR2 assignment
binned<-pr2_rename_taxa(count_nosingleton_wname); head(binned)
data_binned<-data.frame(count_nosingleton_wname, binned) #combined back
head(data_binned[1:2,]);names(data_binned)
unique(data_binned$Taxa)
#
# Remove unwanted samples (due to low sequence number):
data.all.melt<-melt(data_binned) #melt
seq_total<-apply(data_binned[3:25],2,sum)
few1000<-names(subset(seq_total, seq_total<1000))
d<-subset(data.all.melt, !(variable %in% few1000)) #remove those
#
data_filtered <- subset(d, !grepl("Metazoa", d$PR2)) #remove OTUs IDed as metazoan
#
# save(few1000,data_filtered,data_binned, count_nosingleton_wname, sample_info, file="R_obj_rawdata.RData") #save these R objects for later!
# data_filtered - no metazoa and samples that have fewer than 1,000 sequences removed
# count_nosingleton_wname - no global singleton OTUs, with taxonomy names
# data_binned - same as 'count_nosingleton_wname', but with "binned" taxa group names
# sample_info - seq and OTU stats for each sample
# few1000 - sample names that have fewer than 100 sequences
# Whole community plots ----------------------------------------------------------------
setwd("/Users/SarahHu 1/Desktop/Projects/DeepSea_countway/R_analysis_2018/")
load(file="R_obj_rawdata.RData", verbose=T)
#
data.reform<-data_filtered
#Factoring:
tax_order<-c("Diatoms","MAST","Other Stramenopiles","Dinoflagellates","Ciliates","Apicomplexa","Other Alveolates","Cercozoa","Acantharia","RAD (A,B,C)","Polycystines","Foraminifera","Other Rhizaria","Chlorophytes","Other Archaeplastids","Centrohelida","Haptophytes","Crytophytes","Amoebozoa","Excavates","Other/unknown","Opisthokonts","Fungi","Archaea","Bacteria","No blast hit","None","Unassigned")
tax_rename<-c("Diatoms","MAST","Other Stramenopiles","Dinoflagellates","Ciliates","Apicomplexa","Other Alveolates","Cercozoans","Acantharians","RAD (A,B,C)","Polycystines","Foraminifera","Other Rhizarians","Chlorophytes","Other Archaeplastids","Centrohelids","Haptophytes","Crytophytes","Amoebozoans","Excavates","Other/unknown","Opisthokonts","Fungi","Archaea","Bacteria","No blast hit","None","Unassigned")
tax_color<-c("#800026","#cb181d","#e7298a","#df65b0","#fc4e2a","#fd8d3c","#fed976","#c7e9b4","#7fcdbb","#41ae76","#238b45","#006d2c","#00441b","#c6dbef","#6baed6","#1d91c0","#225ea8","#253494","#081d58","#54278f","#8c510a","#bf812d","#dfc27d","#bababa","#4d4d4d","#f0f0f0","#000000","white")
names(tax_color)<-tax_rename
data.reform$tax<-factor(data.reform$Taxa, levels=rev(tax_order), labels=rev(tax_rename))
head(data.reform)
# 
neworder<-c("Sample_12a_Control_1.2cm","Sample_12b_Control_1.2cm","Sample_5a_Control_0.1cm","Sample_5b_Control_0.1cm","Sample_6a_Control_1.2cm","Sample_6b_Control_1.2cm","Sample_8_Control_2.3cm","Sample_11_Edge_0.1cm","Sample_3_Edge_0.1cm","Sample_4a_Edge_1.2cm","Sample_4b_Edge_1.2cm","Sample_9_OrangeMat_0.1cm","Sample_10a_OrangeMat_1.2cm","Sample_10b_OrangeMat_1.2cm","Sample_1a_YellowMat_0.1cm","Sample_1b_YellowMat_0.1cm","Sample_2_YellowMat_1.2cm","Sample_7_YellowMat_0.1cm")
newlabel<-c("Control_1-2cm_orange","Control_1-2cm_orange","Control_0-1cm_yellow","Control_0-1cm_yellow","Control_1-2cm_yellow","Control_1-2cm_yellow","Control_2-3cm_yellow","Edge_0-1cm_orange","Edge_0-1cm_yellow","Edge_1-2cm_yellow","Edge_1-2cm_yellow","Mat_0-1cm_orange","Mat_1-2cm_orange","Mat_1-2cm_orange","Mat_0-1cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_yellow","Mat_0-1cm_yellow2")
#
data.reform$label_order<-factor(data.reform$variable, levels=neworder, labels=newlabel)
head(data.reform[1:2,])
data.reform$label_order<-as.character(data.reform$label_order)
#
#Sum taxa based on sequences:
head(data.reform)
data.agg<-aggregate(data.reform$value, by=list(Taxa=data.reform$Taxa, tax=data.reform$tax,Sample=data.reform$label_order),sum); head(data.agg) #relabeled - sum by these to sum replicates
dim(data.agg)
head(data.agg)
# Re factor to label for plot:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
data.agg$orderall<-factor(data.agg$Sample, levels=order2, labels = label2)
unique(data.agg$orderall)
#
#Figure 2 - Community composition - relative abundance of rRNA reads
plot_binned<-ggplot(data.agg[order(data.agg$tax),], aes(y=x,x=orderall,fill=tax,order=tax))+
  geom_bar(stat="identity", position="fill", color="black")+
  labs(title="PR2 database binned\nMetazoa removed", x="",y="Relative abundance of rRNA reads")+
  theme_bw()+
  theme(legend.position="right",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=9))+
  scale_fill_manual(values=tax_color,guide=guide_legend(reverse=T))+
  theme(axis.text.x = element_text(size=12,angle=0,hjust=1,vjust=0.3,color="black"),axis.text.y = element_text(size=12,color="black"))+ coord_flip()+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=0,hjust = 1,vjust = 1, color = "black"), axis.text.y=element_text(color="black"))
#
plot_binned +theme(legend.position="top")#relative abundance
# W:820 H: 600
#Get shapes for sample names/legend.
# cols.legend=c("darkblue", "lightblue", "orange", "yellow")
plot(0,1)
legend(-1,1.4,c("Mat", "Edge", "Control"), col=c("black"),pt.bg=c("yellow", "yellow", "blue"), pch=c(21,22,24), cex = 2, bty = 'n',pt.cex =5,y.intersp=2)
legend(0,1.4,c("Mat", "Edge", "Control"), col=c("black"),pt.bg=c("orange", "orange", "blue"), pch=c(21,22,24), cex = 2, bty = 'n',pt.cex =5,y.intersp=2)
#
# Now plot OTU richness for each of these groups - Figure S2
# data_filtered - no metazoa and samples that have fewer than 1,000 sequences removed
head(data_filtered) # already melted
tmp0<-data_filtered
tmp0$bin<-ifelse(tmp0$value > 0, 1, 0)
tmp1<-tmp0[c(1,3:10,15:16,18)]
head(tmp1)
tmp2<-aggregate(tmp1$bin, by=list(variable=tmp1$variable, Taxa=tmp1$Taxa),sum)
head(tmp2) 
bin_taxOTUs<-tmp2
#Factoring:
tax_order<-c("Diatoms","MAST","Other Stramenopiles","Dinoflagellates","Ciliates","Apicomplexa","Other Alveolates","Cercozoa","Acantharia","RAD (A,B,C)","Polycystines","Foraminifera","Other Rhizaria","Chlorophytes","Other Archaeplastids","Centrohelida","Haptophytes","Crytophytes","Amoebozoa","Excavates","Other/unknown","Opisthokonts","Fungi","Archaea","Bacteria","No blast hit","None","Unassigned")
tax_rename<-c("Diatoms","MAST","Other Stramenopiles","Dinoflagellates","Ciliates","Apicomplexa","Other Alveolates","Cercozoans","Acantharians","RAD (A,B,C)","Polycystines","Foraminifera","Other Rhizarians","Chlorophytes","Other Archaeplastids","Centrohelids","Haptophytes","Crytophytes","Amoebozoans","Excavates","Other/unknown","Opisthokonts","Fungi","Archaea","Bacteria","No blast hit","None","Unassigned")
tax_color<-c("#800026","#cb181d","#e7298a","#df65b0","#fc4e2a","#fd8d3c","#fed976","#c7e9b4","#7fcdbb","#41ae76","#238b45","#006d2c","#00441b","#c6dbef","#6baed6","#1d91c0","#225ea8","#253494","#081d58","#54278f","#8c510a","#bf812d","#dfc27d","#bababa","#4d4d4d","#f0f0f0","#000000","white")
names(tax_color)<-tax_rename
bin_taxOTUs$tax<-factor(bin_taxOTUs$Taxa, levels=rev(tax_order), labels = rev(tax_rename))
# 
neworder<-c("Sample_12a_Control_1.2cm","Sample_12b_Control_1.2cm","Sample_5a_Control_0.1cm","Sample_5b_Control_0.1cm","Sample_6a_Control_1.2cm","Sample_6b_Control_1.2cm","Sample_8_Control_2.3cm","Sample_11_Edge_0.1cm","Sample_3_Edge_0.1cm","Sample_4a_Edge_1.2cm","Sample_4b_Edge_1.2cm","Sample_9_OrangeMat_0.1cm","Sample_10a_OrangeMat_1.2cm","Sample_10b_OrangeMat_1.2cm","Sample_1a_YellowMat_0.1cm","Sample_1b_YellowMat_0.1cm","Sample_2_YellowMat_1.2cm","Sample_7_YellowMat_0.1cm")
newlabel<-c("Control_1-2cm_orange","Control_1-2cm_orange","Control_0-1cm_yellow","Control_0-1cm_yellow","Control_1-2cm_yellow","Control_1-2cm_yellow","Control_2-3cm_yellow","Edge_0-1cm_orange","Edge_0-1cm_yellow","Edge_1-2cm_yellow","Edge_1-2cm_yellow","Mat_0-1cm_orange","Mat_1-2cm_orange","Mat_1-2cm_orange","Mat_0-1cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_yellow","Mat_0-1cm_yellow2")
bin_taxOTUs$label_order<-factor(bin_taxOTUs$variable, levels=neworder, labels=newlabel)
bin_taxOTUs$label_order<-as.character(bin_taxOTUs$label_order)
#
# Now add replicates:
head(bin_taxOTUs)
bin_taxOTUs.agg<-aggregate(bin_taxOTUs$x, by=list(Taxa=bin_taxOTUs$Taxa, tax=bin_taxOTUs$tax,label_order=bin_taxOTUs$label_order),sum)
head(bin_taxOTUs.agg) #relabeled - sum by these to sum replicates
#
# Re factor to label for plot:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
head(bin_taxOTUs.agg)
bin_taxOTUs.agg$orderall<-factor(bin_taxOTUs.agg$label_order, levels=order2, labels = label2)
#plot
plot_OTUrich<-ggplot(bin_taxOTUs.agg[order(bin_taxOTUs.agg$tax),], aes(y=x,x=orderall,fill=tax,order=tax))+
  geom_bar(stat="identity", position="fill", color="black")+
  # geom_bar(stat="identity", position="stack", color="black")+ #optional
  labs(title="OTU richness\nMetazoa removed", x="",y="Number of OTUs")+
  theme_bw()+theme(legend.position="right",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=9))+ 
  scale_fill_manual(values=tax_color,guide=guide_legend(reverse=T))+
  theme(axis.text.x = element_text(size=12,angle=0,hjust=1,vjust=0.3,color="black"),axis.text.y = element_text(size=12,color="black"))+ coord_flip()+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=0,hjust = 1,vjust = 1, color = "black"), axis.text.y=element_text(color="black"))
#
plot_OTUrich +theme(legend.position="top")
plot_OTUrich %+% subset(bin_taxOTUs.agg, !(tax %in% "Unassigned"))+theme(legend.position="top")
# W:820 H: 600
#
save(bin_taxOTUs, bin_taxOTUs.agg, data.reform, data.agg, file="Data_to_plot.RData")
#
# updated SHu 06-01-2018

# Composition of ciliates --------------------------------------------------------
load("Data_to_plot.RData", verbose=TRUE)
head(data.reform[1:2,])
data_lev4<-aggregate(data.reform$value, by=list(Level4=data.reform$Level4 ,Taxa=data.reform$Taxa, tax=data.reform$tax,label_order=data.reform$label_order),sum)
head(data_lev4)
unique(data_lev4$label_order)
# re-factor:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
data_lev4$orderall<-factor(data_lev4$label_order, levels=order2, labels = label2)
unique(data_lev4$orderall)
head(data_lev4[1:3,])
path_color<-c("#b2df8a","#33a02c","#a6cee3","#1f78b4","#e31a1c","#fed976","#e7298a","#c7e9b4","#f03b20","#fd8d3c","#662506","#9e9ac8","#54278f","#01665e","#80cdc1","#ffffff","#878787","#000000")
#
plot_class<-ggplot(data_lev4, aes(y=x,x=orderall,fill=Level4))+
  geom_bar(stat="identity", position="fill", color="black")+
  labs(title="PR2 database binned", x="",y="Relative abundance of rRNA reads")+theme_bw()+
  coord_flip()+
  theme(legend.position="top",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=12), axis.text.x = element_text(color="black", size=12),axis.text.y = element_text(color="black", size=12))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),strip.background = element_blank(), axis.text.x = element_text(angle=0,hjust = 1,vjust = 1, color = "black"), axis.text.y=element_text(color="black"))
#
plot_class %+% subset(data_lev4, Taxa %in% "Ciliates") +labs(title="Ciliates only")+scale_fill_manual(values=rev(path_color),guide=guide_legend(reverse=T))+theme(legend.position = "top")
# W:820 H: 600
#
# Generate table to compile ciliate numbers
data_lev5<-aggregate(data.reform$value, by=list(Level4=data.reform$Level4, Level5=data.reform$Level5,Taxa=data.reform$Taxa, tax=data.reform$tax,label_order=data.reform$label_order),sum)
head(data_lev5)
ciliates_lev5<-subset(data_lev5, Taxa %in% "Ciliates" & x > 0)
dim(ciliates_lev5)
write.csv(ciliates_lev5, file="Ciliate_seq_count.csv")
#
#

# updated SHu 06-01-2018
#
data.class<-aggregate(data.reform$value, by=list(Level4=data.reform$Level4 ,Taxa=data.reform$Taxa, tax=data.reform$tax,label_order=data.reform$label_order),sum); head(data.agg)
#Labels=data.reform$Labels #relabeled - sum by these to sum replicates
#first bin by taxa category previously given and then plot out the class level. Can change for another level as well.
unique(data.class$label_order) #use to order by habitat:
#Re-factor
all_samples<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
# bymat<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_0-1cm_orange","Control_1-2cm_orange")
data.class$habitat<-factor(data.class$label_order, levels=all_samples)
# data.class$habitat<-factor(data.class$label_order, levels=bymat)
unique(data.class$Level4)
unique(data.class$Taxa)
path_color<-c("#b2df8a","#33a02c","#a6cee3","#1f78b4","#e31a1c","#fed976","#e7298a","#c7e9b4","#f03b20","#fd8d3c","#662506","#9e9ac8","#54278f","#01665e","#80cdc1","#ffffff","#878787","#000000")
plot_class<-ggplot(data.class, aes(y=x,x=habitat,fill=Level4))+
  geom_bar(stat="identity", position="fill", color="black")+labs(title="PR2 database binned", x="",y="Relative abundance of rRNA reads")+
  theme_bw()+ 
  coord_flip()+
  theme(legend.position="top",legend.title = element_blank(),plot.title=element_text(hjust = 0,face='bold',size=12), axis.text.x = element_text(color="black", size=12),axis.text.y = element_text(color="black", size=12))
#
plot_class %+% subset(data.class, Taxa %in% "Ciliates") +labs(title="Ciliates\nonly")+scale_fill_manual(values=rev(path_color))+theme(legend.position = "none")
#
#Re-label for ciliates only!
all_samples_label<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")

plot_class %+% subset(data.class, Taxa %in% "Ciliates") +theme(legend.position="none")+scale_x_discrete(labels=as.character(all_samples_label))+scale_fill_manual(values=rev(path_color))
#Saved as svg: W:870, H:590
#

# Bubble plots - 3 groups ------------------------------------------------
load("Data_to_plot.RData", verbose=TRUE)
head(data.reform[1:2,])
names(data.reform)
#
library(dplyr)
data.relAbun<-data.reform %>% 
  group_by(variable) %>%
  mutate(RelAbun=value/sum(value)) %>%
  data.frame
head(data.relAbun)
# Computed relative abundance for the entire dataset. An alternative would be to calculate relative abundance of classes for each group individually.
#
Rhiz<-subset(data.relAbun, Level2 %in% "Rhizaria")
Cil<-subset(data.relAbun,Level3 %in% "Ciliophora")
Api<-subset(data.relAbun,Level4 %in% "Apicomplexa_X")
#
head(Rhiz[1:2,])
# Overall binning, Labels= new labels!
rhiz<-aggregate(Rhiz$RelAbun, by=list(Taxa=Rhiz$Level2, Labels=Rhiz$label_order,Group=Rhiz$Level4),sum) #sum at level 4
api<-aggregate(Api$RelAbun, by=list(Taxa=Api$Level4, Labels=Api$label_order,Group=Api$Level5),sum)
cil<-aggregate(Cil$RelAbun, by=list(Taxa=Cil$Level3, Labels=Cil$label_order,Group=Cil$Level4),sum) #use below!
combine<-rbind(rhiz,api,cil)
head(combine)
# hist(log(combine$x))
grab<-c("Ciliophora","Apicomplexa_X", "Rhizaria")
combine$order<-factor(combine$Taxa, levels=grab, label=c("Ciliates", "Apicomplexa", "Rhizaria"))
#
# Now we have relative abundances for subsets of families within each of the taxonomic groups.
head(combine[1:3,])
#Relabel:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
combine$Label_final<-factor(combine$Labels, levels=rev(order2), labels = rev(label2))
head(combine[1:4,])
#
label_col=c("yellow","yellow","yellow","orange","orange","yellow","yellow","orange","blue","blue","blue","blue")
length(label_col)
names(label_col)=label2
unique(combine$Group)
api_remove<-c("Novel-Apicomplexa-Class-1","Novel-Apicomplexa-Class-2","Coccidia", "Apicomplexa_XX","") #keep:gregarines and colpodellidae
dim(combine)
combine <- subset(combine, !(Group %in% api_remove) & x>0)
unique(combine$Taxa)
head(combine)
#
plotbubble<-ggplot(combine[order(combine$order),], aes(x=Label_final, y=Group,order=order))+
  geom_point(shape=21,aes(fill=Label_final,size=x),color="black")+
  theme_bw()+ 
  scale_size(range = c(1,9))+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5, color="black"), axis.text.y = element_text(color="black"))+facet_grid(order~.,scales="free", space="free")+
  scale_fill_manual(values=label_col) +
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "black"))+labs(x="", y="", title="")
plotbubble
#H:900, W:650
#
#last updated 06-12-2018
#
# Presence-absence UpsetR -------------------------------------------------
setwd("/Users/SarahHu 1/Desktop/Projects/DeepSea_countway/R_analysis_2018/")
load("Data_to_plot.RData", verbose=TRUE)
head(data.reform[1:2,])
#
tmp0<-data.reform
unique(tmp0$label_order)
names(tmp0)
split<-colsplit(tmp0$label_order, "_", c("Habitat", "Depth", "Color"))
head(split)
split$habitat_wColor<-paste(split$Habitat, split$Color, sep="_")
tmp1<-data.frame(tmp0,split)
head(tmp1)
#
# Sum, change to binary, & cast
Habitat_wCol<-aggregate(tmp1$value, by=list(OTU.ID=tmp1$OTU.ID, habitat_wColor=tmp1$habitat_wColor), sum)
Habitat_wCol$bin<-ifelse(Habitat_wCol$x > 0, 1, 0)
head(Habitat_wCol)
wCol_bin<-dcast(Habitat_wCsol[c(1:2,4)], OTU.ID~habitat_wColor, fill=0)
head(wCol_bin)
###
Habitat_only<-aggregate(tmp1$value, by=list(OTU.ID=tmp1$OTU.ID, Habitat=tmp1$Habitat), sum)
Habitat_only$bin<-ifelse(Habitat_only$x > 0, 1, 0)
head(Habitat_only)
HabOnly_bin<-dcast(Habitat_only[c(1:2,4)], OTU.ID~Habitat, fill=0)
head(HabOnly_bin)
#
# Ciliates only:
tmp2<-subset(tmp1, Taxa %in% "Ciliates")
head(tmp2)
#
Habitat_Ciliateonly<-aggregate(tmp2$value, by=list(OTU.ID=tmp2$OTU.ID, Habitat=tmp2$Habitat), sum)
Habitat_Ciliateonly$bin<-ifelse(Habitat_Ciliateonly$x > 0, 1, 0)
head(Habitat_Ciliateonly)
Habitat_CiliateBin<-dcast(Habitat_Ciliateonly[c(1:2,4)], OTU.ID~Habitat, fill=0)
head(Habitat_CiliateBin)
#
library(UpSetR)
#Look at all
upset(wCol_bin,keep.order=TRUE,order.by=c("freq"))
upset(HabOnly_bin,keep.order=TRUE,order.by=c("freq"))
upset(Habitat_CiliateBin,keep.order=TRUE,order.by=c("freq"))
#
# Last updated 04-27-2018

# OTU richness - Ciliates only --------------------------------------------
# Get OTU richness for ciliates
load("Data_to_plot.RData", verbose=TRUE)
head(data.reform[1:2,])
Cil<-subset(data.reform,Level3 %in% "Ciliophora")
ciliate<-Cil
ciliate$value[ciliate$value !=0 ]<- 1 #change all non zeros in "x" col to 1
head(ciliate)
#
otu_rich<-aggregate(ciliate$value, by=list(Taxa=ciliate$Level3, Labels=ciliate$label_order, Group=ciliate$Level4),sum)
head(otu_rich)
#
#Relabel:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
otu_rich$Label_final<-factor(otu_rich$Labels, levels=rev(order2), labels = rev(label2))
head(otu_rich)
label_col=c("yellow","yellow","yellow","orange","orange","yellow","yellow","orange","blue","blue","blue","blue")
length(label_col)
names(label_col)=label2
head(otu_rich)
#
plotbubble_rich<-ggplot(otu_rich, aes(x=Label_final,y=Group,fill=Label_final,size=x))+
  geom_point(stat="identity",shape=21,color="black")+
  theme_bw()+ 
  scale_size(range = c(3,20))+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5, color="black"), axis.text.y = element_text(color="black"))+
  scale_fill_manual(values=label_col) +
  theme(strip.background = element_blank(),panel.border = element_rect(colour = "black"))+
  labs(x="", y="", title="Ciliate OTU richness")
plotbubble_rich %+% subset(otu_rich, !(Group %in% "" | x==0))
#H:880, W-1000
#
# import OTU distribution information, so pie charts can represent how many of those OTUs are shared vs. unique.
cil_dist<-read.delim("CiliateOTUs_distribution_wtax.txt")
head(otu_rich)
Cil<-subset(data.reform,Level3 %in% "Ciliophora")
ciliate_wcat<-join(Cil, cil_dist, type="left", match="first")
head(ciliate_wcat)
ciliate_wcat$value[ciliate_wcat$value !=0 ]<- 1 #change all non zeros in "x" col to 1
#
otu_rich_wcat<-aggregate(ciliate_wcat$value, by=list(Taxa=ciliate_wcat$Level3, Labels=ciliate_wcat$label_order, Group=ciliate_wcat$Level4, category=ciliate_wcat$Category_binary),sum)
head(otu_rich_wcat)
#
#Relabel:
order2<-c("Mat_0-1cm_yellow2","Mat_1-2cm_yellow","Mat_0-1cm_yellow","Mat_1-2cm_orange","Mat_0-1cm_orange","Edge_1-2cm_yellow","Edge_0-1cm_yellow","Edge_0-1cm_orange","Control_2-3cm_yellow","Control_1-2cm_yellow","Control_0-1cm_yellow","Control_1-2cm_orange")
label2<-c("7 Yellow Mat 0-1cm","2 Yellow Mat 1-2cm","1 Yellow Mat 0-1cm*","10 Orange Mat 1-2cm*","9 Orange Mat 0-1cm","4 Edge 1-2cm*","3 Edge 0-1cm","11 Edge 0-1cm","8 Control 2-3cm","6 Control 1-2cm*","5 Control 0-1cm*","12 Control 1-2cm*")
otu_rich_wcat$Label_final<-factor(otu_rich_wcat$Labels, levels=rev(order2), labels = rev(label2))
head(otu_rich_wcat)
label_col=c("yellow","yellow","yellow","orange","orange","yellow","yellow","orange","blue","blue","blue","blue")
length(label_col)
names(label_col)=label2
head(otu_rich_wcat)
# otu_rich_wcat_no0<-subset(otu_rich_wcat, x>0)
#
# save(otu_rich_wcat, otu_rich, file="plotBubble_Ciliatonly.RData")
load("plotBubble_Ciliatonly.RData")
piebubble<-ggplot(otu_rich_wcat, aes(x=Group, y=x, fill=category))+geom_bar(stat="identity",color=NA,position="fill")+coord_polar(theta='y')+facet_wrap(~Label_final,ncol=12)+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=3),panel.grid.minor = element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), strip.background = element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(), axis.line = element_blank(), strip.text = element_text(size=3))+scale_fill_manual(values=c("#737373", NA))+labs(x="", y="")+ theme(panel.spacing = unit(0, "lines"))
#+theme(legend.position = "none")+theme(plot.margin = unit(c(0,0,0,0), "cm"))
#
piebubble %+% subset(otu_rich_wcat, (Group %in% "Armophorea"))
#+theme(legend.position = "top")
##
#
# plotbubble_rich %+% subset(otu_rich, !(Group %in% "" | x==0))
#
unique(otu_rich_wcat$Group)
library(cowplot)
# plot_grid(piebubble %+% subset(otu_rich_wcat, (Group %in% "Phyllopharyngea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Plagiopylea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Prostomatea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Spirotrichea")), align="hv", ncol=1)
#
# split into several:
unique(otu_rich_wcat$Group)
#
svg(filename="pies_part1.svg", width=10, height=8)
plot_grid(piebubble %+% subset(otu_rich_wcat, (Group %in% "Armophorea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Cariacotrichea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Ciliophora-10")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Ciliophora-5")), align="hv", ncol=1)
dev.off()
#
svg(filename="pies_part2.svg", width=10, height=8)
plot_grid(piebubble %+% subset(otu_rich_wcat, (Group %in% "Ciliophora-7")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Ciliophora-8")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Colpodea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Heterotrichea")), align="hv", ncol=1)
dev.off()
#
svg(filename="pies_part3.svg", width=10, height=8)
plot_grid(piebubble %+% subset(otu_rich_wcat, (Group %in% "Karyorelictea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Litostomatea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Oligohymenophorea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Phyllopharyngea")), align="hv", ncol=1)
dev.off()
#
svg(filename="pies_part4.svg", width=10, height=8)
plot_grid(piebubble %+% subset(otu_rich_wcat, (Group %in% "Plagiopylea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Prostomatea")), piebubble %+% subset(otu_rich_wcat, (Group %in% "Spirotrichea")), align="hv", ncol=1)
dev.off()
#
#
# Last updated SHu - 05/2/2018

# Beta diversity metrics - MDS and ANOSIMS --------------------------------
load(file="R_obj_rawdata.RData", verbose=T)
head(data_filtered[1:2,])

# Convert the raw count data (without metazoa or low sequence smaples) to proportions within a sample
names(data_filtered)
data_tmp<-data_filtered

# Cast wide, only need count information and OTU.ID
names(data_tmp)
data_tmp2<-dcast(data_tmp[c(1,16:17)], OTU.ID~variable)
row.names(data_tmp2)<-data_tmp2$OTU.ID; data_tmp2$OTU.ID<-NULL # set row names to OTU IDs and remove original OTU.ID column

#Transpose and remove the first row
data_tmp3<-as.data.frame(t(data_tmp2))
data_tmp4<-data_tmp3[-1,]
head(data_tmp4[1:2])

# Convert dataframe to numeric
indx <- sapply(data_tmp4, is.factor)
data_tmp4[indx] <- lapply(data_tmp4[indx], function(x) as.numeric(as.character(x)))

#Create relative abundance matrix
data_prop <- data_tmp4/(rowSums(data_tmp4)) #divide each cell by rowSum
head(data_prop[1:2])

#Transform data to test best fit for input data
#4th root transformation
data_prop_4th <- as.data.frame(data_prop^(1/4))

#Sqrt transformation
data_prop_sqrt <- sqrt(data_prop)

#Presence absence
data_prop_pres<-as.matrix(data_prop)
data_prop_pres<-1*(data_prop_pres>0)
data_prop_pres <-(as.data.frame(data_prop_pres))

# NMDS -- Calculate NMDS for each one
NMDS_raw=metaMDS(data_prop,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
head(NMDS_raw$points) #given points for MDS plot
raw_points <- NMDS_raw$points[1:nrow(NMDS_raw$points),] #extract points

NMDS_sqrt=metaMDS(data_prop_sqrt,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
sqrt_points <-NMDS_sqrt$points[1:nrow(NMDS_sqrt$points),]

NMDS_4th=metaMDS(data_prop_4th,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
fourth_points <-NMDS_4th$points[1:nrow(NMDS_4th$points),]

NMDS_pres=metaMDS(data_prop_pres,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
pres_points <-NMDS_pres$points[1:nrow(NMDS_pres$points),]

Stress <- c(NMDS_raw$stress, NMDS_raw$stress, NMDS_sqrt$stress, NMDS_pres$stress) #obtain stress values for each metaMDS analysis
names(Stress) <- c("NMDS_raw","NMDS_4th", "NMDS_sqrt", "NMDS_pres")
print(Stress) #Shows the stress of each NMDS

# write.csv(Stress, file="MDS_Stress_nosingletons.csv")
# Checkpoint to save these dataframes:
save(data_prop, data_prop_4th, data_prop_pres, data_prop_sqrt,raw_points, sqrt_points, fourth_points, pres_points, file="Normed_df_objs.RData")
# Load from previous checkpoint:
# load("Normed_df_objs.RData", verbose=T)

Points <-cbind(raw_points, sqrt_points, fourth_points, pres_points) #gather points together
colnames(Points) <- c("MDS1_raw","MDS2_raw", "MDS1_sqrt","MDS2_sqrt","MDS1_4th","MDS2_4th", "MDS1_pres", "MDS2_pres")
Points <- as.data.frame(Points)
Points$Sample_Name <- row.names(Points)
head(Points) #save points and stresses for each
head(Stress)

# Merge Relative Abundances, Metadata and MDS results
unique(Points$Sample_Name) # Make sure sample names in Points matches samples names in metadata sheet
#
meta <-read.csv("meta_Vent.csv", header=TRUE) #SH
#Blanks in data (which included not detected and not available) were left blank. R filled in as "NA"

head(meta)
data_prop$Sample_Name <- row.names(data_prop)
data_prop_4th$Sample_Name <- row.names(data_prop_4th)

# Merge with meta data
merged_data <- merge(meta, data_prop, by = "Sample_Name", incomparables = NULL) #this example merges the untransformed relative abundance with the meta data - you could merge with any of the transformed data as well. 
merged_data_4th <- merge(meta, data_prop_4th, by = "Sample_Name", incomparables = NULL)

head(merged_data[1:2])
#merge with "Point"
Points$Sample_Name <- row.names(Points)
merged_data_MDS <- merge(merged_data, Points, by = "Sample_Name")
merged_data_4th_MDS <- merge(merged_data_4th, Points, by = "Sample_Name")

#Checkpoint
save(merged_data_4th_MDS, merged_data_MDS, file="Merged_dfs.RData")
load("Merged_dfs.RData", verbose=T)
# write.csv(merged_data_MDS,file = paste ("iTag_18S_MergedData_", Sys.Date(), ".csv", sep=""))
# write.csv(merged_data_4th_MDS,file = paste ("iTag_18S_MergedData_4th_", Sys.Date(), ".csv", sep=""))

#NMDS Only Active Samples Only - 4th root transform
active <- subset(merged_data_MDS, Habitat == "Edge" | Habitat == "Mat") #subset from the merged data. Re-perform MDS with only "active" samples
# head(merged_data_MDS[1:8])
NMDS_4th_active=metaMDS(as.data.frame((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4)),distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
fourth_points_active <-as.data.frame(NMDS_4th_active$points[1:nrow(NMDS_4th_active$points),])
stress_active <- NMDS_4th_active$stress
# write.csv(stress_active, file="MDS_Stress_activeOnly.csv")
# write.csv(fourth_points_active, file="fourth_points_active.csv")

## Last updated 03-28-2018 - SHu

#  NMDS FIGURE -----------------------------------------------------------
# save(merged_data_MDS,file="df_forMDSplot.RData")
load("df_forMDSplot.RData")
head(merged_data_MDS[1:10])
detail<-unique(merged_data_MDS$Detail);detail
length(unique(merged_data_MDS$Detail))
color_detail<-c("orange", "orange", "blue", "yellow", "yellow")
names(color_detail)<-detail
colScale<-scale_color_manual(values=color_detail)
merged_data_MDS$hab<-factor(merged_data_MDS$Detail,levels=detail)
shape_details<-c(21, 22, 24, 21, 22)
head(merged_data_MDS[1:11])
head(merged_data_MDS$Sample_ID)
#
c_wcol <- ggplot(merged_data_MDS, aes(x = MDS1_4th, y = MDS2_4th, shape=hab,fill = hab, label=New_Label), color="black") + geom_text(size=4,vjust=2,hjust="inward") + geom_point(size=6) + theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +scale_fill_manual(values=color_detail)+scale_shape_manual(values = shape_details)
c_wcol+theme(strip.background = element_blank())
#vjust="inward",hjust="inward"
#W:720, H:520
merged_data_MDS[c(1,5,11)]
# last updated - SHu - 03-28-2018




# ANOSIM ------------------------------------------------------------------
# Subsetting merged data for OTU columns of the merged spreadsheet (can also use the unmerged data 'data_prop')
load("Merged_dfs.RData", verbose =T)
bio_data <- merged_data_MDS[,grep("^denovo", names(merged_data_MDS), value=TRUE)]
bio_data_4th <- merged_data_4th_MDS[,grep("^denovo", names(merged_data_4th_MDS), value=TRUE)]

# Calculate ANOSIM
ANOSIM_Habitat_4th=with(merged_data_4th_MDS,anosim(bio_data_4th,Habitat,permutations=999,distance="bray"))
#This example uses 'habitat' as the factor, but you can test anything as long as there are more than two samples per factor and you have that meta-data for all samples.
#Can do this on non-transformed data as well
#ANOSIM_Habitat=with(merged_data_MDS,anosim(bio_data,Habitat,permutations=999,distance="bray"))

#see ANOSIM results
summary(ANOSIM_Habitat_4th)
#Save ANOSIM Output to working directory
sink("summary_ANOSIM_Habitat_4th_nosingletons.txt")
summary(ANOSIM_Habitat_4th)
sink()

#Additional ANOSIMS
#compare ANOSIM R stat 
ANOSIM_Horizon_4th=with(merged_data_4th_MDS,anosim(bio_data_4th,Horizon,permutations=999,distance="bray")) #by habitat
summary(ANOSIM_Horizon_4th)
sink("summary_ANOSIM_Horizon_4th_nosingleton.txt")
summary(ANOSIM_Horizon_4th)
sink()

ANOSIM_MatColor_4th=with(merged_data_4th_MDS,anosim(bio_data_4th,MatColor,permutations=999,distance="bray")) #by mat color
sink("summary_ANOSIM_MatColor_4th_nosingleton.txt")
summary(ANOSIM_MatColor_4th)
sink()


#Subset ANOSIMS - active only (remove habitat = control)
ANOSIM_MatColor_activeOnly_4th=with(active,anosim(as.data.frame((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4)),MatColor,permutations=999,distance="bray"))
sink("summary_ANOSIM_MatColor_activeOnly_4th_nosingleton.txt")
summary(ANOSIM_MatColor_activeOnly_4th)
sink()

ANOSIM_Habitat_activeOnly_4th=with(active,anosim(as.data.frame((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4)),Habitat,permutations=999,distance="bray"))
sink("summary_ANOSIM_Habitat_activeOnly_4th_nosingleton.txt")
summary(ANOSIM_Habitat_activeOnly_4th)
sink()

ANOSIM_Horizon_activeOnly_4th=with(active,anosim(as.data.frame((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4)),Horizon,permutations=999,distance="bray"))
sink("summary_ANOSIM_Horizon_activeOnly_4th_nosingleton.txt")
summary(ANOSIM_Horizon_activeOnly_4th)
sink()


# SIMPER ------------------------------------------
load("Merged_dfs.RData", verbose = T)
#calculate SIMPER according to metadata
SIMPER_Habitat_4th=with(merged_data_MDS,simper(bio_data_4th,Habitat))

#output SIMPER results to .txt file - the output is comparison between each site
options(max.print=999999999)
write.table(capture.output(summary(SIMPER_Habitat_4th)),"SIMPER_Habitat_4th_nosingleton.txt")
# head(summary(SIMPER_Habitat_4th$Mat_Control))

#Additional SIMPER
SIMPER_MatColor_4th=with(merged_data_MDS,simper(bio_data_4th,MatColor))
options(max.print=999999999)
write.table(capture.output(summary(SIMPER_MatColor_4th)),"SIMPER_MatColor_4th.txt")

#Additional SIMPER - Mat color with Active Only
SIMPER_MatColor_activeOnly_4th=with(active,simper((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4),MatColor))
options(max.print=999999999)
write.table(capture.output(summary(SIMPER_MatColor_activeOnly_4th)),"SIMPER_MatColor_activeOnly_4th.txt")

#Additional SIMPER - Sediment Horizon with Active Only
SIMPER_Horizon_activeOnly_4th=with(active,simper((active[,grep("^denovo", names(merged_data_MDS), value=TRUE)])^(1/4),Horizon))
options(max.print=999999999)
write.table(capture.output(summary(SIMPER_Horizon_activeOnly_4th)),"SIMPER_Horizon_activeOnly_4th.txt")

# Last updated- 03-28-2018 SHu


# Alpha diversity ---------------------------------------------------------
load(file="R_obj_rawdata.RData", verbose=T)
head(data_binned[1:2,])
names(data_binned)
# Cast wide and set row.names equal to OTU IDs
data_tmp<-data_binned[c(1,3:25)]
names(data_tmp)  
row.names(data_tmp)<-data_tmp$OTU.ID; data_tmp$OTU.ID<-NULL

# Remove samples with too few sequences
alpha_input<-data_tmp[ ,!(names(data_tmp) %in% few1000)]
dim(data_tmp); dim(alpha_input); length(few1000)

# Find sample with the fewest sequences and randomly subsample
sub<-min(colSums(alpha_input))
sub
rare <- rrarefy(t(alpha_input), sub) #rarefy - subsample to equal depth
subsampled<-as.data.frame(t(rare))
colSums(subsampled)
shannon<-diversity(subsampled,index="shannon",2)
shannon #diversity measurement that accounts for both abundnace and evenness of the species present (both evenness and richness). proportion of species relative to the total multiplied by the ln of the proportion
invsimp<-diversity(subsampled,index="invsimpson",2)
invsimp #evennes in the community
OTU_count<-colSums(subsampled>0)
OTU_count
alpha<-data.frame(shannon, invsimp,OTU_count)
alpha$samples<-row.names(alpha)
alpha.m<-melt(alpha)
head(alpha.m)
var<-colsplit(alpha.m$samples, "_", c("Sample","Num", "Habitat", "Sed_horiz"))
alpha_div<-data.frame(var,alpha.m)
head(alpha_div)
unique(alpha_div$variable)

layout(rbind(c(1:3)))
cols<-c("darkblue", "lightblue", "orange", "yellow")
boxplot(alpha_div$value[alpha_div$variable=="shannon"]~alpha_div$Habitat[alpha_div$variable=="shannon"],main="Shannon diversity index",xaxt='n',col=cols)
boxplot(alpha_div$value[alpha_div$variable=="invsimp"]~alpha_div$Habitat[alpha_div$variable=="invsimp"],main="Inverse Simpsons diversity index",xaxt='n',col=cols)
legend(.9,38,c("Control", "Edge of mat", "Orange mat", "Yellow mat"),col=cols, pch=15,cex=1.5,bty = "n")
boxplot(alpha_div$value[alpha_div$variable=="OTU_count"]~alpha_div$Habitat[alpha_div$variable=="OTU_count"],main="OTU richness",xaxt='n',col=cols)
#
#
#
# Rarefaction curve
head(rare[1:3])
rare<-as.data.frame(rare)
# order.2<-c("X1","X2","X3","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X22","X23") #without low seq samples
# rare$samples<-row.names(rare)
# rare$samples<-factor(rare$samples, levels = order.2)
# rename.2<-c("1_GY08-009_YellowMat_0-1cm","2_GY08-010_YellowMat_0-1cm","3_GY08-012_YellowMat_1-2cm","6_GY08-025_Edgeofmat_0-1cm","7_GY08-027_Edgeofmat_1-2cm","8_GY08-028_Edgeofmat_1-2cm","9_GY08-039_Control-coldmud_0-1cm","10_GY08-040_Control-coldmud_0-1cm","11_GY08-042_Control-coldmud_1-2cm","12_GY08-043_Control-coldmud_1-2cm","13_GY08-054_YellowMat_0-1cm","14_GY08-078_Control-coldmud_2-3cm","15_GY08-212_orangemat_0-1cm","16_GY08-214_orangemat_1-2cm","17_GY08-215_orangemat_1-2cm","18_GY08-221_Edgeofmat_0-1cm","22_GY08-232_Control-coldmud_1-2cm","23_GY08-233_Control-coldmud_1-2cm") #without low seq samples
# rare$label<-factor(rare$samples, labels=rename.2)
# row.names(rare)<-rare$label
dim(rare)
# head(rare[6952:6954])
# rare<-rare[1:6953]
head(rare[1:3])
unique(row.names(rare))
# Place colors in order of the row.names above
col <- c("yellow","yellow","yellow","lightblue","lightblue","lightblue","darkblue","darkblue","darkblue","darkblue","yellow","darkblue","orange","orange","orange","lightblue","darkblue","darkblue")
cols.legend=c("darkblue", "lightblue", "orange", "yellow")
#already randomly subsampled to the same number of sequences:
curves <- rarecurve(rare, step = 20, sample = sub, lwd=2,col=col,lty = "solid", label = FALSE, ylab= "OTUs") #Change to label=TRUE to label actual curves (also check colors/labels)
legend(0,840,c("Control", "Edge of mat", "Orange mat", "Yellow mat"),col=cols.legend, pch=15, cex = 1, bty = 'n') #W:600, H:540
#
cols<-c("darkblue", "lightblue", "orange", "yellow")
boxplot(alpha_div$value[alpha_div$variable=="shannon"]~alpha_div$Habitat[alpha_div$variable=="shannon"],main="Shannon diversity index",xaxt='n',col=cols) #W:500, H:540
#
# Last updated SHu - 08-23-2018


# EXCESS CODE:

