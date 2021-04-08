library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

metric <- args[1]
order <- args[2]
afd_file <- args[3]
metric_file <- args[4]
gff_file <- args[5]
output <- args[6]
chr <- args[7]
window_start <- as.integer(args[8])
window_end <- as.integer(args[9])
plot_start <- as.integer(args[10])
plot_end <- as.integer(args[11])

afd_table=read.table(afd_file,header=T)
metric_table=read.table(metric_file,header=T)
gff_table=read.table(gff_file,header=T)
#metric_table$middle <- ((metric_table$end-metric_table$start)/2)+metric_table$start
png(output, width = 10, height = 5, units = 'cm', res = 300)

genes=subset(gff_table,(gff_table[,2]=='gene'))
genes[genes$orientation == "-", c('start', 'end')] <- genes[genes$orientation == "-", c('end', 'start')]
exons=subset(gff_table,(gff_table[,2]=='exon'))

p <- ggplot()+
  annotate("rect", xmin=window_start, xmax=window_end, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey")+
  geom_point(data=afd_table, aes(x=position, y=afd), size=0.1)+
  geom_line(data=metric_table, aes_string(x="middle", y=tolower(metric)), color='red', size=0.3)+
  coord_cartesian(xlim=c(plot_start, plot_end),ylim=c(-1, 1))+
  scale_x_continuous(name=paste(chr,'(bp)',sep=' '), breaks=seq(plot_start,plot_end,by=10000))+
  scale_y_continuous(name=paste('AFD (', order,')', sep=''), sec.axis = sec_axis(trans=~.*1, name=toupper(metric)))+
  theme_classic()+
  geom_segment(data=genes, mapping=aes(x=genes$start, y=1.05, xend=genes$end, yend=1.05), arrow=arrow(length=unit(0.05, "cm"), type="closed"), size=0.25, color="grey40" )+
  #geom_segment(data=exons, mapping=aes(x=exons$start, y=1.05, xend=exons$end, yend=1.05), size=0.4, color="black")+
  theme(text = element_text(size=5), axis.line=element_line(colour='black',size=0.25), axis.ticks=element_line(colour='black',size=0.25))

print(p)
dev.off()


