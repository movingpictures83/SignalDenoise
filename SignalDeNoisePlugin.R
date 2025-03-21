read_kraken_reports <- function(files, sample_names = NULL, study_name = NULL, min_reads = 2, min_uniq = 2){
  require(data.table)
  require(tibble)
  
  if(is.null(sample_names)){sample_names = files}
  if(is.null(study_name)){study_name = NA}
  if(length(study_name) == 1){study_name = rep(study_name, length(files))}
  
  df = list()
  n = 0
  for(i in 1:length(files)){
    if(round(i/length(files)*100, 2) > n){n = round(i/length(files)*100, 2); cat(paste0('\r',n,'% done   '))}
    x = read.delim(files[i], header = F)
    x$V8 = trimws(x$V8)
    total_reads = x$V2[1] + x$V2[2]
    n_microbiome_reads = sum(x$V2[x$V8 %in% c('Bacteria', 'Fungi', 'Viruses')])
    df[[i]] = data.frame(study = study_name[i], sample = sample_names[i],
                         rank = x$V6, taxid = x$V7, name = x$V8, 
                         reads = x$V2, min = x$V4, uniq = x$V5, 
                         rpm = x$V2/total_reads*10^6,
                         rpmm = x$V2/n_microbiome_reads*10^6)
  }
  df = rbindlist(df) %>% tibble()
  cat('\n')
  
  df = subset(df, reads >= min_reads & uniq >= min_uniq) 
  #df
}

#kr = read_kraken_reports(c("full1.kraken.report.txt", "full2.kraken.report.txt", "full3.kraken.report.txt"))


#dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
#source("RPluMA.R")
source("RIO.R")
require(ggplot2)
require(dplyr)

input <- function(inputfile) {
   reports = readSequential(inputfile)
   kr <<- read_kraken_reports(reports)
#kr = read_kraken_reports(c("HPylori/SRR9713132/output.kraken.report.txt","HPylori/SRR9713133/output.kraken.report.txt","HPylori/SRR9713134/output.kraken.report.txt","HPylori/SRR9713135/output.kraken.report.txt","HPylori/SRR9713136/output.kraken.report.txt","HPylori/SRR9713137/output.kraken.report.txt","HPylori/SRR9713138/output.kraken.report.txt","HPylori/SRR9713139/output.kraken.report.txt"))
}

run <- function() {
#kr = readRDS("exampleData/zhang.reports.RDS")
#kr
# remove taxa detected in < 3 samples
kr <<- kr %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 2) %>%
  select(-nn)

#kr = kr %>%
#	group_by(taxid) #%>%
#        mutate(nn = n())
#	subset(nn > 2)
	
#print(kr[,"nn"])
# run correlations
c2 <<- kr %>%
  subset(rank %in% c('G', 'S')) %>%
  group_by(name) %>%
  summarize(r1 = cor(min,uniq,method='spearman'),
            r2 = cor(min,reads,method='spearman'),
            r3 = cor(reads,uniq,method='spearman'),
            p1 = cor.test(min,uniq,method='spearman')$p.value,
            p2 = cor.test(min,reads,method='spearman')$p.value,
            p3 = cor.test(reads,uniq,method='spearman')$p.value
            )
}

output <- function(outputfile) {
ggplot(c2, aes(r1, r2, color = r3))+
  geom_point() +
  geom_density_2d(size = 0.75, colour = "black") 
saveRDS(kr, paste(outputfile, "kr.rds", sep="/"))
saveRDS(c2, paste(outputfile, "c2.rds", sep="/"))
}
