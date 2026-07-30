library(stringr)
library(dplyr)
library(RODBC)
library(rvest)
library(ggplot2)
library(xtable)
library(grid)
library(gridExtra)

db <- list (essfinder="driver={SQL Server};server=consolework;database=EssFinder;trusted_connection=true",
            cs19="driver={SQL Server};server=consolework;database=CaseStudy19;trusted_connection=true",
            polar="driver={SQL Server};server=consolework;database=Polar;trusted_connection=true",
            polarneu="driver={SQL Server};server=consolework;database=PolarNeu;trusted_connection=true",
            essfinderII="driver={SQL Server};server=consolework;database=EssFinderII;trusted_connection=true")

ds <- data.frame(n=6L,runid=1L, tabname="table06",dbhandle="essfinderII",stringsAsFactors = FALSE)
ds <- rbind(ds,list(7L,2L,"table07","essfinderII"))
ds <- rbind(ds,list(19L,1L,"table19","cs19"))

createData <- function (ds) {
  dbhandle <- odbcDriverConnect(db[ds$dbhandle])
  sql <- str_c("select m.matrixid, m.dimension,m.anzahless,v.supportnummer, v.extendedsupportnummer, isess, grundess from matrizen m inner join vektoren v on m.matrixid=v.matrixid and m.fk_run=v.fk_run where m.fk_run=",ds$runid)
  worktab <<- sqlQuery(dbhandle,sql,stringsAsFactors=FALSE)
  save(worktab,file=str_c(ds$tabname,".rda"))
}

getData <- function(ds) {
  load(str_c(ds$tabname,".rda"))
  worktab <<- worktab
}


plot01 <- function(ds) {
  tab <- worktab %>% select(matrixid, anzahless) %>% distinct() 
  x = ggplot(tab, aes(x=factor(anzahless))) + 
    geom_bar() +
    geom_text(aes(label=..count..), hjust=-0.1,angle=90,stat="bin") +
    #ylim(0,450000) +
    labs(x="Number of ESS", y="Count")
  #ggsave(str_c(ds$tabname,"_1.pdf"), width=30, height=10, units="cm")
  return(x)
}

plot02 <- function(ds) {
  tab <- worktab %>% filter(isess==TRUE) %>% group_by(anzahless,supportnummer) %>% summarise(anzahl=n())
  tab <- tab %>% group_by(anzahless) %>% mutate(percentage=anzahl/sum(anzahl))
  ggplot(tab,aes(x=factor(anzahless),y=percentage, fill=factor(supportnummer))) +
    geom_bar(stat="identity") + 
    guides(fill=guide_legend(reverse=TRUE)) +
    labs(x="Number of ESS", y="Percentage", fill="Supportsize") +
    #scale_fill_brewer(palette="Paired")
  ggsave(str_c(ds$tabname,"_2.pdf"), width=30, height=10, units="cm")
}

plot03 <- function(ds) {
  legend <- c("True: Pure ESS","True: PosDef", "True: CoPos", "True: PosDef","True: CoPos","False","False","False")
  tab <- worktab %>% group_by(supportnummer,grundess) %>% summarise(anzahl=n())
  tab <- tab %>% mutate(candidateisess=legend[grundess]) %>% group_by(supportnummer) %>% mutate(percentage=anzahl/sum(anzahl))
  ggplot(tab,aes(x=factor(supportnummer),y=percentage, fill=candidateisess)) +
    geom_bar(stat="identity") + 
    guides(fill=guide_legend(reverse=TRUE)) +
    labs(x="Supportsize", y="Percentage", fill="Candidate is ESS") + 
    #scale_fill_brewer(palette="Paired")
  ggsave(str_c(ds$tabname,"_3.pdf"), width=30, height=10, units="cm")
}

plotall <- function(ds) {
  plot01(ds)
  plot02(ds)
  plot03(ds)
}

doall <- function() {
  for (i in seq_len(nrow(ds))){
    myds <- ds[i,]
    getData(myds)
    plotall(myds)
  }
}

x <-  by(ds, 1:nrow(ds), function(x) { getData(x);plot01(x) })
str(x[[1]])
do.call(grid.arrange,x)
grid.arrange(x[[1]],x[[2]],x[[3]], ncol = 1, main = "Main title")


#myds <- ds[1,]
#getData(myds)
#plotall(myds)
#createData(myds)
#getData(myds)
#plot01(myds)
#plot02(myds)





# 
# stringsql <- "select AnzahlESS 'ESS', count(*)'Anzahl' from matrizen group by anzahless order by anzahless"
# 
# table19_1 <- sqlQuery(dbhandle19,stringsql,stringsAsFactors=FALSE)
# table19_1 <- table19_1 %>% mutate( RelAnzahl=round(Anzahl/sum(Anzahl),digits=4)) %>%mutate(ESS=as.factor(ESS))
# save(table19_1,file="table19_1.rda")
# load("table19_1.rda")
# 
# table08_1 <- sqlQuery(dbhandle08,stringsql,stringsAsFactors=FALSE)
# table08_1 <- table08_1 %>% mutate( RelAnzahl=round(Anzahl/sum(Anzahl),digits=4)) %>%mutate(ESS=as.factor(ESS))
# save(table08_1,file="table08_1.rda")
# load("table08_1.rda")
# 
# ggplot(table19_1, aes(x=ESS, y=Anzahl)) + geom_bar(stat="identity") + geom_text(aes(label=Anzahl), hjust=-0.1,angle=90) +xlab("Number of ESS") + ylab("Count")+ylim(0,600000)
# #ggplot(table19_1, aes(x=ESS, y=RelAnzahl)) + geom_bar(stat="identity") + geom_text(aes(label=RelAnzahl), hjust=-0.1,angle=90) #+ylim(0,37000)
# ggsave("table19_1.pdf", width=30, height=10, units="cm")
# ggplot(table08_1, aes(x=ESS, y=Anzahl)) + geom_bar(stat="identity") + geom_text(aes(label=Anzahl), hjust=-0.1,angle=90) +xlab("Number of ESS") + ylab("Count")+ylim(0,120000000)
# #ggplot(table08_1, aes(x=ESS, y=RelAnzahl)) + geom_bar(stat="identity") + geom_text(aes(label=RelAnzahl), hjust=-0.1,angle=90) #+ylim(0,37000)
# ggsave("table08_1.pdf", width=30, height=10, units="cm")
# 
# ###########################################################
# stringsql <- "select m.AnzahlESS 'ESS', v.supportnummer 'SupportSize', count(*) 'Anzahl' from matrizen m inner join vektoren v on m.matrixid=v.matrixid and m.fk_run=v.fk_run 
#                 where  v.isess=1  group by m.anzahless, v.supportnummer order by  m.anzahless, v.supportnummer"
# table19_2 <- sqlQuery(dbhandle19,stringsql,stringsAsFactors=FALSE)
# table19_2 <- table19_2 %>% mutate (ESS=as.factor(ESS), SupportSize=as.factor(SupportSize)) %>% 
#   group_by(ESS) %>% 
#   mutate(Percentage=Anzahl/sum(Anzahl))
# save(table19_2,file="table19_2.rda")
# load("table19_2.rda")
# ggplot(table19_2,aes(x=ESS,y=Percentage, fill=SupportSize)) +geom_bar(stat="identity") + guides(fill=guide_legend(reverse=TRUE)) +
#      scale_fill_brewer(palette="Paired") 
# ggsave("table19_2.pdf", width=30, height=10, units="cm")
# 
# ###############################################################
# legend <- c("True: Pure ESS","True: PosDef", "True: CoPos", "True: PosDef","True: CoPos","False","False","False")
# stringsql <- "select v.supportnummer 'SupportSize', v.grundess 'ESS', count(*) 'Anzahl' from matrizen m inner join vektoren v on 
#     m.matrixid=v.matrixid and m.fk_run=v.fk_run group by  v.supportnummer, v.grundess order by  v.supportnummer, v.grundess"
# table19_3 <- sqlQuery(dbhandle19,stringsql,stringsAsFactors=FALSE)
# table19_3 <- table19_3 %>% mutate(SupportSize=as.factor(SupportSize), CandidateIsESS=legend[ESS]) %>%
#     group_by(SupportSize) %>%
#     mutate(Percentage=Anzahl/sum(Anzahl))
# save(table19_3,file="table19_3.rda")
# load("table19_3.rda")
# ggplot(table19_3,aes(x=SupportSize,y=Percentage, fill=CandidateIsESS)) +geom_bar(stat="identity") + guides(fill=guide_legend(reverse=TRUE)) +
#   labs(fill="Candidate is ESS")+scale_fill_brewer(palette="Paired") 
# ggsave("table19_3.pdf", width=30, height=10, units="cm")
