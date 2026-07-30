library(stringr)
library(dplyr)
library(RODBC)
library(rvest)
#library(profr)
library(ggplot2)
library(xtable)
#library(devtools)
library(tidyr)
library(gridExtra)


######################################################CASE STUDY 19
#dbhandle <- odbcDriverConnect('driver={SQL Server};server=consolework;database=EssFinder;trusted_connection=true')

#-------------------------------------------------------------------------------------------------------------------------------------------
stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
where m.fk_run=45  and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
where m.fk_run=45  and m.anzahless=0 order by m.matrixid"


# testtab <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  
# testtab2 <- testtab %>%
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-"))
# 
# t1 <- sqlQuery(dbhandle,"select * from matrizen where fk_run=45 and anzahless!=0",stringsAsFactors=FALSE)
# x <- testtab2 %>% ungroup() %>%select(matrixid) 
# y <- t1 %>% select(matrixid)
# t2 <- setdiff(x,y)
# 
# testtab3 <- sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)  


# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","y","z"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y),z=as.integer(z), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))


#x1 <- table %>% filter(anzahless==9,z==-100)#,75<x,x<90, -89<y,y<(-80),z==-100)

plot_ess <- list()
generate <- function(xx) {
  tab2 <- table %>% filter(z==xx)
 
  plot <- ggplot(tab2, aes(x, y, fill=Number_of_ESS)) +
    geom_raster() + annotate("text", x=-Inf, y=Inf, label=str_c("n=6",", C(x,y,",xx,",y,x,0)"), hjust=-.2, vjust=2)
  
  plot_ess<<-c(plot_ess,list(plot))
  
  #ggsave(str_c("n_6_wert_",xx,"_ESS.pdf"), width=30, height=20, units="cm")
  
#   ggplot(tab2, aes(x, y, fill=Supports)) +
#     geom_raster() + annotate("text", x=-Inf, y=Inf, label=str_c("n=6",", C(x,y,",xx,",y,x,0)"), hjust=-.2, vjust=2)
#   ggsave(str_c("n_6_wert_",xx,"_Supports.pdf"), width=30, height=20, units="cm")   
  
}

lapply(seq(-100,100,by=25),generate)

#multiplot(plot_ess,cols=2)
#do.call(grid.arrange, c(plot_ess,list(nrow=2,ncol=2)))
marrangeGrob(grobs=plot_ess, nrow = 1, ncol = 2)

# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,10,y,10,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_10_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,10,y,10,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_10_2.pdf", width=30, height=20, units="cm")
# 
# #-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
# from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
# where m.fk_run=45 and m.matrixname like '%,-10,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
# where m.fk_run=45 and m.matrixname like '%,-10,%' and m.anzahless=0 order by m.matrixid"
# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-10,y,-10,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg10_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-10,y,-10,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg10_2.pdf", width=30, height=20, units="cm")
# 
# #-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
#   from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
#   where m.fk_run=45 and m.matrixname like '%,50,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
#   where m.fk_run=45 and m.matrixname like '%,50,%' and m.anzahless=0 order by m.matrixid"
# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#    group_by(matrixid, matrixname, anzahless) %>% 
#    summarise(vector=paste(support, collapse="-")) %>%
#    bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#    separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#    mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,50,y,50,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_50_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,50,y,50,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_50_2.pdf", width=30, height=20, units="cm")
# 
# #-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
#   from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
#   where m.fk_run=45 and m.matrixname like '%,-50,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
#   where m.fk_run=45 and m.matrixname like '%,-50,%' and m.anzahless=0 order by m.matrixid"
# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-50,y,-50,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg50_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-50,y,-50,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg50_2.pdf", width=30, height=20, units="cm")
# 
# #-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
# from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
# where m.fk_run=45 and m.matrixname like '%,100,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
# where m.fk_run=45 and m.matrixname like '%,100,%' and m.anzahless=0 order by m.matrixid"
# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,100,y,100,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_100_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,100,y,100,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_100_2.pdf", width=30, height=20, units="cm")
# 
# #-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
# from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
# where m.fk_run=45 and m.matrixname like '%,-100,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
# where m.fk_run=45 and m.matrixname like '%,-100,%' and m.anzahless=0 order by m.matrixid"
# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
# 
# ggplot(table, aes(x, y, fill=Number_of_ESS)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-100,y,-100,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg100_2_Anzahl.pdf", width=30, height=20, units="cm")
# 
# ggplot(table, aes(x, y, fill=Supports)) +
#   geom_raster() + annotate("text", x=-Inf, y=Inf, label="C(x,-100,y,-100,x,0)", hjust=-.2, vjust=2)
# ggsave("table_n6_neg100_2.pdf", width=30, height=20, units="cm")
