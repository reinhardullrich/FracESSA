library(stringr)
library(dplyr)
library(RODBC)
library(rvest)
#library(profr)
library(ggplot2)
library(xtable)
#library(devtools)
library(tidyr)


######################################################CASE STUDY 19
dbhandle <- odbcDriverConnect('driver={SQL Server};server=consolework;database=EssFinder;trusted_connection=true')

#-------------------------------------------------------------------------------------------------------------------------------------------
# stringsql <- "select m.matrixid, m.matrixname, m.anzahless, v.support  
# from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
# where m.fk_run=45 and m.matrixname like '%,10,%' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support"
# stringsql2 <- "select m.matrixid, m.matrixname, m.anzahless from matrizen m 
# where m.fk_run=45 and m.matrixname like '%,10,%' and m.anzahless=0 order by m.matrixid"

toTest <- 
        list(4,0,0,48,"%","x,y,x,0") %>%
#  rbind(list(5,0,0,47,"%","x,y,y,x,0")) 

# toTest <- 
#         list(6,1, 10,45,"10,%,%","10,x,y,x,10,0") %>%
#   rbind(list(6,1,-10,45,"-10,%,%","-10,x,y,x,-10,0")) %>%
#   rbind(list(6,1, 50,45,"50,%,%","50,x,y,x,50,0")) %>%
#   rbind(list(6,1,-50,45,"-50,%,%","-50,x,y,x,-50,0")) %>%
#     #
#   rbind(list(6,2, 10,45,"%,10,%","x,10,y,10,x,0")) %>%
#   rbind(list(6,2,-10,45,"%,-10,%","x,-10,y,-10,x,0")) %>%
#   rbind(list(6,2, 50,45,"%,50,%","x,50,y,-50,x,0")) %>%
#   rbind(list(6,2,-50,45,"%,-50,%","x,-50,y,-50,x,0")) %>%
#     #
#   rbind(list(6,3, 10,45,"%,%,10","x,y,10,y,x,0")) %>%
#   rbind(list(6,3,-10,45,"%,%,-10","x,y,-10,y,x,0")) %>%
#   rbind(list(6,3, 50,45,"%,%,50","x,y,50,y,x,0")) %>%
   rbind(list(6,3,-100,45,"%,%,-50","x,y,-50,y,x,0")) #%>%
#     #
#   rbind(list(7,1, 10,46,"10,%,%","10,x,y,y,x,10,0")) %>%
#   rbind(list(7,1,-10,46,"-10,%,%","-10,x,y,y,x,-10,0")) %>%
#   rbind(list(7,1, 50,46,"50,%,%","50,x,y,y,x,50,0")) %>%
#   rbind(list(7,1,-50,46,"-50,%,%","-50,x,y,y,x,-50,0")) %>%
#   #
#   rbind(list(7,2, 10,46,"%,10,%","x,10,y,y,10,x,0")) %>%
#   rbind(list(7,2,-10,46,"%,-10,%","x,-10,y,y,-10,x,0")) %>%
#   rbind(list(7,2, 50,46,"%,50,%","x,50,y,y,-50,x,0")) %>%
#   rbind(list(7,2,-50,46,"%,-50,%","x,-50,y,y,-50,x,0")) %>%
#   #
#   rbind(list(7,3, 10,46,"%,%,10","x,y,10,10,y,x,0")) %>%
#   rbind(list(7,3,-10,46,"%,%,-10","x,y,-10,-10,y,x,0")) %>%
#   rbind(list(7,3, 50,46,"%,%,50","x,y,50,50,y,x,0")) %>%
#   rbind(list(7,3,-50,46,"%,%,-50","x,y,-50,-50,y,x,0"))
  
toTest <- as.data.frame(toTest)
names(toTest) <- c("n","position","wert","fk_run","sql","name")


generate <- function(xx) {
  
  sql <- str_c("select m.matrixid, m.matrixname, m.anzahless, v.support  
    from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
    where m.fk_run=",xx['fk_run']," and m.matrixname like '",xx['sql'],"' and isess=1 and m.anzahless!=0 order by m.matrixid, v.support")
  
  sql2 <- str_c("select m.matrixid, m.matrixname, m.anzahless from matrizen m 
    where m.fk_run=",xx['fk_run']," and m.matrixname like '",xx['sql'],"' and m.anzahless=0 order by m.matrixid")
  print(sql)
  
  if (xx['position']==1) {
    switcher <- c("nix","x","y")
  } else if (xx['position']==2) {
    switcher <- c("x","nix","y")
  } else if (xx['position']==3) {
    switcher <- c("x","y","nix")
  } else if (xx['position']==0) {
    switcher <- c("x","y")
  }
  
  table <- sqlQuery(dbhandle,sql,stringsAsFactors=FALSE)  %>% 
    group_by(matrixid, matrixname, anzahless) %>% 
    summarise(vector=paste(support, collapse="-")) %>%
    bind_rows(sqlQuery(dbhandle,sql2,stringsAsFactors=FALSE)) %>% 
    separate(matrixname,switcher,sep="\\,",remove=FALSE) %>% 
    mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
  
  #tt<<- table
  #tt1 <- tt %>% filter(anzahless==9)
  
  

  ggplot(table, aes(x, y, fill=Number_of_ESS)) +
    geom_raster() + annotate("text", x=-Inf, y=Inf, label=str_c("n=",xx['n'],", C(",xx['name'],")"), hjust=-.2, vjust=2)
  ggsave(str_c("table_",xx['n'],"_",xx['position'],"_",xx['wert'],"_ESS.pdf"), width=30, height=20, units="cm")
  
  ggplot(table, aes(x, y, fill=Supports)) +
    geom_raster() + annotate("text", x=-Inf, y=Inf, label=str_c("n=",xx['n'],", C(",xx['name'],")"), hjust=-.2, vjust=2)
    ggsave(str_c("table_",xx['n'],"_",xx['position'],"_",xx['wert'],"_Supports.pdf"), width=30, height=20, units="cm")
  #return(sql)
  
}

#generate("lalalal")
#fk_run <- c(1,2,3)
#sql <- c("sdlkf","sdklfj","cxxxsdklfj")
#y <- data.frame(fk_run,sql)
apply(toTest,1,generate)

# 
# table <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)  %>% 
#   group_by(matrixid, matrixname, anzahless) %>% 
#   summarise(vector=paste(support, collapse="-")) %>%
#   bind_rows(sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)) %>%
#   separate(matrixname,c("x","nix","y"),sep="\\,",remove=FALSE) %>% 
#   mutate(x=as.integer(x),y=as.integer(y), Number_of_ESS=as.factor(anzahless),Supports=as.factor(paste(anzahless,vector,sep="...")))
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
