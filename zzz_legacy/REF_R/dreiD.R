library(stringr)
library(dplyr)
library(RODBC)
library(ggplot2)
library(tidyr)
library(SphericalCubature)
library(RColorBrewer)

db <- list (essfinder="driver={SQL Server};server=consolework;database=EssFinder;trusted_connection=true",
            cs19="driver={SQL Server};server=consolework;database=CaseStudy19;trusted_connection=true",
            polar="driver={SQL Server};server=consolework;database=Polar;trusted_connection=true",
            polarneu="driver={SQL Server};server=consolework;database=PolarNeu;trusted_connection=true",
            essfinderII="driver={SQL Server};server=consolework;database=EssFinderII;trusted_connection=true")

ds <- data.frame(n=4L,dbs="essfinder",runid=48L, ispolar=FALSE, dof=2L, points=NA_integer_, pattern=NA_character_, slices=NA_integer_, zTransform=FALSE,
                 stringsAsFactors = FALSE)#1
ds <- rbind(ds,list(5L,"essfinder",47L,FALSE,2L,NA,NA,NA,FALSE))#2
ds <- rbind(ds,list(6L,"essfinder",45L,FALSE,3L,NA,NA,NA,FALSE))#3
ds <- rbind(ds,list(7L,"essfinder",46L,FALSE,3L,NA,NA,NA,FALSE))#4
ds <- rbind(ds,list(15L,"essfinder",94L,FALSE,3L,NA,NA,NA,FALSE))#5
ds <- rbind(ds,list(19L,"cs19",1L,FALSE,3L,NA,NA,NA,FALSE))#6
ds <- rbind(ds,list(4L,"polar",1L,TRUE,2L,NA,NA,NA,FALSE))#7
ds <- rbind(ds,list(5L,"polar",2L,TRUE,2L,NA,NA,NA,FALSE))#8
ds <- rbind(ds,list(6L,"polar",3L,TRUE,3L,NA,NA,NA,FALSE))#9
ds <- rbind(ds,list(7L,"polar",4L,TRUE,3L,NA,NA,NA,FALSE))#10
ds <- rbind(ds,list(13L,"polar",7L,TRUE,3L,NA,NA,NA,FALSE))#11
ds <- rbind(ds,list(15L,"polar",5L,TRUE,3L,NA,NA,NA,FALSE))#12
ds <- rbind(ds,list(19L,"polar",6L,TRUE,3L,NA,NA,NA,FALSE))#13
ds <- rbind(ds,list(8L,"polar",8L,TRUE,4L,NA,NA,NA,FALSE))#14
ds <- rbind(ds,list(9L,"polar",9L,TRUE,4L,NA,NA,NA,FALSE))#15
ds <- rbind(ds,list(17L,"polar",14L,TRUE,2L,20001,"{0, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 1, 1}",NA,FALSE))#16
ds <- rbind(ds,list(17L,"polar",15L,TRUE,4L,34476,"{0, 1, 2, 3, 1, 3, 4, 4, 2, 2, 4, 4, 3, 1, 3, 2, 1}",25,FALSE))#17
ds <- rbind(ds,list(10L,"polar",17L,TRUE,3L,501501,"{0, 1, 2, 1, 2, 3, 2, 1, 2, 1}",NA,FALSE))#18
ds <- rbind(ds,list(14L,"polar",16L,TRUE,3L,501501,"{0, 1, 2, 1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 1}",NA,FALSE))#19
ds <- rbind(ds,list(14L,"polar",21L,TRUE,4L,262701,"{0, 1, 2, 3, 3, 2, 1, 4, 1, 2, 3, 3, 2, 1}",50,FALSE))#20
ds <- rbind(ds,list(16L,"polar",18L,TRUE,2L,100001,"{0, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1}",NA,FALSE))#21
ds <- rbind(ds,list(16L,"polar",19L,TRUE,4L,262701,"{0, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1}",50,FALSE))#22
ds <- rbind(ds,list(7L,"polar",27L,TRUE,3L,20301,NA,NA,TRUE))#23
ds <- rbind(ds,list(7L,"polar",31L,TRUE,3L,20301,NA,NA,TRUE))#24
ds <- rbind(ds,list(13L,"polar",32L,TRUE,3L,20301,"{0, 1, 2, 2, 3, 1, 3, 1, 3, 2, 2, 1}",NA,TRUE))#25
ds <- rbind(ds,list(7L,"polarneu",1L,TRUE,3L,20301,NA,NA,TRUE))#26
ds <- rbind(ds,list(13L,"polarneu",2L,TRUE,3L,20301,"{0, 1, 2, 2, 3, 1, 3, 1, 3, 2, 2, 1}",NA,TRUE))#27
ds <- rbind(ds,list(13L,"polarneu",12L,TRUE,3L,20301,"{0, 1, 2, 2, 3, 1, 3, 1, 3, 2, 2, 1}",NA,TRUE))#28
ds <- rbind(ds,list(19L,"polarneu",15L,TRUE,3L,NA,NA,NA,TRUE))#29
ds <- rbind(ds,list(19L,"polarneu",22L,TRUE,3L,NA,NA,NA,TRUE))#30

getfilename <- function(x) {
  x <- x %>% str_replace("\\{","") %>% str_replace("\\}","") %>% str_replace_all(",","") %>% str_replace_all(" ","")
  return(paste(x,collapse = "_"))
}

#gets data from database and writes table in a file given by "getfilename", takes one row of ds as input
generatefile <- function(ds) {
  str(ds)
  dbhandle <- odbcDriverConnect(db[[ds$dbs]])
  
  stringsql <- str_c("select m.matrixid, m.matrixname, m.anzahless, v.support  
                     from matrizen m inner join vektoren v on v.matrixid=m.matrixid and m.fk_run=v.fk_run 
                     where m.fk_run=",ds$runid,"  and isess=1 and m.anzahless!=0 order by m.matrixid, v.support")
  stringsql2 <- str_c("select m.matrixid, m.matrixname, m.anzahless from matrizen m 
                      where m.fk_run=",ds$runid,"  and m.anzahless=0 order by m.matrixid")
  
  tab1 <- sqlQuery(dbhandle,stringsql,stringsAsFactors=FALSE)
  tab2 <- sqlQuery(dbhandle,stringsql2,stringsAsFactors=FALSE)
  table <- tab1  %>% 
    group_by(matrixid, matrixname, anzahless) %>% 
    summarise(vector=paste(str_replace_all(support,"\\,",""), collapse=".")) %>%
    bind_rows(tab2)
  
  write.table(table,getfilename(ds))  
}

# takes 1 row of ds as input, loads the date via "getfilename" and does all calculations
dodata <- function(ds) {
  #ds <- ds[6,]
  str(ds)
  table <- read.table(getfilename(ds),header=TRUE,stringsAsFactors=FALSE)
  U <- matrix(c(-sqrt(3),1,sqrt(2),sqrt(3),1,sqrt(2),0,-2,sqrt(2)),nrow=3,ncol=3)
  
  essdigits <- nchar(max(table$anzahless),type="chars")

  if (ds$ispolar){
    table <- table %>% separate(matrixname,c("phi","theta","theta2"),sep="\\;",remove=FALSE) %>%
      mutate(phi=as.numeric(str_replace(phi,",",".")),#phitemp
             #phi=ifelse(phitemp>1,phitemp-2,phitemp),
             theta=as.numeric(str_replace(theta,",",".")), 
             theta2=as.numeric(str_replace(theta2,",",".")), 
             Number_of_ESS=as.factor(anzahless),
             Supports=as.factor(str_c(formatC(anzahless, width = essdigits, format = "d", flag = "0"),"-",vector))) %>%
      select(-anzahless,-vector)
    
    
#     if (ds$dof==3L) {
#       
#       if (!ds$zTransform) {
#         dmatrix <- t(data.matrix(table %>% select(theta,phi)))
#         dmatrix2 <- polar2rect(rep(1,nrow(table)),dmatrix*pi)
#         dmatrix3 <- dmatrix2[c(2,3,1),]
#         dmatrix4 <- (U %*% dmatrix3)[c(3,1,2),]
#         dmatrix5 <- t((rect2polar(dmatrix4)$phi)/pi)
#         colnames(dmatrix5) <- c("z_theta","z_phi")
#         table <- cbind(table, dmatrix5)
#       }
#     }
    
  } else {
    table <- table %>% separate(matrixname,c("x","y","z"),sep="\\,",remove=FALSE) %>%
      mutate(Number_of_ESS=as.factor(anzahless),Supports=as.factor(str_c(formatC(anzahless, width = essdigits, format = "d", flag = "0"),"-",vector))) %>%
      select(-anzahless,-vector)
    
#     if (ds$dof==2L) {
#       dmatrix <- t(data.matrix(table %>% select(x,y)))
#       dmatrix2 <- t((rect2polar(dmatrix)$phi)/pi)
#       colnames(dmatrix2) <- c("phi")
#       table <- cbind(table, dmatrix2)
#     }
#     
#     if (ds$dof==3L) {
#       dmatrix <- t(data.matrix(table %>% select(x,y,z)))
#       dmatrix2 <- t((rect2polar(dmatrix)$phi)/pi)
#       colnames(dmatrix2) <- c("theta","phi")
#       table <- cbind(table, dmatrix2)
#       dmatrixx <- (U %*% dmatrix3) #[c(3,1,2),]
#       dmatrix5 <- t((rect2polar(dmatrix4)$phi)/pi)
#       colnames(dmatrix5) <- c("z_theta","z_phi")
#       table <- cbind(table, dmatrix5)
#      }
  }
  write.table(table,str_c(getfilename(ds),"_2_")) 
}

# plotting the function
doplot <- function(ds) {
  table <- read.table(str_c(getfilename(ds),"_2_"),header=TRUE,stringsAsFactors=FALSE) %>%
    mutate(Number_of_ESS=as.factor(Number_of_ESS),Supports=as.factor(Supports))

  getPalette = colorRampPalette(brewer.pal(9, "RdGy"))
  
  if (ds$ispolar) {
    if (ds$dof==2L) {
      
      ggplot(table, aes(phi, Number_of_ESS)) + geom_point(aes(color=Number_of_ESS),size=1) + guides(fill=guide_legend(ncol=100)) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_ESS"), hjust=-.2, vjust=2) +
        scale_colour_manual(values=getPalette(length(unique(table$Number_of_ESS)))) +
        guides(colour = guide_legend(override.aes = list(size=5,linetype=0)))
      ggsave(str_c(getfilename(ds),"_ESS.pdf"), width=30, height=20, units="cm")
      
      ggplot(table, aes(phi,Supports )) + geom_point(aes(color=Supports),size=1) + guides(fill=guide_legend(ncol=100)) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_Supports"), hjust=-.2, vjust=2) +
        scale_colour_manual(values=getPalette(length(unique(table$Supports)))) +
      guides(colour = guide_legend(override.aes = list(size=5,linetype=0)))
      ggsave(str_c(getfilename(ds),"_Supports.pdf"), width=30, height=20, units="cm")
      
    }
    if (ds$dof==3L) {
      
      ggplot(table, aes(phi, theta)) + geom_raster(aes(fill=Number_of_ESS)) + #guides(fill=guide_legend(ncol=100)) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_ESS"), hjust=-.2, vjust=2) +
        scale_fill_manual(values=getPalette(length(unique(table$Number_of_ESS))))
      ggsave(str_c(getfilename(ds),"_ESS.pdf"), width=30, height=20, units="cm")
      
      ggplot(table, aes(phi, theta)) + geom_raster(aes(fill=Supports)) + #guides(fill=guide_legend(ncol=100)) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_Supports"), hjust=-.2, vjust=2) +
        scale_fill_manual(values=getPalette(length(unique(table$Supports))))
      ggsave(str_c(getfilename(ds),"_Supports.pdf"), width=30, height=20, units="cm")
      
    }  
    if (ds$dof==4L) {
      
      generateSlices <- function(value) {
        tab2 <- table %>% filter(abs(theta2-value)<1e-05)
        #print(nrow(tab2))
        plot <- ggplot(tab2, aes(phi, theta)) + geom_raster(aes(fill=Number_of_ESS)) + #guides(fill=guide_legend(ncol=100)) + 
          annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_ESS, theta2=",value,"*pi"), hjust=-.2, vjust=2) +
          scale_fill_manual(values=getPalette(length(unique(tab2$Number_of_ESS))))
        sup <- ggplot(tab2, aes(phi, theta)) + geom_raster(aes(fill=Supports)) + #guides(fill=guide_legend(ncol=100)) + 
          annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_Supports"), hjust=-.2, vjust=2) +
          scale_fill_manual(values=getPalette(length(unique(tab2$Supports))))
        return(list(plot,sup))
      }
      
      result <- lapply(seq(0,1,by=1/25),generateSlices)
      plot_ess <- lapply(result,function(x) return(x[[1]]))
      plot_sup <- lapply(result,function(x) return(x[[2]]))
      
      #pdf(str_c(getfilename(ds),"_ESS.pdf"),width=11.7, height=8.3); marrangeGrob(plot_ess, nrow = 1, ncol = 1) ; dev.off()
      pdf(str_c(getfilename(ds),"_ESS.pdf"),width=11.7, height=8.3)
      lapply(plot_ess,plot)
      dev.off()
      #pdf(str_c(getfilename(ds),"_Supports.pdf"),width=11.7, height=8.3)
      #lapply(plot_sup,plot)
      #dev.off()
      #pdf(str_c(getfilename(ds),"_Supports.pdf"),width=11.7, height=8.3); marrangeGrob(plot_sup, nrow = 1, ncol = 1) ; dev.off()
    }
    
#   } else {
#     if (ds$dof==2L) {
#       
#       ggplot(table, aes(phi, Number_of_ESS)) + geom_point(aes(color=Number_of_ESS)) + #guides(fill=guide_legend(ncol=100)) + 
#         annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_ESS"), hjust=-.2, vjust=2) +
#         scale_colour_manual(values=getPalette(length(unique(table$Number_of_ESS)))) 
#       ggsave(str_c(getfilename(ds),"_ESS.pdf"), width=30, height=20, units="cm")
#       
#       ggplot(table, aes(phi, Supports)) + geom_point(aes(color=Supports)) + #guides(fill=guide_legend(ncol=100)) + 
#         annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_Supports"), hjust=-.2, vjust=2)+
#         scale_colour_manual(values=getPalette(length(unique(table$Supports)))) 
#       ggsave(str_c(getfilename(ds),"_Supports.pdf"), width=30, height=20, units="cm")
#       
#     }
#     if (ds$dof==3L) {
#       
#       ggplot(table, aes(phi, theta)) + geom_point(aes(color=Number_of_ESS),size=0.5) + 
#         annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_ESS"), hjust=-.2, vjust=2) +
#         scale_colour_manual(values=getPalette(length(unique(table$Number_of_ESS)))) +
#         guides(colour = guide_legend(override.aes = list(size=5,linetype=0)))  
#       ggsave(str_c(getfilename(ds),"_ESS.pdf"), width=30, height=20, units="cm")
#     }     
  }
}

#plot after z-transformation
doplot_z <- function(ds) {
  table <- read.table(str_c(getfilename(ds),"_FERTIG"),header=TRUE,stringsAsFactors=FALSE) %>%
    mutate(Number_of_ESS=as.factor(Number_of_ESS),Supports=as.factor(Supports))
  getPalette = colorRampPalette(brewer.pal(9, "RdGy"))
  
  if (ds$ispolar) {
    if (ds$dof==3L) {
      ggplot(table, aes(z_phi, z_theta)) + geom_point(aes(color=Number_of_ESS),size=0.5) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_z_ungenau_ESS"), hjust=-.2, vjust=2) +
        scale_colour_manual(values=getPalette(length(unique(table$Number_of_ESS)))) +
        guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) 
      ggsave(str_c(getfilename(ds),"_z_ungenau_ESS.pdf"), width=30, height=20, units="cm")
      
      ggplot(table, aes(z_phi, z_theta)) + geom_point(aes(color=Supports),size=0.5) + 
        annotate("text", x=-Inf, y=Inf, label=str_c(getfilename(ds),"_z_ungenau_Supports"), hjust=-.2, vjust=2) +
        scale_colour_manual(values=getPalette(length(unique(table$Supports)))) +
        guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) 
      ggsave(str_c(getfilename(ds),"_z_ungenau_Supports.pdf"), width=30, height=20, units="cm")
    }
  }
}

# all together
doall <- function(ds) {
  generatefile(ds)
  dodata(ds)
  doplot(ds)
}


ds.input <- ds[29,]#[9:12,]

for (i in seq_len(nrow(ds.input))) {
  row <- ds.input[i,]
  print(getfilename(row))
  generatefile(row)
  dodata(row)
  doplot(row)
}


           