require(dplyr)
#enn <-       as.integer(c(1, 2, 3, 4, 5, 6,  7,  8,  9, 10, 11, 12,  13,  14,  15,  16,  17,  18,   19,   20))
#foundEss <-  as.integer(c(1, 2, 3, 4, 5, 6, 14, 20, 27, 40, 33, 96, 143, 196, 360, 256, 306, 630, 1444, 1024))
#suppsize <- as.integer(c (1, 1, 1, 2, 2, 4,  3,  4,  5,  5,  7,  6,   6,   7,   7,   8,   8,   9,    9,   10))

#all <-  list()
# dimension, anzahl_ess, supportsize

tab1 <- data.frame(t(data.frame(
  #c(1,1,1),
  #c(2,1,2),
  c(2,2,1),
  #c(3,1,3),
  c(3,3,1),
  #c(4,1,4),
  c(4,4,2),
  #c(5,1,5),
  c(5,5,3),
  #c(6,1,6),
  c(6,6,4),
  c(6,8,3),
  c(6,9,2),
  #c(7,1,7),
  c(7,7,5),
  c(7,14,3),
  #c(8,1,8),
  c(8,8,6),
  c(8,10,4),# ist auch drinnen suppsize=5 8x und suppsize=4 2x !!!!
  c(8,20,4),
  #c(9,1,9),
  c(9,9,7),
  c(9,18,5),
  c(9,30,3),
  #c(10,1,10),
  c(10,10,8), # auch ess=12 mit suppsize=7 10x und suppsize=5 x2 -------> aber 10,10,8 liefert besseres ergebnis!!
  c(10,25,6),
  c(10,32,5),
  c(10,40,3),
  #c(11,1,11),
  c(11,11,9),
  c(11,33,7),
  c(11,66,5),
  #c(12,1,12),
  c(12,14,9),
  c(12,36,8),
  c(12,96,6),
  #c(13,1,13),
  c(13,13,11),
  c(13,39,9),
  c(13,130,7),
  c(13,143,6),
  c(14,14,11),
  c(14,35,9),# gibts nicht, ist eine obere schranke
  c(14,49,8),
  c(14,196,7),
  c(15,15,12),
  c(15,20,11),
  c(15,30,10),
  c(15,125,9),
  
  c(15,360,7),
  c(16,256,8),
  c(17,306,8),
  c(18,630,9),
  c(19,1444,9),
  c(20,1024,10)
)))

tab2 <- tab1 %>%transmute(n2=X1,k2=X2,y2=X3)
tab1 <- tab1 %>%transmute(n1=X1,k1=X2,y1=X3)

tab <-  merge(tab1,tab2,all=TRUE)
# 
# tab <- data.frame()
# x=20
# for (n in 2:x) {
#   for (m in 2:x) {
#     tab <- rbind(tab,expand.grid(n=n,m=m,k=foundEss[m],y=suppsize[n],z=foundEss[n]))
#   }
# }
 tab <- tab %>%
   mutate(result=(k1*k2^y1)^(1/(n1*n2)), neues.n=n1*n2, anzahl.ess=k1*k2^y1)
 
 tab <- tab %>% arrange(desc(result))