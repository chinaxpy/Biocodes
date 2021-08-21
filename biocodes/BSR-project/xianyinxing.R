library("vcfR")
setwd("H:/")
vcf <- read.vcfR("./mk.vcf")
gt <- extract.gt(vcf)
ad <- extract.gt(vcf, "AD")
head(gt)
head(ad)
mask0 <- which( gt[,"pa1K"] == "1/1"  &  gt[,"pa2G"] == "0/0" & gt[,"po6G"] != "NA" & gt[,"po6K"] != "NA" ) 
ad_flt0 <- as.matrix(ad[mask0,c("po6G", "po6K")])
colnames(ad_flt0) <- c("A_Pool","a_Pool")


ED_list0 <- apply(ad_flt0, 1, function(x){

  count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
  depth1 <- count[1] + count[2]
  depth2 <- count[3] + count[4]
  ED_0 <- count[2] / depth1
  ED_1 <- count[4] / depth2
  ED =  ED_1 -ED_0  
  return(ED)
})



mask1 <- which( gt[,"pa1K"] == "0/0"  &  gt[,"pa2G"] == "1/1" & gt[,"po6G"] != "NA" & gt[,"po6K"] != "NA" ) 
ad_flt1 <- as.matrix(ad[mask1,c("po6G", "po6K")])
colnames(ad_flt1) <- c("A_Pool","a_Pool")


ED_list1 <- apply(ad_flt1, 1, function(x){
  
  count <- as.numeric(unlist(strsplit(x, ",",fixed = TRUE,useBytes = TRUE)))
  depth1 <- count[1] + count[2]
  depth2 <- count[3] + count[4]
  ED_0 <- count[1] / depth1
  ED_1 <- count[3] / depth2
  ED =  ED_1 -ED_0  
  return(ED)
})

ED_list <- c(ED_list0,ED_list1)



par(mfrow = c(3,7))

for (i in c("Chr1A","Chr2A","Chr3A","Chr4A","Chr5A","Chr6A","Chr7A",
            "Chr1B","Chr2B","Chr3B","Chr4B","Chr5B","Chr6B","Chr7B",
            "Chr1D","Chr2D","Chr3D","Chr4D","Chr5D","Chr6D","Chr7D") ){
  # i<-"Chr3A"
  ED_flt <- ED_list[grepl(i,names(ED_list))]
  pos <- as.numeric(substring(names(ED_flt), 7))
  plot(pos, ED_flt,
       ylim=c(-1,1),
       pch = 20, cex = 0.2,
       xlab = i,
       ylab = "delt SNP index")
}


#####



calcValueByWindow <- function(pos, value,
                              window_size = 2000000,
                              step_size =   1000000){
  # get the max position in the postion
  max_pos <- max(pos)
  
  # construct the window
  window_start <- seq(0, max_pos + window_size,step_size)
  window_end <- window_start + window_size
  mean_value <- vector(mode = "numeric", length = length(window_start))
  
  # select the value inside the window
  for (j in seq_along(window_start)){
    
    pos_in_window <- which(pos > window_start[j] &   pos < window_end[j])
    value_in_window <- value[pos_in_window]
    mean_value[j] <- mean(value_in_window)
    
  }
  # remove the Not A Number position
  nan_pos <-  is.nan(mean_value)
  mean_value <- mean_value[! nan_pos]
  window_pos <- ((window_start + window_end)/ 2)[!nan_pos]
  df <- data.frame(pos   = window_pos,
                   value = mean_value)
  return(df)
}



#################
par(mfrow = c(3,7))

for (i in c("Chr1A","Chr2A","Chr3A","Chr4A","Chr5A","Chr6A","Chr7A",
            "Chr1B","Chr2B","Chr3B","Chr4B","Chr5B","Chr6B","Chr7B",
            "Chr1D","Chr2D","Chr3D","Chr4D","Chr5D","Chr6D","Chr7D") ){  
  # i<-"chr1A"
  ED_flt <- ED_list[grepl(i,names(ED_list))]
  pos <- as.numeric(substring(names(ED_flt), 7))
  # bin
  df <- calcValueByWindow(pos = pos, value =  ED_flt )
  
  plot(x = pos, y =ED_flt,
       ylim = c(-1,1),
       pch = 20, cex = 0.2,
       xlab = i,
       ylab = ("SNP index"))
  lines(x = df$pos, y = df$value, col = "red")
}

##################



for (i in c("Chr1A","Chr2A","Chr3A","Chr4A","Chr5A","Chr6A","Chr7A",
            "Chr1B","Chr2B","Chr3B","Chr4B","Chr5B","Chr6B","Chr7B",
            "Chr1D","Chr2D","Chr3D","Chr4D","Chr5D","Chr6D","Chr7D") ){  
  # i<-"chr1A"
  ED_flt <- ED_list[grepl(i,names(ED_list))]
  pos <- as.numeric(substring(names(ED_flt), 7))
  df <- calcValueByWindow(pos = pos, value =  ED_flt )
  
  plot(x = pos/1e6, y =ED_flt,
       ylim = c(-1,1),
       pch = 20, cex = 0.2,
       xlab = i,
       ylab = ("SNP index"))
  lines(x = df$pos/1e6, y = df$value, col = "#d7191c") # window红色
  # abline(lm(df$value ~ df$pos), lwd=4, col="blue",type = "s")  #  直线
  # lines(predict(loess(df$value ~ df$pos)),lwd=4, col="blue")
  abline(h=0.8,col="#1a9641") #绿色
  abline(h=0.5,col="#fdae61") #橘黄色
}


############
ED_flt1 <- ED_list[grepl("Chr",names(ED_list))]
pos1 <- as.numeric(substring(names(ED_flt1), 7))
chr1 <- substring(names(ED_flt1), 1,5)
data1<- data.frame(pos= pos1,ED_flt=ED_flt1,chr=chr1)


####
df_window <- data.frame(pos="",value="",chr="")
for (i in c("Chr1A","Chr2A","Chr3A","Chr4A","Chr5A","Chr6A","Chr7A",
            "Chr1B","Chr2B","Chr3B","Chr4B","Chr5B","Chr6B","Chr7B",
            "Chr1D","Chr2D","Chr3D","Chr4D","Chr5D","Chr6D","Chr7D") ){  
  # i<-"chr1A"
  ED_flt_chr <- ED_list[grepl(i,names(ED_list))]
  pos <- as.numeric(substring(names(ED_flt_chr), 7))
  df <- calcValueByWindow(pos = pos, value =  ED_flt_chr )
  df$chr <- rep(i,length(df$pos))
  df_window <- rbind(df_window,df)
}
df_window<- df_window[df_window$chr!="",]
df_window$value <- as.numeric(df_window$value) 
df_window$pos <- as.numeric(df_window$pos) 
###
library(dplyr)
data2 <- full_join(data1,df_window)
data2 <- data2 %>% filter(chr != "ChrUn")
library(ggplot2)
p<- ggplot(data2,aes(pos/1e6, ED_flt) )+ 
  geom_point( cex = 0.2,color="blue") +
  labs(y='Delt SNP index',x='Chrome (MB) ')+
  #geom_smooth(method = "loess", se = T,color = "red",size = 0.5) +
  geom_path(aes(y = value),color = "red",size = 0.5) +
  facet_wrap(. ~ chr,scales="free_x", ncol = 6)+
  # geom_hline(aes(yintercept=0.75), colour="#1a9641", linetype="dashed",size=0.25)+
  geom_hline(aes(yintercept=0.5), colour="#fdae61", linetype="dashed",size=0.25)

p

################

