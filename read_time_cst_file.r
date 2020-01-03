-library(EMD)
library(FWHVA)
library(igraph)
library(latex2exp)
library(infotheo)
library(mpmi)

MyPath="G:/eresearch/tmp/dipole/"

#bleedingdata<- read.delim(file=paste(MyPath,"110411_bleeding_zero_0.txt",sep=""), header=TRUE, sep="\t")

Status<-0

#conn <- file(paste(MyPath,"110411_bleeding_zero_125.txt",sep=""))
#conn <- file("data/110411_bleeding_zero_0_4.txt")
inputDegreeSeq=array(0,c(16,16,2048))
inputValue=array(0,c(16,16,2048))
b_inputDegreeSeq=array(0,c(16,16,2048))
b_inputValue=array(0,c(16,16,2048))
c_inputDegreeSeq=array(0,c(16,16,2048))
c_inputValue=array(0,c(16,16,2048))

for (myk in 1:2){
    if (myk==1)
      conn <- file("data/zuba_bleeding_3.txt")
      #conn <- file("data/p111716_bleeding2.txt")
    else 
      conn <- file("data/zuba_clot_3.txt")
      #conn <- file("data/p111716_clot2.txt")

    bleedingdata <-readLines(conn)
    close(conn)
    
    myLine<-0;
    #S_Len<-107371-101059
    S_Len<-79491-74819
    
    Status<-0
    O1_1=rep(0,16)
    O2_1=rep(0,16)
    O_loc=array(0,c(16,16))
    Oij=array(0,c(16,16,S_Len))
    First_line<-1# 74819  #o11,79495  #776219

    k61<-1
    InStr1<-"o2,1/"
    InStr61<-"o6,1/"
    InStr<-"o1,1/"
    k<-1
    k21<-1
    j<-1
    for (i in  1:(length(bleedingdata)) ) {
      temp<-grep(InStr,bleedingdata[i])
      if (length(temp)>=1) {
        O_loc[(j),k]<-i;
        #          cat (InStr,"\n")
        #InStr<-paste0("o",((j) %% 16)+1 ,',',(k),'/')
        j<-j+1
        InStr<-paste0("o",k ,',',(j),'/')
        if (j>16) {
          j=1
          k<-k+1
          #InStr<-paste0("o",(j %% 16)+1 ,',',(k),'/')
          InStr<-paste0("o",k ,',',j,'/')
          #            cat (InStr,"--\n")
        }
        k21<-k21+1
      }        
    }
    
    if (k21<256){
      k<-1
      j<-1
      k21<-1
      InStr<-"o1,1/"
      for (i in  1:(length(bleedingdata)) ) {
        temp<-grep(InStr,bleedingdata[i])
        if (length(temp)>=1) {
          O_loc[(j),k]<-i;
          cat (InStr,"\n")
          #InStr<-paste0("o",((j) %% 16)+1 ,',',(k),'/')
          j<-j+1
          InStr<-paste0("o",j ,',',(k),'/')
          if (j>16) {
            j=1
            k<-k+1
            #InStr<-paste0("o",(j %% 16)+1 ,',',(k),'/')
            InStr<-paste0("o",j ,',',k,'/')
            cat (InStr,"--\n")
          }
          k21<-k21+1
        }   
        
      }
    }
    

    Data_Skip<-100
    S_Len<-O_loc[2,1]-O_loc[1,1]-80-Data_Skip #4
    C_O1_1=array(0,c(16,16))
    C_MO1_1=array(0,c(16,16))
    #B_O1_1=array(0,c(16,16))
    #B_MO1_1=array(0,c(16,16))
    #N_O1_1=array(0,c(16,16))
    #N_MO1_1=array(0,c(16,16))
    #B_O1_1=array(0,c(16,16))
    #B_MO1_1=array(0,c(16,16))
    
    
    
    #First_line<-101059  #113691(3)  ## 138955   (1)
    #par(mfrow=c(7,1), mar=c(2,1,2,1))
    for (k in 1:nrow(O_loc))  
      for (j in 1:nrow(O_loc)) {
        #for (k in 1:length(O_loc))  
        #for (k in 1:16) {
        #First_line<-O1_1[k]+2220
        First_line<-O_loc[k,j]+Data_Skip
        in_x<-1
        Input_T<-rep(0,S_Len)
        Input_D<-rep(0,S_Len)
        for (i in First_line:(First_line+S_Len)) {
          #print(linn[i])
          mystr<-unlist(strsplit(bleedingdata[i], "    "))
          loc<-1
          for (x in mystr){
            if (nchar(x)>3) {
              if (loc==1) Input_T[in_x]=as.numeric(x)       
              else if (loc==2) Input_D[in_x]=as.numeric(x)       
              loc<-loc+1
            }
          }
          #cat(Input_T[in_x],"\t",Input_D[in_x],"\n")
          in_x=in_x+1
        }
        ind1<-max(which(Input_T<3))  #3
        ind2<-max(which(Input_T<12))  #12
        inputData<-Input_D[ind1:ind2]
        #approx(x, inputData, method = "constant")
        x <- Input_T[ind1:ind2]
        outy<-approx(x, inputData, method = "constant",n=2048)      
        if (myk==1)
          b_inputValue[k,j,]=outy$y
        else
          c_inputValue[k,j,]=outy$y
        
        adj2=getHVGAdjMatrix(outy$y,1)
        adj2[1,1]=0
        adj2[2,2]=0
        adj2[3,3]=0
        g2 <- graph_from_adjacency_matrix(adj2, mode="undirected", weighted=NULL)
        #      plot(Input_T, Input_D, type="l")
        myDegree<-degree(g2)
        if (myk==1)
          c_inputDegreeSeq[k,j,]<-myDegree
        else
          b_inputDegreeSeq[k,j,]<-myDegree
        #      plot(myDegree,type='l')
        C_O1_1[k,j]<-max(myDegree)
        C_MO1_1[k,j]<-mean(myDegree)      
        #     AO1_1[k]<-mean(myDegree)
      }
}

