#add noise
xloc=1:2048
SNR = c(-40,-25,-10)
snr = 10^ (0.1* SNR);
noise <- rnorm(2048)*snr
b_noise=array(0,c(16,16,length(xloc)))
c_noise=array(0,c(16,16,length(xloc)))
for (i in 1:16)
  for (j in 1:16){
    mymean<-mean(b_inputValue[i,i,xloc])
    mysd<-sd(b_inputValue[i,i,xloc])
    b_noise[i,j,xloc]= b_inputValue[i,i,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*b_inputValue[i,j,xloc]
    mymean<-mean(b_inputValue[i,i,xloc])
    mysd<-sd(c_inputValue[i,i,xloc])
    c_noise[i,j,xloc]= c_inputValue[i,i,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*c_inputValue[i,j,xloc]
  }

mymean<-mean(b_inputValue[13,13,xloc/2])
mysd<-sd(b_inputValue[13,13,xloc/2])
b_1313_noise= b_inputValue[13,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*b_inputValue[13,13,xloc]
c_1313_noise= c_inputValue[13,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*c_inputValue[13,13,xloc]

mymean<-mean(b_inputValue[9,13,xloc/2])
mysd<-sd(b_inputValue[9,13,xloc/2])
b_913_noise= b_inputValue[9,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*b_inputValue[9,13,xloc]
c_913_noise= c_inputValue[9,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*c_inputValue[9,13,xloc]

mymean<-mean(b_inputValue[5,13,xloc/2])
mysd<-sd(b_inputValue[5,13,xloc/2])
b_513_noise= b_inputValue[5,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*b_inputValue[5,13,xloc]
c_513_noise= c_inputValue[5,13,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*c_inputValue[5,13,xloc]

mymean<-mean(b_inputValue[1,9,xloc/2])
mysd<-sd(b_inputValue[1,9,xloc/2])
b_19_noise= b_inputValue[1,9,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*b_inputValue[1,9,xloc]
c_19_noise= c_inputValue[1,9,xloc]+rnorm(length(xloc), mean = 0, sd = 1)*snr[3]*c_inputValue[1,9,xloc]
