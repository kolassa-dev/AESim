#' @title Simulating correllated efficacy and severity data
#' @description \code{onesample} generates a data frame containing one efficacy variable and a number of adverse event duration and severity variables.
#' @param nn The number of samples.
#' @param eff A cutoff creating a binary efficacy variable from a convolution of exponentials.
#' @param aeprob A table whose columns represent adverse events, and whose rows represent probabilities of various severity levels.
#' @param alpha a multiplier for the shared component between efficacy, severity, duration.
#' @param durexp Overall duration expectation
#' @param severrel a multiplier for the shared component between severity, duration.
#' @param override If true fixes expected turation to 8.
#' @return A data frame with the above columns.
#' @examples
#' #An example with adverse event duration 8 days, low correlation among
#' #efficacy, severity, and duration.  alpha=1
#' aaa<-onesample(durexp=8)
#' #An example with adverse event duration 8 days, high correlation among
#' #efficacy, severity, and duration.
#' aaa<-onesample(durexp=8,alpha=7)
#' #An example with adverse event duration 8 days, no correlation among
#' #efficacy, severity, and duration, extra correlation between severity and
#' #duration.  Setting severrel=1 isn't supported yet.
#' aaa<-onesample(durexp=8,severrel=1,alpha=0)
#' @importFrom stats rexp
#' @export
onesample<-function(nn=200,eff=.4,
   aeprob=cbind( c(.80,.05,.05,.05,.05), c(.85,.05,.05,.03,.02),
      c(.90,.05,.03,.02,.01), c(.95,.02,.01,.01,.01)),alpha=1,durexp=NULL,
     severrel=0,override=FALSE){
   if(is.null(durexp)){
      durexp<-rep(1,dim(aeprob)[1])
   }
   if(length(durexp)==1){
      durexp<-rep(durexp,dim(aeprob)[1])
   }
   myp<-function(z,av){
      av<-rev(sort(av))
#     alpha<-av[1]; beta<-av[2]; gamma<-av[3]
      if((av[3]==av[2])&(av[2]!=0)) message("hurt 2")
      if(av[1]==av[2]){
         if(av[3]==0){
            out <- 1 - (av[1] + z)/(av[1]*exp(z/av[1]))
         }else{
            out <- (av[1]^ 2*(1 - exp(-z/av[1])) - 2*av[1]*av[3] - 
            (av[1]*(-2*av[3] + z))/exp(z/av[1]) + 
            av[3]*(av[3] - av[3]/exp(z/av[3]) + z/exp(z/av[1])))/ (av[1] - av[3])^2
         }
      }else{
         if(av[3]==0){
            if(av[2]==0){
               out <- 1-exp(-z/av[1])
            }else{
               out <- (av[1] - av[1]/exp(z/av[1]) + av[2]*(-1 + exp(-z/av[2])))/(av[1] - av[2])
            }
         }else{
            out<-(av[2] - av[2]/exp(z/av[1]) + 
               (av[2]^2*(exp(-z/av[1]) - exp(-z/av[2])))/(-av[1] + av[2]) - 
               (av[3]*(av[1] - av[1]/exp(z/av[1]) + (-1 + exp(-z/av[3]))*av[3]))/
               (av[1] - av[3]))/(av[2] - av[3])
         }
      }
      return(out)
   }
# Extract the number of columlns from the third argument, giving the number
# of adverse events.
   nae<-dim(aeprob)[2]
# Create data frames to contain severity (aes) and duration (aed) values.
   aes<-aed<-data.frame(array(NA,c(nn,nae)))
# Give names to these columns.
   names(aes)<-paste("Sev",seq(nae),sep="")
   names(aed)<-paste("Dur",seq(nae),sep="")
# Create a common variable contributing to efficacy, AE serverity, and duration.
   bad<-alpha*rexp(nn)## To adjust correlation, multiply this by a constant.  Sepcifically, to lower correlation, multiply by something less than one.
   extra<-severrel*rexp(nn)
# Create a continuous efficacy variable, consisting of the shared component calculated above, and a new independent component.  This is now Gamma(2).
   ra1<-rexp(nn)+bad
# Convert this efficacy variable into a binary variable, by cutting off by the quantile of the gamma(2) distributon.
#  cat("About to call myp with paramters",c(alpha,1,0),"\n")
   ef1<-as.numeric(myp(ra1,c(alpha,1,0))>(1-eff)) ## As written, eff reflects the probability of a positive efficacy outcome.
#  browser()
   temp<-data.frame(array(NA,c(nn,nae)))
   names(temp)<-paste("T",seq(nae),sep="")
   for(jj in seq(nae)){ 
# Generate an independent component for each ae, and convert to a probability.
#     cat("About to call myp with paramters",c(alpha,1,severrel),"\n")
      temp[,jj]<-myp(rexp(nn)+bad+extra,c(alpha,1,severrel))
# Use this to select among ordered categories.
      aes[,jj]<-apply(outer(cumsum(aeprob[,jj]),temp[,jj],"<"),2,"sum")
# Generate duration for AE.
      aed[,jj]<-if(override) rexp(nn)+7*rexp(nn) else bad+extra+(durexp[jj]-alpha-severrel)*rexp(nn)
   }
   return(data.frame(cbind(ef1,aes,aed)))
}
#set.seed(5)
#An example with adverse event duration 8 days, low correlation among
#efficacy, severity, and duration.  alpha=1
#aaa<-onesample(durexp=8)
#An example with adverse event duration 8 days, correlation among severity and duration but not efficacy.
#aaa<-onesample(durexp=8,alpha=0,severrel=1)
#An example with adverse event duration 8 days, high correlation among
#efficacy, severity, and duration.
#aaa<-onesample(durexp=8,alpha=7)
#An example with adverse event duration 8 days, no correlation among
#efficacy, severity, and duration, extra correlation between severity and
#duration.  Setting severrel=1 isn't supported yet.
#aaa<-onesample(durexp=8,severrel=2,alpha=0)
#' @title Adding a new efficacy  variable.
#' @description \code{addeff} generates a second efficacy variable dependent on first.
#' @param ds input data set
#' @param pv response probabilities for initial nonresponders and initial responders.
#' @return A data frame like input with an additional efficacy value.
#' @examples
#' #Create a treatment set
#' ct<-onesample(eff=0.10,aeprob=cbind(c(0.965,0.020,0.015,0,0),
#'    c(0.95,0.03,0.02,0,0), c(0.98,0.01,0.01,0,0)))
#' ct$group<-0
#' ct<-addeff(ct,c(0.07,0.90))
#' tr<-onesample(eff=0.36,aeprob=cbind(c(0.945,0.030,0.020,0.005,0),
#'    c(0.90,0.05,0.04,0.01,0),c(0.90,0.05,0.04,00.01,0)))
#' tr$group<-1
#' tr<-addeff(tr,c(0.18,0.90))
#' all<-rbind(ct,tr)
#' @importFrom stats rbinom
#' @export
addeff<-function(ds,pv){
   ds$ef2<-NA
   ds$ef2[ds$ef1==1]<-rbinom(sum(ds$ef1==1),prob=pv[2],size=1)
   ds$ef2[ds$ef1==0]<-rbinom(sum(ds$ef1==0),prob=pv[1],size=1)
   return(ds)
}







