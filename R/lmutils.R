#
### Simple utility functions for building linear models.
#
lmforwardsequentialsidak <- function(y,object,alpha=0.05) {
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  while(TRUE) {
    pval <- rep(NA,ncol(object$u))
    lm1 <- lm(as.formula(paste("y~",if(length(included)) paste(paste(paste("V_",included,sep=""),collapse="+"),sep="") else "1",sep="")), data=object)
    for(i in candidates) {
      lm2 <- lm(as.formula(paste("y~",if(length(included)) paste(paste(paste("V_",included,sep=""),collapse="+"),"+",sep="") else "","V_",i,sep="")), data=object)
      aovcomp <- anova(lm1,lm2)
      pval[i] <- 1-(1-aovcomp[["Pr(>F)"]][2L])^(length(candidates)-length(included))
    }
    if(min(pval,na.rm=TRUE) < alpha) {
      included <- c(included,candidates[candidates==which.min(pval)])
      candidates <- candidates[candidates!=which.min(pval)]
    } else {
      aovlm1 <- anova(lm1)
      lm1[["SidakP"]] <- c(NA,1-(1-aovlm1[["Pr(>F)"]][-c(1,length(aovlm1$Df))])^(ncol(object$u):(ncol(object$u)-sum(aovlm1$Df[-length(aovlm1$Df)])+2L)))
      return(lm1)
    }
  }
}
#
lmforwardsequentialAICc <- function(y,object) {
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  while(TRUE) {
    AICc2 <- rep(NA,ncol(object$u))
    lm1 <- lm(as.formula(paste("y~",if(length(included)) paste(paste(paste("V_",included,sep=""),collapse="+"),sep="") else "1",sep="")), data=object)
    k1 <- length(lm1$coef) ; AICc1 <- AIC(lm1) + (2*k1*(k1+1)/(length(y)-k1-1))
    for(i in candidates) {
      lm2 <- lm(as.formula(paste("y~",if(length(included)) paste(paste(paste("V_",included,sep=""),collapse="+"),"+",sep="") else "","V_",i,sep="")), data=object)
      k2 <- length(lm2$coef) ; AICc2[i] <- AIC(lm2) + (2*k2*(k2+1)/(length(y)-k2-1))
    }
    if(min(AICc2,na.rm=TRUE) < AICc1) {
      included <- c(included,candidates[candidates==which.min(AICc2)])
      candidates <- candidates[candidates!=which.min(AICc2)]
    } else {
      lm1$AICc <- AICc1
      return(lm1)
    }
  }
}
#
