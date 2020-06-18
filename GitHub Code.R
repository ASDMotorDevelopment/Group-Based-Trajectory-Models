#Add the required packages
library(lcmm)
library(knitr)

#The following code comes from LCTMtools package
#See references for details
LCTMtoolkit <- function(model) {
  n <- nrow(model$pprob)
  K <- ncol(model$pprob) - 2
  p <- model$pprob[, c(paste0("prob", 1:K, sep = ""))]
  
  
  if (class(model$call) == "SAS") {
    PI <- os$PI/100
  } else {
    PI <- exp(c(model$best[1:(K - 1)], 0))/(sum(exp(c(model$best[1:(K - 1)],
                                                      0))))
  }
  
  
  outputs <- matrix(0, nrow = 3, ncol = K)
  colnames(outputs) <- paste0("Class_", 1:K, sep = "")
  rownames(outputs) <- c("APPA", "OCC", "Mismatch")
  outputs[1, ] <- appa(p)
  outputs[2, ] <- occ(p, PI)
  outputs[3, ] <- mismatch(p, PI)
  Recommendation <- c("Greater than 0.7", "Greater than 5", "Close to zero")
  outputs <- data.frame(round(outputs, 3), Recommendation)
  
  ep  <- entropy(p)
  rep <- relative_entropy(p)
  outputs1 <- t(data.frame(Entropy = ep,
                           Relative_entropy = rep,
                           BIC = model$BIC,
                           AIC = model$AIC))
  Recommendation <- c("Close to zero", "Close to 1", "-", "-")
  outputs1 <- round(outputs1, 3)
  outputs1 <- data.frame(Model = outputs1, Recommendation = Recommendation)
  
  outputs2 <- list(outputs, outputs1)
  names(outputs2) <- c("Class-specific", "Model-specific")
  
  print(outputs2)
  
  outputs3 <- list(appa=outputs[1,1:K],
                   occ=outputs[2,1:K],
                   mismatch=outputs[3,1:K],
                   entropy=ep,
                   relativeentropy=rep,
                   BIC=model$BIC,
                   AIC=model$AIC
  )
  return(outputs3)
}
#occ
occ <- function(p, pi) {
  app <- appa(p)
  numerator <- app/(1 - app)
  denominator <- pi/(1 - pi)
  occ <- numerator/denominator
  return(occ)
}
#Appa
appa <- function(p) {
  # determine class size
  K <- ncol(p)
  
  # determine class
  group <- class_assignment(p)
  
  # save vector for appa values
  app <- rep(NA, times = K)
  
  # Compute average posterior probabilites
  for (i in 1:K) {
    classp <- p[group == i, i]
    if (length(classp) != 0) 
      app[i] <- mean(classp)
  }
  return(app)
}
#mismatch
mismatch <- function(p, pi) {
  # determine class
  group <- class_assignment(p)
  K <- ncol(p)
  
  # Tabulate the actual proportions
  actprop <- tabulate(group, nbins = K)/nrow(p)
  
  # Compute the mismatch (mismatch=actual-estimated)
  return(actprop - pi)
}
#entropy
entropy <- function(p) {
  ent <- -1 * sum(p * log(p), na.rm = TRUE)
  
  return(ent)
}
#relative entropy
relative_entropy <- function(p) {
  K <- ncol(p)
  n <- nrow(p)
  
  relEntropy <- 1 + (sum(p * log(p), na.rm = TRUE)/(n * log(K)))
  
  return(relEntropy)
}
#Class assignmnet
class_assignment <- function(p) {
  return(as.numeric(apply(p, 1, which.max)))
}

#######

#Start with the analysis of gross motor scores
#Model A: Linear
#k=1
modelA.1=hlme(fixed = Y ~ 1 + X,
              random = ~-1,
              ng = 1, 
              data = data1, 
              subject = "ID")
linA=c("Linear",1,modelA.1$BIC,NA, modelA.1$conv)
#k=2
modelA.2=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 2, 
              data = data1, 
              subject = "ID")
toolkitA.2=LCTMtoolkit(modelA.2)
linA=rbind(linA,c("Linear",2,modelA.2$BIC,toolkitA.2$relativeentropy, modelA.2$conv))
#k=3
modelA.3=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 3, 
              data = data1, 
              subject = "ID")
toolkitA.3=LCTMtoolkit(modelA.3)
linA=rbind(linA, c("Linear",3, modelA.3$BIC,toolkitA.3$relativeentropy, modelA.3$conv))
#k=4
modelA.4=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 4, 
              data = data1, 
              subject = "ID")
toolkitA.4=LCTMtoolkit(modelA.4)
linA=rbind(linA, c("Linear",4, modelA.4$BIC,toolkitA.4$relativeentropy, modelA.4$conv))
#k=5
modelA.5=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 5, 
              data = data1, 
              subject = "ID")
toolkitA.5=LCTMtoolkit(modelA.5)
linA=rbind(linA, c("Linear",5, modelA.5$BIC,toolkitA.5$relativeentropy, modelA.5$conv))
#k=6
modelA.6=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 6, 
              data = data1, 
              subject = "ID")
toolkitA.6=LCTMtoolkit(modelA.6)
linA=rbind(linA, c("Linear",6, modelA.6$BIC,toolkitA.6$relativeentropy, modelA.6$conv))
#k=7
modelA.7=hlme(fixed = Y ~ 1 + X,
              mixture = ~1 + X,
              random = ~-1,
              ng = 7, 
              data = data1, 
              subject = "ID")
toolkitA.7=LCTMtoolkit(modelA.7)
linA=rbind(linA, c("Linear",7, modelA.7$BIC,toolkitA.7$relativeentropy, modelA.7$conv))

#Model B: Quadratic
#k=1
modelB.1=hlme(fixed = Y ~ 1 + X + I(X^2),
              random = ~-1,
              ng = 1, 
              data = data1, 
              subject = "ID")
linB=c("Quadratic",1,modelB.1$BIC,NA, modelB.1$conv)
#k=2
modelB.2=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 2, 
              data = data1, 
              subject = "ID")
toolkitB.2=LCTMtoolkit(modelB.2)
linB=rbind(linB,c("Quadratic",2,modelB.2$BIC,toolkitB.2$relativeentropy, modelB.2$conv))
#k=3
modelB.3=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 3, 
              data = data1, 
              subject = "ID")
toolkitB.3=LCTMtoolkit(modelB.3)
linB=rbind(linB, c("Quadratic",3,modelB.3$BIC,toolkitB.3$relativeentropy, modelB.3$conv))
#k=4
modelB.4=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 4, 
              data = data1, 
              subject = "ID")
toolkitB.4=LCTMtoolkit(modelB.4)
linB=rbind(linB, c("Quadratic",4,modelB.4$BIC,toolkitB.4$relativeentropy, modelB.4$conv))
#k=5
modelB.5=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 5, 
              data = data1, 
              subject = "ID")
toolkitB.5=LCTMtoolkit(modelB.5)
linB=rbind(linB, c("Quadratic",5,modelB.5$BIC,toolkitB.5$relativeentropy, modelB.5$conv))
#k=6
modelB.6=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 6, 
              data = data1, 
              subject = "ID")
toolkitB.6=LCTMtoolkit(modelB.6)
linB=rbind(linB, c("Quadratic",6,modelB.6$BIC,toolkitB.6$relativeentropy, modelB.6$conv))
#k=7
modelB.7=hlme(fixed = Y ~ 1 + X + I(X^2),
              mixture = ~1 + X + I(X^2),
              random = ~-1,
              ng = 7, 
              data = data1, 
              subject = "ID")
toolkitB.7=LCTMtoolkit(modelB.7)
linB=rbind(linB, c("Quadratic",7,modelB.7$BIC,toolkitB.7$relativeentropy, modelB.7$conv))

#Output for all models
kable(rbind(linA,linB), col.names = c("Model Type","Number of Groups", "BIC","Entropy", "Convergence"), row.names = FALSE, align = "c")

#Now look at the output for the choosen model
#APPA and OCC
LCTMtoolkit(modelB.4)
#Group Classification
postprob(modelB.4)
#Model Summary
summary(modelB.4)

#######
#Repeat for the analysis of fine motor scores
#Model A: Linear
#k=1
mA.1=hlme(fixed = Z ~ 1 + X,
          random = ~-1,
          ng = 1, 
          data = data1, 
          subject = "ID")
linA=c("Linear",1,mA.1$BIC,NA, mA.1$conv)
#k=2
mA.2=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 2, 
          data = data1, 
          subject = "ID")
toolkitA.2=LCTMtoolkit(mA.2)
linA=rbind(linA,c("Linear",2,mA.2$BIC,toolkitA.2$relativeentropy, mA.2$conv))
#k=3
mA.3=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 3, 
          data = data1, 
          subject = "ID")
toolkitA.3=LCTMtoolkit(mA.3)
linA=rbind(linA, c("Linear",3, mA.3$BIC,toolkitA.3$relativeentropy, mA.3$conv))
#k=4
mA.4=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 4, 
          data = data1, 
          subject = "ID")
toolkitA.4=LCTMtoolkit(mA.4)
linA=rbind(linA, c("Linear",4, mA.4$BIC,toolkitA.4$relativeentropy, mA.4$conv))
#k=5
mA.5=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 5, 
          data = data1, 
          subject = "ID")
toolkitA.5=LCTMtoolkit(mA.5)
linA=rbind(linA, c("Linear",5, mA.5$BIC,toolkitA.5$relativeentropy, mA.5$conv))
#k=6
mA.6=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 6, 
          data = data1, 
          subject = "ID")
toolkitA.6=LCTMtoolkit(mA.6)
linA=rbind(linA, c("Linear",6, mA.6$BIC,toolkitA.6$relativeentropy, mA.6$conv))
#k=7
mA.7=hlme(fixed = Z ~ 1 + X,
          mixture = ~1 + X,
          random = ~-1,
          ng = 7, 
          data = data1, 
          subject = "ID")
toolkitA.7=LCTMtoolkit(mA.7)
linA=rbind(linA, c("Linear",7, mA.7$BIC,toolkitA.7$relativeentropy, mA.7$conv))

#Model B: Quadratic
#k=1
mB.1=hlme(fixed = Z ~ 1 + X + I(X^2),
          random = ~-1,
          ng = 1, 
          data = data1, 
          subject = "ID")
linB=c("Quadratic",1,mB.1$BIC,NA, mB.1$conv)
#k=2
mB.2=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 2, 
          data = data1, 
          subject = "ID")
toolkitB.2=LCTMtoolkit(mB.2)
linB=rbind(linB,c("Quadratic",2,mB.2$BIC,toolkitB.2$relativeentropy, mB.2$conv))
#k=3
mB.3=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 3, 
          data = data1, 
          subject = "ID")
toolkitB.3=LCTMtoolkit(mB.3)
linB=rbind(linB, c("Quadratic",3,mB.3$BIC,toolkitB.3$relativeentropy, mB.3$conv))
#k=4
mB.4=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 4, 
          data = data1, 
          subject = "ID")
toolkitB.4=LCTMtoolkit(mB.4)
linB=rbind(linB, c("Quadratic",4,mB.4$BIC,toolkitB.4$relativeentropy, mB.4$conv))
#k=5
mB.5=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 5, 
          data = data1, 
          subject = "ID")
toolkitB.5=LCTMtoolkit(mB.5)
linB=rbind(linB, c("Quadratic",5,mB.5$BIC,toolkitB.5$relativeentropy, mB.5$conv))
#k=6
mB.6=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 6, 
          data = data1, 
          subject = "ID")
toolkitB.6=LCTMtoolkit(mB.6)
linB=rbind(linB, c("Quadratic",6,mB.6$BIC,toolkitB.6$relativeentropy, mB.6$conv))
#k=7
mB.7=hlme(fixed = Z ~ 1 + X + I(X^2),
          mixture = ~1 + X + I(X^2),
          random = ~-1,
          ng = 7, 
          data = data1, 
          subject = "ID")
toolkitB.7=LCTMtoolkit(mB.7)
linB=rbind(linB, c("Quadratic",7,mB.7$BIC,toolkitB.7$relativeentropy, mB.7$conv))

#Model C: Cubic
#k=1
mC.1=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 1, 
          data = data1, 
          subject = "ID")
linC=c("Cubic",1,mC.1$BIC,NA, mC.1$conv)
#k=2
mC.2=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 2, 
          data = data1, 
          subject = "ID")
toolkitC.2=LCTMtoolkit(mC.2)
linC=rbind(linC,c("Cubic",2,mC.2$BIC,toolkitC.2$relativeentropy, mC.2$conv))
#k=3
mC.3=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 3, 
          data = data1, 
          subject = "ID")
toolkitC.3=LCTMtoolkit(mC.3)
linC=rbind(linC, c("Cubic",3,mC.3$BIC,toolkitC.3$relativeentropy, mC.3$conv))
#k=4
mC.4=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 4, 
          data = data1, 
          subject = "ID")
toolkitC.4=LCTMtoolkit(mC.4)
linC=rbind(linC, c("Cubic",4,mC.4$BIC,toolkitC.4$relativeentropy, mC.4$conv))
#k=5
mC.5=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 5, 
          data = data1, 
          subject = "ID")
toolkitC.5=LCTMtoolkit(mC.5)
linC=rbind(linC, c("Cubic",5,mC.5$BIC,toolkitC.5$relativeentropy, mC.5$conv))
#k=6
mC.6=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2)+ I(X^3),
          random = ~-1,
          ng = 6, 
          data = data1, 
          subject = "ID")
toolkitC.6=LCTMtoolkit(mC.6)
linC=rbind(linC, c("Cubic",6,mC.6$BIC,toolkitC.6$relativeentropy, mC.6$conv))
#k=7
mC.7=hlme(fixed = Z ~ 1 + X + I(X^2) + I(X^3),
          mixture = ~1 + X + I(X^2) + I(X^3),
          random = ~-1,
          ng = 7, 
          data = data1, 
          subject = "ID")
toolkitC.7=LCTMtoolkit(mC.7)
linC=rbind(linC, c("Cubic",7,mC.7$BIC,toolkitC.7$relativeentropy, mC.7$conv))

#Output for all models
kable(rbind(linA,linB,linC), col.names = c("Model Type","Number of Groups", "BIC","Entropy", "Convergence"), row.names = FALSE, align = "c")

#Now look at the output for the choosen model
#APPA and OCC
LCTMtoolkit(mC.4)
#Group Classification
postprob(mC.4)
#Model Summary
summary(mC.4)

#References
#Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of Extended Mixed Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#Journal of Statistical Software, 78(2). https://doi.org/10.18637/jss.v078.i02
