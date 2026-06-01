
# set directory for Design matrix and spiking data
setwd("E:/Analysis_current/01_Behav+imaging/Version_4/DK103/2018-02-22_lobv/GLM/glmRfiles/CorData_small/")

# Read in spiking data
Ycor <- as.matrix(read.csv("YcorSmall.csv", header=FALSE))
YcorBin <- as.matrix((Ycor > 0)+0)

# Read in design matrices (normal and PCA)
DM <- as.matrix(read.csv("DMcorSmall.csv", header=FALSE))

DMpca <- as.matrix(read.csv("DMcorSmallPCAscore.csv", header=FALSE))
nPCs <- as.numeric(read.csv("nPCvar.csv", header = FALSE))
DMnPCs <- DMpca[,1:nPCs]

# Load in xvalfold info
iFold <- as.matrix(read.csv("iFoldCsmall.csv", header = FALSE));
iHold <- as.matrix(which(iFold == max(iFold)))
iTest <- as.matrix(which(iFold < max(iFold)))
nFold <- as.matrix(iFold[iTest])

# Make test and holdout matrices
Ycor_test <- Ycor[iTest,]
Ycor_hold <- Ycor[iHold,]
YcorBin_test <- YcorBin[iTest,]
YcorBin_hold <- YcorBin[iHold,]

DM_test <- DM[iTest,]
DM_hold <- DM[iHold,]
DMnPCs_test <- DMnPCs[iTest,]
DMnPCs_hold <- DMnPCs[iHold,]


# Test good cell (i = 184) - gaussian
i <- 184
cvG1 = cv.glmnet(DM_test,Ycor_test[,i],alpha = 1, foldid = nFold, standardize = FALSE)
cvG.9 = cv.glmnet(DM_test,Ycor_test[,i],alpha = .9, foldid = nFold, standardize = FALSE)
cvG.5 = cv.glmnet(DM_test,Ycor_test[,i],alpha = 0.5, foldid = nFold, standardize = FALSE)

Ycor_pred1 <- predict(cvG1, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.9 <- predict(cvG.9, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.5 <- predict(cvG.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,2))
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred1[1:500], col = 'red')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.9[1:500], col = 'green')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.5[1:500], col = 'blue')


cvGB1 = cv.glmnet(DM_test,YcorBin_test[,i],alpha = 1, foldid = nFold, standardize = FALSE)
cvGB.9 = cv.glmnet(DM_test,YcorBin_test[,i],alpha = .9, foldid = nFold, standardize = FALSE)
cvGB.5 = cv.glmnet(DM_test,YcorBin_test[,i],alpha = 0.5, foldid = nFold, standardize = FALSE)

YcorB_pred1 <- predict(cvGB1, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.9 <- predict(cvGB.9, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.5 <- predict(cvGB.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,3))
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred1[1:500], col = 'red')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.9[1:500], col = 'red')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.5[1:500], col = 'red')
plot(coef(cvGB1 , s = "lambda.min"))
plot(coef(cvGB.9 , s = "lambda.min"))
plot(coef(cvGB.5 , s = "lambda.min"))


cvGnPCsB1 = cv.glmnet(DMnPCs_test,YcorBin_test[,i],alpha = 1, foldid = nFold, standardize = FALSE)
cvGnPCsB.9 = cv.glmnet(DMnPCs_test,YcorBin_test[,i],alpha = .9, foldid = nFold, standardize = FALSE)
cvGnPCsB.5 = cv.glmnet(DMnPCs_test,YcorBin_test[,i],alpha = 0.5, foldid = nFold, standardize = FALSE)

YcorB_pred1 <- predict(cvGnPCsB1, newx = DMnPCs_hold, type = "response", s = "lambda.min")
YcorB_pred.9 <- predict(cvGnPCsB.9, newx = DMnPCs_hold, type = "response", s = "lambda.min")
YcorB_pred.5 <- predict(cvGnPCsB.5, newx = DMnPCs_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,3))
plot(YcorBin_hold[1:100,i],type = 'l')
lines(YcorB_pred1[1:100], col = 'red')
plot(YcorBin_hold[1:100,i],type = 'l')
lines(YcorB_pred.9[1:100], col = 'red')
plot(YcorBin_hold[1:100,i],type = 'l')
lines(YcorB_pred.5[1:100], col = 'red')
plot(coef(cvGnPCsB1 , s = "lambda.min"))
plot(coef(cvGnPCsB.9 , s = "lambda.min"))
plot(coef(cvGnPCsB.5 , s = "lambda.min"))







# Test good cell (i = 184) - Poisson
i <- 184
cvP1 = cv.glmnet(DM_test,Ycor_test[,i], family = "poisson", alpha = 1, foldid = nFold, standardize = FALSE)
cvP.9 = cv.glmnet(DM_test,Ycor_test[,i], family = "poisson", alpha = .9, foldid = nFold, standardize = FALSE)
cvP.5 = cv.glmnet(DM_test,Ycor_test[,i], family = "poisson", alpha = 0.5, foldid = nFold, standardize = FALSE)

Ycor_pred1 <- predict(cvP1, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.9 <- predict(cvP.9, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.5 <- predict(cvP.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,2))
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred1[1:500], col = 'red')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.9[1:500], col = 'green')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.5[1:500], col = 'blue')


cvPB1 = cv.glmnet(DM_test,YcorBin_test[,i], family = "poisson", alpha = 1, foldid = nFold, standardize = FALSE)
cvPB.9 = cv.glmnet(DM_test,YcorBin_test[,i], family = "poisson", alpha = .9, foldid = nFold, standardize = FALSE)
cvPB.5 = cv.glmnet(DM_test,YcorBin_test[,i], family = "poisson", alpha = 0.5, foldid = nFold, standardize = FALSE)

YcorB_pred1 <- predict(cvPB1, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.9 <- predict(cvPB.9, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.5 <- predict(cvPB.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,3))
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred1[1:500], col = 'red')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.9[1:500], col = 'green')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.5[1:500], col = 'blue')
plot(coef(cvPB1 , s = "lambda.min"))
plot(coef(cvPB.9 , s = "lambda.min"))
plot(coef(cvPB.5 , s = "lambda.min"))

# Test good cell (i = 184) - Bernoulli
i <- 184
cvB1 = cv.glmnet(DM_test,Ycor_test[,i], family = "binomial", alpha = 1, foldid = nFold, standardize = FALSE)
cvB.9 = cv.glmnet(DM_test,Ycor_test[,i], family = "binomial", alpha = .9, foldid = nFold, standardize = FALSE)
cvB.5 = cv.glmnet(DM_test,Ycor_test[,i], family = "binomial", alpha = 0.5, foldid = nFold, standardize = FALSE)

Ycor_pred1 <- predict(cvB1, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.9 <- predict(cvB.9, newx = DM_hold, type = "response", s = "lambda.min")
Ycor_pred.5 <- predict(cvB.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,2))
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred1[1:500], col = 'red')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.9[1:500], col = 'green')
plot(Ycor_hold[1:500,i],type = 'l')
lines(Ycor_pred.5[1:500], col = 'blue')


cvBB1 = cv.glmnet(DM_test,YcorBin_test[,i], family = "binomial", alpha = 1, foldid = nFold, standardize = FALSE)
cvBB.9 = cv.glmnet(DM_test,YcorBin_test[,i], family = "binomial", alpha = .9, foldid = nFold, standardize = FALSE)
cvBB.5 = cv.glmnet(DM_test,YcorBin_test[,i], family = "binomial", alpha = 0.5, foldid = nFold, standardize = FALSE)

YcorB_pred1 <- predict(cvBB1, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.9 <- predict(cvBB.9, newx = DM_hold, type = "response", s = "lambda.min")
YcorB_pred.5 <- predict(cvBB.5, newx = DM_hold, type = "response", s = "lambda.min")

par(mfrow = c(2,2))
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred1[1:500], col = 'red')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.9[1:500], col = 'green')
plot(YcorBin_hold[1:500,i],type = 'l')
lines(YcorB_pred.5[1:500], col = 'blue')












for (i in 1:nRois)
{
  cvG <- cv.glmnet(DM_test,5)
  
}


# Run an actual GLM!!
ytest12 <- YtrialTest[,12]
yhold12 <- YtrialHold[,12]

fit12 <- glmnet(DMsmallTest, ytest12, alpha = 0.5)

# test different alphas and plot
cvfit12_1 <- cv.glmnet(DMsmallTest, ytest12, alpha = 1 ,foldid = iFoldTest)
cvfit12_09 <- cv.glmnet(DMsmallTest, ytest12, alpha = 0.9 ,foldid = iFoldTest)
cvfit12_05 <- cv.glmnet(DMsmallTest, ytest12, alpha = 0.5 ,foldid = iFoldTest)
cvfit12_0 <- cv.glmnet(DMsmallTest, ytest12, alpha = 0 ,foldid = iFoldTest)

par(mfrow = c(2,2))
plot(cvfit12_1); plot(cvfit12_09); plot(cvfit12_05); plot(cvfit12_0)

# test poisson glm
fitP12 <- glmnet(DMsmallTest, ytest12, family = "poisson", alpha = 0.5)
cvfitP12_09 <- cv.glmnet(DMsmallTest, ytest12,  family = "poisson", alpha = 0.9 ,foldid = iFoldTest)
cvfitP12_09 <- cv.glmnet(DMsmallTest, ytest12,  family = "poisson", alpha = 0.9 ,foldid = iFoldTest, type.measure = "mse")


# test bernoulli glm







