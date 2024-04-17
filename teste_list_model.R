
eval_myBiomodModelOut3["ROC","sensitivity","RF",,]


meanSenRF3 <- mean(eval_myBiomodModelOut3["ROC","Sensitivity","RF",,])
desvSenRF3 <- sd(eval_myBiomodModelOut3["ROC","Sensitivity","RF",,])

eval_myBiomodModelOut3["ROC","Specificity","RF",,]
meanSpeRF3 <- mean(eval_myBiomodModelOut3["ROC","Specificity","RF",,])
desvSpeRF3 <- sd(eval_myBiomodModelOut3["ROC","Specificity","RF",,])


str(myBiomodModelOut3)

myBiomodModelOut3@
