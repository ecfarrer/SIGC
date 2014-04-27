##Split data
sampling.function<-function(species.info,trainProp){
  train<-sample(1:dim(species.info)[1],floor(dim(species.info)[1]*trainProp))
  test<-setdiff(1:dim(species.info)[1],  train)
  return(list(train=species.info[train,],test=species.info[test,]))
}


##Create list with training and evaluations matrices for our interaction comparisons
make.training.eval<-function(AB,Ao,Bo,oo,trainProp){
  
  #remove sites that have NAs for Worldclim bioclim values
  AB<-AB[is.na(AB[,1])==FALSE,]
  Ao<-Ao[is.na(Ao[,1])==FALSE,]
  Bo<-Bo[is.na(Bo[,1])==FALSE,]
  oo<-oo[is.na(oo[,1])==FALSE,]
  
  PseudoA <- rbind(Bo,oo) 
  PseudoB <- rbind(Ao,oo)
  
  #sample pseudoabsence points from p00
  set.seed(0)
  aaa<-sort(sample(nrow(PseudoA), nrow(rbind(AB,Ao))))
  bbb<-sort(sample(nrow(PseudoB), nrow(rbind(AB,Bo))))
  
  backgrA<-PseudoA[aaa,]#backgrA <- randomPoints(predictors.crop, dim(spApres)[1])
  backgrB<-PseudoB[bbb,]#backgrB <- randomPoints(predictors.crop, dim(spBpres)[1])
    
  #split the data to the training and evaluation sets
  one<-sampling.function(AB,trainProp)
  two<-sampling.function(Ao,trainProp)
  thr<-sampling.function(Bo,trainProp)
  fur<-sampling.function(backgrA,trainProp)
  fiv<-sampling.function(backgrB,trainProp)
  
  #determine which of the pseudoabsence are points of occurrence of other species
  sixT<-which(fur[[1]][,22] %in% Bo[,22] && fiv[[1]][,21] %in% Bo[,21])
  sixE<-which(fur[[2]][,22] %in% Bo[,22] && fur[[2]][,21] %in% Bo[,21])
  sevT<-which(fiv[[1]][,22] %in% Ao[,22] && fiv[[1]][,21] %in% Ao[,21])
  sevE<-which(fiv[[2]][,22] %in% Ao[,22] && fiv[[2]][,21] %in% Ao[,21])
  
  ##set up Presence/Absence column for each species for each set
  #AB, A0 == 1; 00 == 0
  APAtrain<-c(rep(1,c(dim(one[[1]])[1]+dim(two[[1]])[1])),
                rep(0,dim(fur[[1]])[1]))
  APAeval<-c(rep(1,c(dim(one[[2]])[1]+dim(two[[2]])[1])),
               rep(0,dim(fur[[2]])[1]))
  #AB, 0B == 1; 00 == 0
  BPAtrain<-c(rep(1,c(dim(one[[1]])[1]+dim(thr[[1]])[1])),
                rep(0,dim(fiv[[1]])[1]))
  BPAeval<-c(rep(1,c(dim(one[[2]])[1]+dim(thr[[2]])[1])),
#                 rep(0,c(dim(two[[2]])[1]+dim(fur[[2]])[1])))
               rep(0,dim(fiv[[2]])[1]))
  
  #Set up Other species column
  APAtrainB<-c(rep(1,dim(one[[1]])[1]),
               rep(0,dim(two[[1]])[1]),
               rep(1,length(sixT)),
               rep(0,c(dim(fur[[1]])[1]-length(sixT))))
  APAevalB<-c(rep(1,dim(one[[2]])[1]),
              rep(0,dim(two[[2]])[1]),
              rep(1,length(sixE)),
              rep(0,c(dim(fur[[2]])[1]-length(sixE))))
  #AB, 0B == 1; 00 == 0
  BPAtrainA<-c(rep(1,dim(one[[1]])[1]),
               rep(0,dim(thr[[1]])[1]),
               rep(1,length(sevT)),
               rep(0,c(dim(fiv[[1]])[1]-length(sevT))))
  BPAevalA<-c(rep(1,dim(one[[2]])[1]),
              rep(0,dim(thr[[2]])[1]),
              rep(1,length(sevE)),
              rep(0,c(dim(fiv[[2]])[1]-length(sevE))))
  
  #create training and evaluation env. variable set for species A
  trainA<-rbind(one[[1]],two[[1]],fur[[1]])
  evalA<-rbind(one[[2]],two[[2]],fur[[2]])
  
  #standardize (mean =0, var=1) based upon the values from the training set
  mean.trainA<-apply(trainA,2,mean,na.rm=TRUE)
  sd.trainA<-apply(trainA,2,sd,na.rm=TRUE)
  
  trainA <- trainA - kronecker(matrix(rep(1,dim(trainA)[1]), ncol=1), 
                              matrix(mean.trainA,nrow=1));
  trainA <- trainA / kronecker(matrix(rep(1,dim(trainA)[1]), ncol=1), 
                               matrix(sd.trainA,nrow=1));
  evalA <- evalA - kronecker(matrix(rep(1,dim(evalA)[1]), ncol=1), 
                             matrix(mean.trainA,nrow=1));
  evalA <- evalA / kronecker(matrix(rep(1,dim(evalA)[1]), ncol=1), 
                             matrix(sd.trainA,nrow=1));
  
  #create training and evaluation env. variable set for species B
  trainB<-rbind(one[[1]],thr[[1]],fiv[[1]])
  evalB<-rbind(one[[2]],thr[[2]],fiv[[2]])
  
  #standardize (mean =0, var=1) based upon the values from the training set
  mean.trainB<-apply(trainB,2,mean,na.rm=TRUE)
  sd.trainB<-apply(trainB,2,sd,na.rm=TRUE)
  
  trainB <- trainB - kronecker(matrix(rep(1,dim(trainB)[1]), ncol=1), 
                               matrix(mean.trainB,nrow=1));
  trainB <- trainB / kronecker(matrix(rep(1,dim(trainB)[1]), ncol=1), 
                               matrix(sd.trainB,nrow=1));
  evalB <- evalB - kronecker(matrix(rep(1,dim(evalB)[1]), ncol=1), 
                             matrix(mean.trainB,nrow=1));
  evalB <- evalB / kronecker(matrix(rep(1,dim(evalB)[1]), ncol=1), 
                             matrix(sd.trainB,nrow=1));
  
  #bind the pres/abs data to the environmental data
  trainA<-cbind(trainA,APAtrain,APAtrainB)
  trainB<-cbind(trainB,BPAtrain,BPAtrainA)
  
  evalA<-cbind(evalA,APAeval,APAevalB)
  evalB<-cbind(evalB,BPAeval,BPAevalA)
  
  return(list(trainA=trainA,
              evalA=evalA,
              trainB=trainB,
              evalB=evalB))
  
}