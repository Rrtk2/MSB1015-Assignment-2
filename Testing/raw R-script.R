#-----------------------------------------------------------------------------------------------------#
#		Block 00		GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#

# Copyright statement comment:
#   All rights reserved.
#
# Author comment:
#   Rick Reijnders
#   Script version: 29-09-2019

# File description:
#	Name
#	  raw R-script.R
#
#   Purpose
#     
#   Inputs
#	  

#-----------------------------------------------------------------------------------------------------#
#		Block 01		(Install &) Load packages
#-----------------------------------------------------------------------------------------------------#
# Install packages if needed and load, or just load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

# Insert all packages in requiredpackages
requiredpackages <-
  c("WikidataQueryServiceR","ggplot2","rJava","rcdk","pls","e1071","neuralnet","randomForest",
	"crayon")
	
for (i in requiredpackages) {
	if (!requireNamespace(i, quietly = TRUE))
		BiocManager::install(i, ask = F, dependencies = TRUE)
	require(as.character(i), character.only = TRUE)
	print(i)
}

# Set message styles
Confirm.style <- green 
Warning.style <- yellow 
Error.style <- red


cat(Warning.style("Done running block 01\n"))
cat(Confirm.style("Starting...\n"))

#-----------------------------------------------------------------------------------------------------#
#		Block 02		SETTINGS
#-----------------------------------------------------------------------------------------------------#
# General settings
options(stringsAsFactors 	= F)
Verbose 					= 3		#0-3:	0= no feedback  1= prints results  2=prints results + feedback  3= prints all

# Select 80% for train, 20% for test
Randomfactor = 0.8

# Setting seed to keep consistent results 
set.seed(1)

if(Verbose>=3) cat(Warning.style("Done running block 02a\n"))

#-----------------------------------------------------------------------------------------------------#
#		Block 03		FUNCTIONS
#-----------------------------------------------------------------------------------------------------#

# Root Mean Squared Error
RMSE = function(yact, ypred){
  sqrt(mean((yact - ypred)^2))
}

# Mean Absolute Error 
MAE = function(yact, ypred){
  mean(abs(yact - ypred))
}


#-----------------------------------------------------------------------------------------------------#
#		Block 03		Get query results
#-----------------------------------------------------------------------------------------------------#
endpoint = "https://query.wikidata.org/bigdata/namespace/wdq/sparql"

query = 'SELECT ?comp ?compLabel ?bp ?bpUnitLabel ?CC WHERE {
  ?comp wdt:P31/wdt:P279* wd:Q41581 ;
        p:P2102 [
          ps:P2102 ?bp ;
          psv:P2102/wikibase:quantityUnit  ?bpUnit
        ] .
		?comp wdt:P233 ?CC .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
'

dataObj = query_wikidata(query)

#-----------------------------------------------------------------------------------------------------#
#		Block 03		bp check ; to kelvin
#-----------------------------------------------------------------------------------------------------#
#C to Kelvin:
# 0°C + 273.15 = 273,15K
idC = which(dataObj$bpUnitLabel=="degree Celsius")
for( i in 1:length(idC)){
dataObj[idC[i],]$bp = (dataObj[idC[i],]$bp + 273.15)
dataObj[idC[i],]$bpUnitLabel = "kelvin"
}

#F to Kelvin:
# (0°F − 32) × 5/9 + 273.15 = 255,372K
idF = which(dataObj$bpUnitLabel=="degree Fahrenheit")
for( i in 1:length(idF)){
dataObj[idF[i],]$bp= (dataObj[idF[i],]$bp - 32) * (5/9) + 273.15
dataObj[idF[i],]$bpUnitLabel = "kelvin"
}

#additional failsafe if other metrics than Celsius and Fahrenheit are added (other metrics are removed)
dataObj = dataObj[which(dataObj$bpUnitLabel=="kelvin"),]

# Make backup
dataObjBACKUP = dataObj

# Finalize data structure
dataObj = data.frame(dataObjBACKUP$compLabel, dataObjBACKUP$bp, dataObjBACKUP$CC)
names(dataObj) = c("Comp","bp","CC")
#-----------------------------------------------------------------------------------------------------#
# 		Block 03		rcdk data extraction see:https://cran.r-project.org/web/packages/rcdk/vignettes/molform.html
#-----------------------------------------------------------------------------------------------------#

sp <- get.smiles.parser()
dc <- get.desc.categories()

for( i in 1:length(dataObj$CC)){

	# Get smiles from result query and add information
	molecule <- parse.smiles(dataObj$CC[i])[[1]]
	convert.implicit.to.explicit(molecule)
	#formula <- get.mol2formula(molecule,charge=0)


	# Store the found information
	#dataObj$formula[i] = {formula} # S4 object cannot be transferred nicely
	#dataObj$mass[i] = formula@mass
	#dataObj$string[i] = formula@string
	#dataObj$charge[i] = formula@charge

	# M/z values 
	#dataObj$isotopes[i] = {get.isotopes.pattern(formula,minAbund=0.1)}

	# Fingerprint? values
	#dataObj$fingerprint[i] = {get.fingerprint(molecule = molecule)}

	# Create a dataframe which contains all info possible to extract using the descriptors
	datafr = dataObj$Comp[i]
	for(o in 1:5){
		dn <- get.desc.names(dc[o])
		datafr = cbind(datafr, eval.desc(molecule, dn))
	}
	
	# This is done now as it will break if descriptors change (amount)
	if(exists("mydata")){
		mydata[i,] = datafr
		}else{
		mydata = datafr
	}

}

mydataBACKUP = mydata
mydata = mydataBACKUP
descs = mydata[,-1]

if(T){
	# crude way to extract 'important' features based on correlation
	descs <- descs[, !apply(descs, 2, function(x) any(is.na(x)) )]
	descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )]
	r2 <- which(cor(descs)^2 > .9, arr.ind=TRUE)
	r2 <- r2[ r2[,1] > r2[,2] , ]
	descs <- descs[, -unique(r2[,2])]
}

my_data = descs
#-----------------------------------------------------------------------------
#Finish prepare data
#-----------------------------------------------------------------------------
#my_data = my_data[,-1]  # Remove ID
n_col <- ncol(my_data)
BoilPoint = n_col+1
my_data[,BoilPoint] = dataObj$bp
names(my_data)[BoilPoint] = 'BoilPoint'



#-----------------------------------------------------------------------------------------------------#
# 				Data subsetting; Train; Test
#-----------------------------------------------------------------------------------------------------#
# ceiling(dim(my_data)[1]0.8)
#make data objects
samples.upper = sample(length(my_data[,1]), floor(length(my_data[,1])*Randomfactor))  #get all unique samples (not frequecies) and sample 80%
samples.total = (1:length(my_data[,1])) # all unique samples (not frequecies)
samples.lowerl = samples.total[!samples.total %in% samples.upper]  #which samples are sampled
data.train = my_data[samples.upper,]#contains all upper logical
data.train = as.data.frame(data.train)

#data.train = t(data.train)
data.test = my_data[samples.lowerl,]#contains all lower logical
data.test = as.data.frame(data.test)
#data.test = t(data.test)

# NOTE: annotations are within location [,32007]
data.train[,BoilPoint] = data.train[,BoilPoint]
data.test[,BoilPoint] =data.test[,BoilPoint]

#remove nas
data.train=data.train[!is.na(data.train[,1]),]
data.test=data.test[!is.na(data.test[,1]),]
set.seed(123)


xNN = data.train[,-(BoilPoint)]
yNN = data.train[,BoilPoint]

dat = data.frame(xNN, y = yNN)
yactual = data.test[,BoilPoint]



#-----------------------------------------------------------------------------------------------------#
#							Linear Regression
#-----------------------------------------------------------------------------------------------------#

# S3 method for formula
GLMfit = glm(y ~ ., dat,family = "gaussian")

#so put data.predicted and actual data.test together
ypredGLM = predict(GLMfit, as.data.frame(data.test[,-(BoilPoint)]),type="response")

# Root mean squared error
RMSE(yactual, ypredGLM)
#	Results in: 36.76

# Mean absolute error
MAE(yactual, ypredGLM)
#	Results in: 22.11

plot(dat$y, predict(GLMfit, xNN),
     xlab="Observed BP", ylab="Predicted BP",
     pch=19, xlim=c(100, 700), ylim=c(100, 700))
abline(0,1, col='red')

# check for under/overfitting
 ypredGLM = predict(GLMfit, as.data.frame(data.train[,-(BoilPoint)]),type="response")
# Root mean squared error
RMSE(data.train[,(BoilPoint)], ypredGLM)
#	Results in: 23.32

# Mean absolute error
MAE(data.train[,(BoilPoint)], ypredGLM)
#	Results in: 7.844

#-----------------------------------------------------------------------------------------------------#
#							Linear Regression
#-----------------------------------------------------------------------------------------------------#


# S3 method for formula
PLSfit = plsr(y ~ ., ncomp = 20, data = dat, validation = "LOO")

# Select amount of components
ncomp.onesigma <- selectNcomp(PLSfit, method = "onesigma", plot = TRUE)
ncomp.permut <- selectNcomp(PLSfit, method = "randomization", plot = TRUE)

# If methods doubt between 5 or 7 components, take the highest rounded mean (6).
PLSfit = plsr(y ~ ., ncomp = ceiling(mean(ncomp.onesigma,ncomp.permut)), data = dat, validation = "LOO")

#so put data.predicted and actual data.test together
ypredPLS = predict(PLSfit, as.data.frame(data.test[,-(BoilPoint)]),type="response")

# Root mean squared error
RMSE(yactual, ypredPLS)
#	Results in: 53.59

# Mean absolute error
MAE(yactual, ypredPLS)
#	Results in: 39.25

#-----------------------------------------------------------------------------------------------------#
# 										RF
#-----------------------------------------------------------------------------------------------------#


prot.rf <- randomForest(y ~ ., importance=TRUE, proximity=TRUE,data=dat, ntree=200)

# after building the random forest, now we apply it on the test dataset
ypredRF <- predict(prot.rf, data.test[,-(BoilPoint)])


# Root mean squared error
RMSE(yactual, ypredRF)
#	Results in: 20.26

# Mean absolute error
MAE(yactual, ypredRF)
#	Results in: 7.63




