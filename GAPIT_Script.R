#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
biocLite("multtest")
install.packages("gplots")
install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d
install.packages("Matrix", dependencies = TRUE)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")
source("gapit_functions.txt")

##########################Load files#################################################################
myY  <- read.csv(file.choose(), head = TRUE)
str(myY)
myG <- read.table(file.choose(), head = FALSE)


##############L########### Statistical distributions of phenotype#######
myPhenotype<-GAPIT.Phenotype.View(
  myY=myY
)
########Models###########
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3,
  model="FarmCPU"
)


