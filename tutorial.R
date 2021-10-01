library(tercen)
library(dplyr)
library(flowCore)
library(FlowSOM)
library(devtools)
library(CytoNorm)

#Script taken from the cytonorm git page 
#Serve as a tutorial


dir <- system.file("extdata", package = "CytoNorm")
files <- list.files(dir, pattern = "fcs$")
data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Type = stringr::str_match(files, "_([12]).fcs")[, 2],
                   Batch = stringr::str_match(files, "PTLG[0-9]*")[, 1],
                   stringsAsFactors = FALSE)

data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]

train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")

ff <- flowCore::read.FCS(data$Path[1])
channels <- flowCore::colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
                                     40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
                                     49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30, 
                                     28, 29, 25, 35)]
transformList <- flowCore::transformList(channels,
                                         cytofTransform)
transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)

fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

cvs <- testCV(fsom,
              cluster_values = c(5, 10, 15)) 

cvs$pctgs$`10`

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 6000, 
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 10,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)

CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized_tutorial",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)

df_norm<-paste("./Normalized_tutorial/",list.files(path="./Normalized_tutorial", pattern="Norm_"), sep="")
flow_norm <- flowCore::read.FCS(df_norm[1])

