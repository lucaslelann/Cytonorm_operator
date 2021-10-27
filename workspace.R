#install_github('saeyslab/CytoNorm')

library(tercen)
library(dplyr)
library(flowCore)
library(FlowSOM)
library(devtools)
library(CytoNorm)
#system.file("extdata", package = "CytoNorm")

#docID
options("tercen.workflowId" = "835f1113e61a613dedcbbaf7640313ed")
options("tercen.stepId"     = "108ed769-3b37-4134-b4a7-aa48d6afb196")

getOption("tercen.workflowId")
getOption("tercen.stepId")

############################### FUNCTION
# fcs_to_data
# input filename of fcs data
# output dataframe of the fcs data
fcs_to_data = function(filename) {
  data_fcs = read.FCS(filename, transformation = FALSE)
  names_parameters = data_fcs@parameters@data$desc
  data = as.data.frame(exprs(data_fcs))
  col_names = colnames(data)
  names_parameters = ifelse(is.na(names_parameters),col_names,names_parameters)
  colnames(data) = names_parameters
  data %>%
    mutate_if(is.logical, as.character) %>%
    mutate_if(is.integer, as.double) %>%
    mutate(.ci = rep_len(0, nrow(.))) %>%
    mutate(filename = rep_len(basename(filename), nrow(.)))
}

############################## read FCS files

# get the input from tercen
ctx <- tercenCtx()
task<-ctx$task
nclust <- as.double(ctx$op.value('cluster'))
ncells <- as.double(ctx$op.value('number_of_cells'))

#option set workflow default option set in tutorial
nclust <- 10
ncells <- 6000

data_all <-as.matrix(ctx) %>% t()
colnames(data_all) <- ctx$rselect()[[1]]
data_all <-cbind(data_all, ctx$cselect())

chan_nb <- length(ctx$rselect()[[1]])

colnames(data_all)[grep("[T,t]ype",colnames(data_all))]<-"type"
colnames(data_all)

colnames(data_all)[grep("[F,f]ilename",colnames(data_all))]<-"filename"
colnames(data_all)

colnames(data_all)[grep("[B,b]atch",colnames(data_all))]<-"batch"
colnames(data_all)

train_data <- data_all[data_all["type"]== "Train",]
validate_data <- data_all[data_all["type"]== "validate",]
batch_train_data <- unique(data_all[data_all["type"]== "Train",]$batch)
batch_validate_data <- unique(data_all[data_all["type"]== "validate",]$batch)

#create temporary file 

dir.create("train")
for (filename in unique(train_data$"filename"))     {
  tmp_train_file_data <- train_data[train_data["filename"] == filename,]
  
  #write in the data file the channels without the annotation columns but with the time 
  flow.dat <- flowCore::flowFrame(as.matrix(tmp_train_file_data[c(1:(chan_nb),(chan_nb+5))]))
  outfile<-paste("train/",filename, sep="")
  write.FCS(flow.dat, outfile)
}

dir.create("validate")
for (filename in unique(validate_data$"filename"))     {
  tmp_val_file_data <- validate_data[validate_data["filename"] == filename,]
  #write in the data file the channels without the annotation columns but with the time 
  flow.dat <- flowCore::flowFrame(as.matrix(tmp_val_file_data[c(1:(chan_nb),(chan_nb+5))]))
  outfile<-paste("validate/",filename, sep="")
  write.FCS(flow.dat, outfile)
}

# Open 1 of the flow cytometry files to get the channels and set up the tranformation
#channels <- flowCore::colnames(train_data)[c(1:chan_nb)]

#channels <- colnames(train_data)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
                                     # 40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
                                     # 49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30, 
                                     # 28, 29, 25, 35)-3]
#channels <- colnames(train_data)[c(52,26,34,53,12,48,46,8,33,43,13,18,30,23,47,16,28,22,45,41,14,21,15,37,39,17,38,42,50,31,10)]
#channels <- colnames(validate_data)[c(52,26,34,53,12,48,46,8,33,43,13,18,30,23,47,16,28,22,45,41,14,21,15,37,39,17,38,42,50,31,10)]
#channels <- colnames(validate_data)[c(1:20)]
channels <- colnames(train_data)[c(0:chan_nb)]

transformList <- flowCore::transformList(channels,cytofTransform)
transformList.reverse <- flowCore::transformList(channels,cytofTransform.reverse)

############################ Preparation 
list_train<-list.files("train",full.names = TRUE)
list_validate<-list.files("validate",full.names = TRUE)

fsom <- prepareFlowSOM(list_train,
                       channels,
                       nCells = ncells,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = nclust,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

########################### Model determination
model <- CytoNorm.train(files = list_train,
                        labels =  batch_train_data,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = ncells, 
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = nclust,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)
########################## Application of CytoNorm

CytoNorm.normalize(model = model,
                   files = list_validate,
                   labels = batch_validate_data,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)

############################# Output

f.names<- paste("./Normalized/",list.files(path="./Normalized", pattern="Norm_"), sep="")
test.fun<-f.names%>%
  lapply(function(filename){
    data = fcs_to_data(filename)
    if (!is.null(task)) {
      # task is null when run from RStudio
      actual <-get("actual",  envir = .GlobalEnv) + 1
      assign("actual", actual, envir = .GlobalEnv)
      evt = TaskProgressEvent$new()
      evt$taskId = task$id
      evt$total = length(f.names)
      evt$actual = actual
      evt$message = paste0('processing FCS file ' , filename)
      ctx$client$eventService$sendChannel(task$channelId, evt)
    } else {
      cat('processing FCS file ' , filename)
    }
    data
  }) %>%
  bind_rows() 


test.fun%>%
  ctx$addNamespace()  %>%
  ctx$save()

unlink("train",recursive = TRUE)
unlink("validate",recursive = TRUE)
unlink("Normalized",recursive = TRUE)

