#install_github('saeyslab/CytoNorm')

library(tercen)
library(dplyr)
library(flowCore)
library(FlowSOM)
library(devtools)
library(CytoNorm)
library(shiny)

#docID
#options("tercen.workflowId" = "c760adb979f713659d9e4826b004a670")
#options("tercen.stepId"     = "9c30d12f-7040-4ca1-a4a9-14ba31986a38")

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
#data <- ctx$select()

data_all <-as.matrix(ctx) %>% t()
colnames(data_all) <- ctx$rselect()[[1]]
data_all <-cbind(data_all, ctx$cselect())

chan_nb <- length(ctx$rselect()[[1]])

train_data <- data_all[data_all["js0.Type"]== "Train",]
validate_data <- data_all[data_all["js0.Type"]== "control",]

#create temporary file 

dir.create("train")
for (js0.filename in unique(train_data$"js0.filename"))     {
  tmp_file_data <- train_data[train_data["js0.filename"] == js0.filename,]
  flow.dat <- flowCore::flowFrame(as.matrix(tmp_file_data[1:chan_nb]))
  
  outfile<-paste("train/",js0.filename, sep="")
  write.FCS(flow.dat, outfile)
}

dir.create("validate")
for (js0.filename in unique(validate_data$"js0.filename"))     {
  tmp_file_data <- train_data[validate_data["js0.filename"] == js0.filename,]
  flow.dat <- flowCore::flowFrame(as.matrix(tmp_file_data[1:chan_nb]))
  
  outfile<-paste("validate/",js0.filename, sep="")
  write.FCS(flow.dat, outfile)
}

# Open 1 of the flow cytometry files to get the channels and set uo the tranformation

channels <- flowCore::colnames(train_data)[c(1:chan_nb)]
transformList <- flowCore::transformList(channels,
                                         cytofTransform)
transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)

############################ Preparation 
list_train<-list.files("train",full.names = TRUE)
list_validate<-list.files("validate",full.names = TRUE)

fsom <- prepareFlowSOM(list_train,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

########################### Model determination
model <- CytoNorm.train(files = list_train,
                        labels = list_train,
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
########################## Application of CytoNorm

CytoNorm.normalize(model = model,
                   files = list_validate,
                   labels = list_validate,
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
      actual = get("actual",  envir = .GlobalEnv) + 1
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


# library(tercen)
# library(dplyr)
# library(flowCore)
# library(FlowSOM)
# library(devtools)
# library(CytoNorm)
# 
# 
# ### use for worspace
# #options("tercen.workflowId" = "c760adb979f713659d9e4826b004a670")
# #options("tercen.stepId"     = "5deafcbd-8d1a-4a40-b14c-7e8b8b4795b6")
# #getOption("tercen.workflowId")
# #getOption("tercen.stepId")
# 
# ############################### FUNCTION
# # fcs_to_data
# # input filename of fcs data
# # output dataframe of the fcs data
# fcs_to_data = function(filename, which.lines) {
#   data_fcs = read.FCS(filename, which.lines, transformation = FALSE)
#   names_parameters = data_fcs@parameters@data$desc
#   data = as.data.frame(exprs(data_fcs))
#   col_names = colnames(data)
#   names_parameters = ifelse(is.na(names_parameters),col_names,names_parameters)
#   colnames(data) = names_parameters
#   data %>%
#     mutate_if(is.logical, as.character) %>%
#     mutate_if(is.integer, as.double) %>%
#     mutate(.ci = rep_len(0, nrow(.))) %>%
#     mutate(filename = rep_len(basename(filename), nrow(.)))
# }
# 
# 
# 
# ############################## read FCS files
# 
# # get the input from tercen
# ctx <- tercenCtx()
# df <- ctx$cselect()
# 
# which.lines <- NULL
# if(!is.null(ctx$op.value('which.lines')) && !ctx$op.value('which.lines') == "NULL") which.lines <- as.integer(ctx$op.value('which.lines'))
# f.names<-NULL
# docId = df$J.Documentid[1]
# 
# #create temporary file 
# 
# for (docId in  df$J.Documentid){
#   doc = ctx$client$fileService$get(substr(docId,1,nchar(docId)-1))
#   filename = tempfile(fileext = ".fcs")
#   writeBin(ctx$client$fileService$download(substr(docId,1,nchar(docId)-1)), filename)
#   on.exit(unlink(filename))
#   
#   # unzip if archive
#   if(length(grep(".zip", doc$name)) > 0) {
#     tmpdir <- tempfile( fileext = ".fcs")
#     unzip(filename, exdir = tmpdir)
#     f.names <- list.files(tmpdir, full.names = TRUE)
#   } else {
#     f.names <- c(f.names,filename)
#   }
#   
#   # check FCS
#   if(any(!isFCSfile(f.names))) stop("Not all imported files are FCS files.")
#   
#   assign("actual", 0, envir = .GlobalEnv)
#   task = ctx$task
#   
# }
# ############################## CytoNorm
# 
# data.input <- data.frame(File = f.names,
#                          #Path = file.path(dir, files),
#                          Type = df$J.type,
#                          Batch = df$J.filename,
#                          stringsAsFactors = FALSE)
# 
# 
# # separation of the data in training and validation dataset
# train_data <- dplyr::filter(data.input , Type == "train")
# validation_data <- dplyr::filter(data.input , Type == "validate")
# 
# 
# # Open 1 of the flow cytometry files to get the channels and set uo the tranformation
# ff <- flowCore::read.FCS(data.input $File[1])
# channels <- flowCore::colnames(ff)[c(3:55)]
# transformList <- flowCore::transformList(channels,
#                                          cytofTransform)
# transformList.reverse <- flowCore::transformList(channels,
#                                                  cytofTransform.reverse)
# ############################ Preparation 
# fsom <- prepareFlowSOM(train_data$File,
#                        channels,
#                        nCells = 6000,
#                        FlowSOM.params = list(xdim = 5,
#                                              ydim = 5,
#                                              nClus = 10,
#                                              scale = FALSE),
#                        transformList = transformList,
#                        seed = 1)
# 
# #cvs <- testCV(fsom, cluster_values = c(5, 10, 15)) 
# #cvs$pctgs$`10`
# ########################### Model determination
# model <- CytoNorm.train(files = train_data$File,
#                         labels = train_data$Batch,
#                         channels = channels,
#                         transformList = transformList,
#                         FlowSOM.params = list(nCells = 6000, 
#                                               xdim = 5,
#                                               ydim = 5,
#                                               nClus = 10,
#                                               scale = FALSE),
#                         normMethod.train = QuantileNorm.train,
#                         normParams = list(nQ = 101,
#                                           goal = "mean"),
#                         seed = 1,
#                         verbose = TRUE)
# 
# 
# ########################## Application of CytoNorm
# 
# CytoNorm.normalize(model = model,
#                    files = validation_data$File,
#                    labels = validation_data$Batch,
#                    transformList = transformList,
#                    transformList.reverse = transformList.reverse,
#                    normMethod.normalize = QuantileNorm.normalize,
#                    outputDir = "Normalized",
#                    prefix = "Norm_",
#                    clean = TRUE,
#                    verbose = TRUE)
# 
# ############################# Output
# 
# f.names<- paste("./Normalized/",list.files(path="./Normalized", pattern="Norm_"), sep="")
# test.fun<-f.names%>%
#   lapply(function(filename){
#     data = fcs_to_data(filename, which.lines)
#     if (!is.null(task)) {
#       # task is null when run from RStudio
#       actual = get("actual",  envir = .GlobalEnv) + 1
#       assign("actual", actual, envir = .GlobalEnv)
#       evt = TaskProgressEvent$new()
#       evt$taskId = task$id
#       evt$total = length(f.names)
#       evt$actual = actual
#       evt$message = paste0('processing FCS file ' , filename)
#       ctx$client$eventService$sendChannel(task$channelId, evt)
#     } else {
#       cat('processing FCS file ' , filename)
#     }
#     data
#   }) %>%
#   bind_rows() 
# 
# # replace the documenId for the file names 
# for (i in c(1:length(data.input$File))){
#   test.fun$filename[grepl(strsplit(data.input$File,"/")[[i]][4], test.fun$filename)] <- data.input$Batch[[i]]
# } 
# 
# 
# test.fun%>%
#   
#   ctx$addNamespace()  %>%
#   ctx$save()
# 
# unlink("Normalized",recursive = TRUE)
