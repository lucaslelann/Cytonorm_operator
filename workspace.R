library(tercen)
library(dplyr)
library(devtools)

library(flowCore)
library(FlowSOM)
install_github('saeyslab/CytoNorm')


options("tercen.workflowId" = "c760adb979f713659d9e4826b004a670")
options("tercen.stepId"     = "53508f67-4c40-4bfb-9c5d-c50a28a6b129")

getOption("tercen.workflowId")
getOption("tercen.stepId")
ctx <- tercenCtx()
(ctx = tercenCtx()) %>% 
  select(.ga.gs0.variable)
(ctx = tercenCtx())  %>% 
  select(.y, .ci, .ri) %>% 
  group_by(.ci, .ri) %>%
  summarise(median = median(.y)) %>%
  ctx$addNamespace() %>%
  ctx$save()

dir <- "/home/rstudio/projects/CytoNorm/Cytonorm_operator"
files <- list.files(dir, pattern = "fcs$")
data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Type = stringr::str_match(files, "_([12]).fcs")[, 2],
                   Batch = stringr::str_match(files, "PTLG[0-9]*")[, 1],
                   stringsAsFactors = FALSE)
data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
# 
# train_data <- dplyr::filter(data, Type == "Train")
# validation_data <- dplyr::filter(data, Type == "Validation")
# 
# ff <- flowCore::read.FCS(data$Path[1])
# channels <- flowCore::colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
#                                      40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
#                                      49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30, 
#                                      28, 29, 25, 35)]
# transformList <- flowCore::transformList(channels,
#                                          cytofTransform)
# transformList.reverse <- flowCore::transformList(channels,
#                                                  cytofTransform.reverse)
# fsom <- prepareFlowSOM(train_data$Path,
#                        channels,
#                        nCells = 6000,
#                        FlowSOM.params = list(xdim = 5,
#                                              ydim = 5,
#                                              nClus = 10,
#                                              scale = FALSE),
#                        transformList = transformList,
#                        seed = 1)
# 
# cvs <- testCV(fsom,
#               cluster_values = c(5, 10, 15)) 
# 
# cvs$pctgs$`10`
