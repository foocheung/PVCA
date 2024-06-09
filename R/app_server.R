#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#'
pvca <- function (abatch,expInfo) {
  theDataMatrix=abatch
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)

  ########## Center the data (center rows) ##########
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)

  ########## Compute correlation matrix &  Obtain eigenvalues ##########

  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues /eigenValuesSum

  ##===========================================
  ##  Getting the experimental information
  ##===========================================

  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)

  ########## Merge experimental file and eigenvectors for n components ##########

  ## pc_n is the number of principal components to model

  ## Use fixed pc_n
  pc_n=5

  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
  mycounter = 0 
  for (i in 1:pc_n){
    for (j in 1:expDesignRowN){
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
    }
  }

  AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
  Data <- cbind(AAA,pc_data_matrix)


  ####### Edit these variables according to your factors #######

  variables <-c (colnames(exp_design))
  for (i in 1:length(variables))
  {
    Data$variables[i] <- as.factor(Data$variables[i] )
  }


  ########## Mixed linear model ##########
  op <- options(warn = (-1))
  effects_n = expDesignColN+1 #effects size without interaction terms
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

  ##============================#
  ##  Get model functions
  ##============================#
  model.func <- c()
  index <- 1

  ##  level-1
  for (i in 1:length(variables))
  {
    mod = paste("(1|", variables[i], ")",   sep="")
    model.func[index] = mod
    index = index + 1
  }

  function.mods <- paste (model.func , collapse = " + ")

  ##============================#
  ##  Get random effects  #
  ##============================#

  for (i in 1:pc_n){
    y = (((i-1)*expDesignRowN)+1)
    funct <- paste ("pc_data_matrix", function.mods, sep =" ~ ")
    Rm1ML <- lme4::lmer( funct ,
                   Data[y:(((i-1)*expDesignRowN)+expDesignRowN),],
                   REML = TRUE, verbose = FALSE, na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i,] <- c(unlist(lme4::VarCorr(Rm1ML)),resid=sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML,"cnms")),"resid")

  ########## Standardize Variance ##########
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    mySum = sum(randomEffectsMatrix[i,])
    for (j in 1:effects_n){
      randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum
    }
  }

  ########## Compute Weighted Proportions ##########

  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n){
      randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
    }
  }

  ########## Compute Weighted Ave Proportions ##########

  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

  for (j in 1:effects_n){
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat=randomEffectsMatrixWtAveProp, label=effectsNames))
}

app_server <- function(input, output, session) {
  # Your application server logic

  observeEvent(input$file_exprs, {
    withProgress(message = "Generating PCA plot...", value = 0, {
      exprs <- readr::read_tsv(input$file_exprs$datapath)
      aaa<<-exprs
      row.names(exprs)<-exprs[[1]]
      exprs<-exprs[-1]
    })
  })

  observeEvent(input$file_metadata, {
    withProgress(message = "Generating PCA plot...", value = 0, {
      metadata <- readr::read_tsv(input$file_metadata$datapath)
      bbb<<-metadata
      row.names(metadata)<-metadata[[1]]
      metadata <- metadata[-1]
    })
  })

  # Generate PCA plot
  output$pca_plot <- renderPlot({
    req(input$submit)
    withProgress(message = "Generating PVCA plot...", value = 0, {
      metadata <- readr::read_tsv(input$file_metadata$datapath)
      exprs <- readr::read_tsv(input$file_exprs$datapath)
      aa<<-exprs
      row.names(exprs)<-exprs[[1]]
      exprs<-exprs[-1]

      bb<<-metadata
      row.names(metadata)<-metadata[[1]]
      metadata <- metadata[-1]

      e<<-exprs
      m<<-metadata
      pvcaObj <- pvca(exprs, metadata)
      par(mar = c(9, 6, 4, 1))
      p <- barplot(pvcaObj$dat, ylab = "Proportion variance explained",
                   ylim = c(0, 1.1), col = "blue", las = 2, cex.axis = 1, cex.names = 1)
      axis(1, at = p, labels = pvcaObj$label, cex.axis = 1, las = 2)
      p
    })
  })




  }
