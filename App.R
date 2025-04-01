# Auteur : Axelle Bourez
# Date : 2025-03-26
# Licence : MIT (voir fichier LICENSE pour les détails)
# Description : Analyses Multiblocks Interactives

library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(reshape2)
library(ropls)
library(mixOmics)
library(MBAnalysis)
library(dplyr)
library(plotly)
library(ggrepel)
library(openxlsx)
library(ConsensusOPLS)
library(RVAideMemoire)
library(abind)
library(visNetwork)
library(igraph)


# cleanNumericData()
cleanNumericData <- function(df, idCol = 1) {
  df <- as.data.frame(df)
  rownames(df) <- as.character(df[[idCol]])
  df <- df[, -idCol, drop = FALSE]
  for (j in seq_along(df)) {
    df[[j]] <- as.numeric(as.character(df[[j]]))
  }
  
  # Retire les colonnes dont la variance est nulle
  var_ok <- sapply(df, function(x) var(x, na.rm = TRUE) > 0)
  df <- df[, var_ok, drop = FALSE]
  
  # Retire les lignes qui contiennent des NA
  df <- df[complete.cases(df), ]
  
  return(df)
}

# UI 

ui <- dashboardPage(
  dashboardHeader(title = "Analyses Multiblocks Interactives"),
  dashboardSidebar(
    fileInput("proteoFile", "Fichier protéomique (CSV)", accept = ".csv"),
    fileInput("metaboFile", "Fichier métabolomique (CSV)", accept = ".csv"),
    fileInput("metadataFile", "Fichier metadata (CSV)", accept = ".csv"),
    actionButton("goData", "Charger les données")
  ),
  dashboardBody(
    tabBox(
      width = 12,
      
      # Distribution
      tabPanel("Distribution",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goDistribution", "Lancer Distribution"),
                          br(), br(),
                          plotlyOutput("densiteNonScaled"),
                          br(),
                          plotlyOutput("densiteScaled"),
                          br(), br(),
                          downloadButton("download_Distribution", "Exporter Distribution")
                 )
               )
      ),
      
      # PCA Proteomique
      tabPanel("PCA Proteomique",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPCAProteo", "Lancer PCA Proteomique"),
                          br(), br(),
                          plotlyOutput("pcaProteoScree"),
                          br(),
                          plotlyOutput("pcaProteoScores12"),
                          br(),
                          plotlyOutput("pcaProteoScores34"),
                          br(),
                          plotlyOutput("pcaProteoLoadings"),
                          br(), br(),
                          downloadButton("download_PCAProteo", "Exporter PCA Proteomique")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("pcaProteoSummary")
                 )
               )
      ),
      
      # PCA Métabolomique
      tabPanel("PCA Métabolomique",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPCAMetabo", "Lancer PCA Métabolomique"),
                          br(), br(),
                          plotlyOutput("pcaMetaboScree"),
                          br(),
                          plotlyOutput("pcaMetaboScores12"),
                          br(),
                          plotlyOutput("pcaMetaboScores34"),
                          br(),
                          plotlyOutput("pcaMetaboLoadings"),
                          br(), br(),
                          downloadButton("download_PCAMetabo", "Exporter PCA Métabolomique")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("pcaMetaboSummary")
                 )
               )
      ),
      tabPanel("PLS‑DA Proteomique",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPLSDAProteo", "Lancer PLS‑DA Protéomique (mixOmics)"),
                          br(), br(),
                          plotlyOutput("plsdaProteoScores"),
                          br(),
                          plotlyOutput("plsdaProteoLoadingsBar1"),
                          br(),
                          plotlyOutput("plsdaProteoVIP"),
                          br(),
                          plotlyOutput("plsdaProteoLoadingsVIP"),
                          br(), br(),
                          downloadButton("download_PLSDAProteo", "Exporter PLS‑DA Protéomique"),
                          hr(),
                          actionButton("goOPLSProteo", "Lancer OPLS Protéomique (ropls)"),
                          br(), br(),
                          plotlyOutput("oplsProteoScores"),
                          br(),
                          plotlyOutput("oplsProteoLoadings")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("plsdaProteoSummary")
                 )
               )
      ),
      
      # PLS‑DA Métabolomique
      tabPanel("PLS‑DA Métabolomique",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPLSDAMetabo", "Lancer PLS‑DA Métabolomique (mixOmics)"),
                          br(), br(),
                          plotlyOutput("plsdaMetaboScores"),
                          br(),
                          plotlyOutput("plsdaMetaboLoadingsBar1"),
                          br(),
                          plotlyOutput("plsdaMetaboVIP"),
                          br(),
                          plotlyOutput("plsdaMetaboLoadingsVIP"),
                          br(), br(),
                          downloadButton("download_PLSDAMetabo", "Exporter PLS‑DA Métabolomique"),
                          hr(),
                          actionButton("goOPLSMetabo", "Lancer OPLS Métabolomique (ropls)"),
                          br(), br(),
                          plotlyOutput("oplsMetaboScores"),
                          br(),
                          plotlyOutput("oplsMetaboLoadings")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("plsdaMetaboSummary")
                 )
               )
      ),
      
      # PLS Regression
      
      tabPanel("PLS Regression",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPLSReg", "Lancer PLS Regression"),
                          br(), br(),
                          plotlyOutput("plsRegScoresProteo"),
                          br(),
                          plotlyOutput("plsRegScoresMetabo"),
                          br(),
                          plotlyOutput("plsRegLoadingsProteo"),
                          br(),
                          plotlyOutput("plsRegLoadingsMetabo"),
                          br(), br(),
                          downloadButton("download_PLSReg", "Exporter PLS Regression")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("plsRegSummary")
                 )
               )
      ),
      
      # PLS Canonical
      
      tabPanel("PLS Canonical",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goPLSCano", "Lancer PLS Canonical"),
                          br(), br(),
                          plotlyOutput("plsCanoScoresProteo"),
                          br(),
                          plotlyOutput("plsCanoScoresMetabo"),
                          br(),
                          plotlyOutput("plsCanoLoadingsProteo"),
                          br(),
                          plotlyOutput("plsCanoLoadingsMetabo"),
                          br(), br(),
                          downloadButton("download_PLSCano", "Exporter PLS Canonical")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("plsCanoSummary")
                 )
               )
      ),
      
      # ComDim
      
      tabPanel("ComDim",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          actionButton("goComDim", "Lancer ComDim"),
                          br(), br(),
                          plotlyOutput("comdimSaliences"),
                          br(),
                          plotlyOutput("comdimScores"),
                          br(),
                          plotlyOutput("comdimLoadings"),
                          br(), br(),
                          downloadButton("download_ComDim", "Exporter ComDim")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("comdimSummary")
                 )
               )
      ),
      
      # Block PLS‑DA
      
      tabPanel("Block-PLS‑DA",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          selectInput("groupVarBlock", "Variable de regroupement:", choices = NULL),
                          numericInput("ncompBlock", "Nombre de composantes:", value = 2, min = 1),
                          actionButton("goBlockPLSDA", "Lancer block PLS‑DA"),
                          br(), br(),
                          plotlyOutput("blockPLS_expl"),
                          br(),
                          plotlyOutput("blockPLS_scores"),
                          br(),
                          plotlyOutput("blockPLS_weighted"),
                          br(),
                          plotlyOutput("blockPLS_corrCircle"),
                          br(),
                          plotlyOutput("blockPLS_loadings"),
                          br(),
                          verbatimTextOutput("blockPLS_permtest"),
                          br(), br(),
                          downloadButton("download_BlockPLSDA", "Exporter block PLS‑DA")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("blockPLS_summary")
                 )
               )
      ),
      tabPanel("Réseau de corrélation",
               fluidRow(
                 box(width = 4,
                     sliderInput("networkCutoff", "Choisir le cutoff :", 
                                 min = 0, max = 1, value = 0.5, step = 0.01),
                     actionButton("goNetwork", "Générer le Network")
                 ),
                 box(width = 8,
                     visNetworkOutput("networkPlot"), 
                     downloadButton("download_Network", "Exporter")
                 )
               )
      ),
      
      # Consensus OPLS
      
      tabPanel("Consensus OPLS",
               tabsetPanel(
                 tabPanel("Graphiques",
                          br(),
                          selectInput("groupVarCons", "Variable de regroupement:", choices = NULL),
                          actionButton("goConsensusOPLS", "Lancer Consensus OPLS"),
                          br(), br(),
                          plotlyOutput("consOPLS_Q2"),
                          br(),
                          plotlyOutput("consOPLS_DQ2"),
                          br(),
                          plotlyOutput("consOPLS_R2Y"),
                          br(),
                          plotlyOutput("consOPLS_contrib"),
                          br(),
                          plotlyOutput("consOPLS_scores"),
                          br(),
                          plotlyOutput("consOPLS_loadings"),
                          br(),
                          plotlyOutput("consOPLS_VIP"),
                          br(),
                          plotlyOutput("consOPLS_VIP_Loadings"),
                          br(), br(),
                          downloadButton("download_ConsensusOPLS", "Exporter Consensus OPLS")
                 ),
                 tabPanel("Résumé",
                          verbatimTextOutput("consOPLS_summary")
                 )
               )
      )
    )
  )
)


# Serveur 

server <- function(input, output, session) {
  
  # Chargement des données 
  proteoData <- eventReactive(input$goData, {
    req(input$proteoFile)
    fread(input$proteoFile$datapath, sep = ";", header = TRUE)
  })
  
  metaboData <- eventReactive(input$goData, {
    req(input$metaboFile)
    fread(input$metaboFile$datapath, sep = ";", header = TRUE)
  })
  
  metadataData <- eventReactive(input$goData, {
    req(input$metadataFile)
    md <- fread(input$metadataFile$datapath, sep = ";", header = TRUE)
    
    if (!("sample_name" %in% colnames(md))) {
      md$sample_name <- as.character(md[[1]])
    }
    md$sample_name <- trimws(as.character(md$sample_name))
    md
  })
  
  # Mise à jour des variables de regroupement
  observeEvent(input$goData, {
    req(metadataData())
    md <- metadataData()
    groupChoices <- names(md)[sapply(md, function(x) is.character(x) || is.factor(x))]
    if (length(groupChoices) > 0) {
      updateSelectInput(session, "groupVarBlock", choices = groupChoices, selected = groupChoices[1])
      updateSelectInput(session, "groupVarCons", choices = groupChoices, selected = groupChoices[1])
    }
  })
  
  # Distribution
  
  densiteNonScaledData <- eventReactive(input$goDistribution, {
    req(proteoData())
    dt <- as.data.table(proteoData())
    melted <- melt(dt, id.vars = names(dt)[1],
                   variable.name = "variable",
                   value.name = "value")
    melted$value <- as.numeric(as.character(melted$value))
    melted
  })
  
  densiteScaledData <- eventReactive(input$goDistribution, {
    req(proteoData())
    dt <- as.data.frame(proteoData())
    proteo_numeric <- dt[, -1, drop = FALSE]
    for (j in seq_along(proteo_numeric)) {
      proteo_numeric[[j]] <- as.numeric(as.character(proteo_numeric[[j]]))
    }
    scaled_proteo <- scale(proteo_numeric)
    scaled_df <- as.data.frame(scaled_proteo)
    scaled_df$sample <- rownames(scaled_df)
    melted <- melt(scaled_df, id.vars = "sample",
                   variable.name = "variable",
                   value.name = "value")
    melted$value <- as.numeric(as.character(melted$value))
    melted
  })
  
  output$densiteNonScaled <- renderPlotly({
    req(densiteNonScaledData())
    p <- ggplot(densiteNonScaledData(), aes(x = value, color = variable)) +
      geom_density() +
      guides(color = "none") +
      labs(title = "Densité (non mises à l'échelle)", x = "Valeur", y = "Densité")
    ggplotly(p)
  })
  
  output$densiteScaled <- renderPlotly({
    req(densiteScaledData())
    p <- ggplot(densiteScaledData(), aes(x = value, color = variable)) +
      geom_density() +
      guides(color = "none") +
      labs(title = "Densité (mises à l'échelle)", x = "Valeur mise à l'échelle", y = "Densité")
    ggplotly(p)
  })
  
  # Download distribution
  output$download_Distribution <- downloadHandler(
    filename = function() { paste("Distribution_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      addWorksheet(wb, "NonScaled")
      writeData(wb, "NonScaled", densiteNonScaledData())
      addWorksheet(wb, "Scaled")
      writeData(wb, "Scaled", densiteScaledData())
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # PCA Proteomique
  
  pcaProteoPrcomp <- eventReactive(input$goPCAProteo, {
    req(proteoData())
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    prcomp(proteo_filtered, scale. = TRUE)
  })
  
  pcaProteoOpls <- eventReactive(input$goPCAProteo, {
    req(proteoData())
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    opls(x = proteo_filtered, predI = 9, crossvalI = 7, permI = 100)
  })
  
  output$pcaProteoScree <- renderPlotly({
    req(pcaProteoPrcomp())
    df <- data.frame(variance = pcaProteoPrcomp()$sdev^2)
    df$Component <- factor(seq_len(nrow(df)))
    p <- ggplot(df, aes(x = Component, y = variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_light() +
      labs(title = "Scree Plot - PCA Proteomique", x = "Composantes", y = "Variance expliquée")
    ggplotly(p)
  })
  
  output$pcaProteoScores12 <- renderPlotly({
    req(pcaProteoOpls(), pcaProteoPrcomp(), metadataData())
    scores <- as.data.frame(pcaProteoOpls()@scoreMN)
    scores$sample_name <- rownames(scores)
    
    # merge avec metadata
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    pc_var <- (pcaProteoPrcomp()$sdev^2) / sum(pcaProteoPrcomp()$sdev^2) * 100
    p <- ggplot(merged_scores, aes(x = p1, y = p2,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PCA Scores (Dim.1 vs Dim.2)",
           x = sprintf("PC1 (%.1f%%)", pc_var[1]),
           y = sprintf("PC2 (%.1f%%)", pc_var[2])) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$pcaProteoScores34 <- renderPlotly({
    req(pcaProteoOpls(), pcaProteoPrcomp(), metadataData())
    scores <- as.data.frame(pcaProteoOpls()@scoreMN)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    pc_var <- (pcaProteoPrcomp()$sdev^2) / sum(pcaProteoPrcomp()$sdev^2) * 100
    p <- ggplot(merged_scores, aes(x = p3, y = p4,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PCA Scores (Dim.3 vs Dim.4)",
           x = sprintf("PC3 (%.1f%%)", pc_var[3]),
           y = sprintf("PC4 (%.1f%%)", pc_var[4])) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$pcaProteoLoadings <- renderPlotly({
    req(pcaProteoPrcomp())
    loadings <- as.data.frame(pcaProteoPrcomp()$rotation)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = PC1, y = PC2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PCA Loadings - Proteomique", x = "PC1", y = "PC2") +
      theme_light()
    ggplotly(p)
  })
  
  output$pcaProteoSummary <- renderPrint({
    req(pcaProteoPrcomp())
    summary(pcaProteoPrcomp())
  })
  
  # Download PCA Proteo
  output$download_PCAProteo <- downloadHandler(
    filename = function() { paste("PCA_Proteomique_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      # Scree
      scree_df <- data.frame(variance = pcaProteoPrcomp()$sdev^2)
      scree_df$Component <- factor(seq_len(nrow(scree_df)))
      addWorksheet(wb, "Scree")
      writeData(wb, "Scree", scree_df)
      
      # Loadings
      loadings <- as.data.frame(pcaProteoPrcomp()$rotation)
      loadings$variable <- rownames(loadings)
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      # Scores via opls
      scores <- as.data.frame(pcaProteoOpls()@scoreMN)
      scores$sample_name <- rownames(scores)
      scores$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      # Résumé
      res_text <- capture.output(summary(pcaProteoPrcomp()))
      addWorksheet(wb, "Résumé")
      writeData(wb, "Résumé", data.frame(Text = res_text))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # PCA Métabolomique
  
  pcaMetaboResult <- eventReactive(input$goPCAMetabo, {
    req(metaboData(), metadataData())
    dt <- metaboData()
    metabo_filtered <- cleanNumericData(dt, idCol = 1)
    opls(x = metabo_filtered, predI = 9, crossvalI = 7, permI = 100)
  })
  output$pcaMetaboScree <- renderPlotly({
    req(pcaMetaboResult())
    df <- data.frame(variance = pcaMetaboResult()@pcaVarVn)
    df$Dim <- rownames(df)
    p <- ggplot(df, aes(x = Dim, y = variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_light() +
      labs(title = "Scree Plot - PCA Métabolomique", x = "Composantes", y = "Variance expliquée")
    ggplotly(p)
  })
  output$pcaMetaboLoadings <- renderPlotly({
    req(pcaMetaboResult())
    loadings <- as.data.frame(pcaMetaboResult()@loadingMN)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = p1, y = p2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PCA Loadings - Métabolomique", x = "p1", y = "p2") +
      theme_light()
    ggplotly(p)
  })
  output$pcaMetaboSummary <- renderPrint({
    req(pcaMetaboResult())
    capture.output(summary(pcaMetaboResult()))
  })
  output$pcaMetaboScores12 <- renderPlotly({
    req(pcaMetaboResult(), metadataData())
    scores <- as.data.frame(pcaMetaboResult()@scoreMN)
    scores$sample_name <- rownames(scores)
    scores$subclass1 <- metadataData()$subclass1
    scores$p1 <- as.numeric(scores$p1)
    scores$p2 <- as.numeric(scores$p2)
    var_vec <- pcaMetaboResult()@pcaVarVn
    percent1 <- round((var_vec["p1"] / sum(var_vec)) * 100, 1)
    percent2 <- round((var_vec["p2"] / sum(var_vec)) * 100, 1)
    p <- ggplot(scores, aes(x = p1, y = p2, color = subclass1, text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PCA Scores (Dim.1 vs Dim.2)",
           x = paste0("p1 (", percent1, "%)"),
           y = paste0("p2 (", percent2, "%)")) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  output$pcaMetaboScores34 <- renderPlotly({
    req(pcaMetaboResult(), metadataData())
    scores <- as.data.frame(pcaMetaboResult()@scoreMN)
    scores$sample_name <- rownames(scores)
    scores$subclass1 <- metadataData()$subclass1
    scores$p3 <- as.numeric(scores$p3)
    scores$p4 <- as.numeric(scores$p4)
    var_vec <- pcaMetaboResult()@pcaVarVn
    percent3 <- round((var_vec["p3"] / sum(var_vec)) * 100, 1)
    percent4 <- round((var_vec["p4"] / sum(var_vec)) * 100, 1)
    p <- ggplot(scores, aes(x = p3, y = p4, color = subclass1, text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PCA Scores (Dim.3 vs Dim.4)",
           x = paste0("p3 (", percent3, "%)"),
           y = paste0("p4 (", percent4, "%)")) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  output$download_PCAMetabo <- downloadHandler(
    filename = function() { paste("PCA_Metabolomique_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      scree_df <- data.frame(variance = pcaMetaboResult()@pcaVarVn)
      scree_df$Dim <- rownames(scree_df)
      addWorksheet(wb, "Scree")
      writeData(wb, "Scree", scree_df)
      
      loadings <- as.data.frame(pcaMetaboResult()@loadingMN)
      loadings$variable <- rownames(loadings)
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      scores <- as.data.frame(pcaMetaboResult()@scoreMN)
      scores$sample_name <- rownames(scores)
      scores$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      res_text <- capture.output(summary(pcaMetaboResult()))
      addWorksheet(wb, "Résumé")
      writeData(wb, "Résumé", data.frame(Text = res_text))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  
  #  PLS‑DA Proteomique
  
  plsdaProteoResult <- eventReactive(input$goPLSDAProteo, {
    req(proteoData(), metadataData())
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    group <- as.factor(metadataData()$subclass1)  # variable de regroupement
    plsda(proteo_filtered, group, ncomp = 2)
  })
  
  #  Calcul de la performance via validation croisée
  plsdaProteo_perf <- reactive({
    req(plsdaProteoResult())
    perf(plsdaProteoResult(), 
         validation = "Mfold", 
         folds = min(7, min(table(metadataData()$subclass1))), 
         nrepeat = 10, 
         progressBar = FALSE)
  })
  
  output$plsdaProteoSummary <- renderPrint({
    req(plsdaProteoResult())
    print(plsdaProteoResult())
  })
  
  output$plsdaProteoPerfPlot <- renderPlot({
    req(plsdaProteo_perf())
    perf_res <- plsdaProteo_perf()
    q2_vals <- as.numeric(as.character(perf_res$Q2))
    q2_df <- data.frame(Comp = seq_along(q2_vals), Q2 = q2_vals)
    ggplot(q2_df, aes(x = Comp, y = Q2)) +
      geom_line() +
      geom_point() +
      labs(title = "Performance PLS‑DA (Q²)",
           x = "Composante", y = "Q²") +
      theme_minimal()
  })
  
  output$plsdaProteoScores <- renderPlotly({
    req(plsdaProteoResult(), metadataData())
    res <- plsdaProteoResult()
    scores <- as.data.frame(res$variates$X)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(res$explained_variance)) {
      ev <- res$explained_variance$X
      percent1 <- round(ev[1] * 100, 1)
      ylab <- paste0("Comp1 (", percent1, "%)")
    } else {
      ylab <- "Comp1"
    }
    p <- ggplot(merged_scores, aes(x = sample_name, y = comp1,
                                   fill = .data[["subclass1"]], text = sample_name)) +
      geom_bar(stat = "identity") +
      labs(title = "PLS‑DA Scores (Proteomique)",
           x = "Sample", y = ylab) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggplotly(p, tooltip = "text")
  })
  
  output$plsdaProteoLoadingsBar1 <- renderPlotly({
    req(plsdaProteoResult())
    res <- plsdaProteoResult()
    load_df <- as.data.frame(res$loadings$X)
    load_df$variable <- rownames(load_df)
    load_df <- load_df[order(load_df$comp1), ]
    load_df$rank <- seq_len(nrow(load_df))
    
    p <- ggplot(load_df, aes(x = rank, y = comp1, label = variable)) +
      geom_bar(stat = "identity") +
      geom_text(aes(angle = 90,
                    y = ifelse(comp1 > 0, comp1 + 0.002, comp1 - 0.002)),
                size = 2, hjust = ifelse(load_df$comp1 < 0, 1, 0)) +
      labs(title = "Loadings plot on Dim.1 (PLS‑DA Proteo)") +
      theme_light()
    ggplotly(p, tooltip = "label")
  })
  
  # PLS‑DA Protéomique via opls (ropls)
  
  # 1. Définition du modèle OPLS
  oplsProteoResult <- eventReactive(input$goOPLSProteo, {
    req(proteoData(), metadataData())
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    group <- as.factor(metadataData()$subclass1)
    # Calcul du modèle OPLS (nombre de composantes prédictives estimé automatiquement et 100 permutations)
    opls(x = proteo_filtered, y = group, predI = NA, permI = 100)
  })
  
  # 2. Affichage Q²
  output$oplsProteoPerfPlot <- renderPlotly({
    req(oplsProteoResult())
    res <- oplsProteoResult()
    q2_value <- as.numeric(res@summaryDF["Q2(cum)", ])
    q2_df <- data.frame(Metric = "Q2(cum)", Value = q2_value)
    
    p <- ggplot(q2_df, aes(x = Metric, y = Value)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Performance OPLS - Q2(cum) (Protéomique)",
           x = "", y = "Q2") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$oplsProteoScores <- renderPlotly({
    req(oplsProteoResult())
    res <- oplsProteoResult()
    scores <- as.data.frame(res@scoreMN)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    p <- ggplot(merged_scores, aes(x = p1, y = p2, color = .data[["subclass1"]], text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "Scores OPLS - Protéomique", x = "p1", y = "p2") +
      theme_light()
    
    ggplotly(p, tooltip = "text")
  })
  
  output$oplsProteoLoadings <- renderPlotly({
    req(oplsProteoResult())
    res <- oplsProteoResult()
    loadings <- as.data.frame(res@loadingMN)
    loadings$variable <- rownames(loadings)
    
    p <- ggplot(loadings, aes(x = p1, y = p2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "Loadings OPLS - Protéomique", x = "p1", y = "p2") +
      theme_light()
    
    ggplotly(p)
  })
  
  # Plot VIP
  output$plsdaProteoVIP <- renderPlotly({
    req(plsdaProteoResult())
    res <- plsdaProteoResult()
    vip_scores <- vip(res)[, 1]
    vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
    vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), ]
    vip_df$rank <- seq_len(nrow(vip_df))
    
    p <- ggplot(vip_df, aes(x = rank, y = VIP, label = variable)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "VIP Plot - PLS‑DA Proteomique", x = "Rank", y = "VIP") +
      theme_light()
    ggplotly(p)
  })
  
  # Plot loadings vs VIP
  output$plsdaProteoLoadingsVIP <- renderPlotly({
    req(plsdaProteoResult())
    res <- plsdaProteoResult()
    loadings <- data.frame(variable = rownames(res$loadings$X),
                           loading = res$loadings$X[, 1])
    vip_scores <- vip(res)[, 1]
    if (is.null(names(vip_scores)) || length(vip_scores) == 0) {
      names(vip_scores) <- rownames(res$loadings$X)
    }
    vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
    merged_df <- merge(loadings, vip_df, by = "variable")
    
    p <- ggplot(merged_df, aes(x = loading, y = VIP, label = variable)) +
      geom_point(color = "steelblue", size = 3) +
      geom_text_repel() +
      labs(title = "Loadings vs VIP - PLS‑DA Proteomique",
           x = "Loadings (Comp1)", y = "VIP") +
      theme_light()
    ggplotly(p)
  })
  
  # 9. Téléchargement des résultats
  output$download_PLSDAProteo <- downloadHandler(
    filename = function() { paste("PLSDA_Proteomique_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- plsdaProteoResult()
      
      # Scores
      scores <- as.data.frame(res$variates$X)
      scores$sample_name <- rownames(scores)
      scores$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      # VIP
      vip_scores <- vip(res)[, 1]
      vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
      vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), ]
      vip_df$rank <- seq_len(nrow(vip_df))
      addWorksheet(wb, "VIP")
      writeData(wb, "VIP", vip_df)
      
      # Loadings
      loadings <- data.frame(variable = rownames(res$loadings$X),
                             loading = res$loadings$X[, 1])
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      # VIP vs Loadings
      merged_df <- merge(loadings, vip_df, by = "variable")
      addWorksheet(wb, "VIP_vs_Loadings")
      writeData(wb, "VIP_vs_Loadings", merged_df)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # 8) PLS‑DA Métabolomique
  
  plsdaMetaboResult <- eventReactive(input$goPLSDAMetabo, {
    req(metaboData(), metadataData())
    metabo_filtered <- cleanNumericData(metaboData(), idCol = 1)
    group <- as.factor(metadataData()$subclass1)
    plsda(metabo_filtered, group, ncomp = 2)
  })
  
  output$plsdaMetaboSummary <- renderPrint({
    req(plsdaMetaboResult())
    capture.output(print(plsdaMetaboResult()))
  })
  
  output$plsdaMetaboScores <- renderPlotly({
    req(plsdaMetaboResult(), metadataData())
    res <- plsdaMetaboResult()
    scores <- as.data.frame(res$variates$X)
    scores$sample_name <- rownames(scores)
    
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(res$explained_variance)) {
      ev <- res$explained_variance$X
      percent1 <- round(ev[1] * 100, 1)
      ylab <- paste0("Comp1 (", percent1, "%)")
    } else {
      ylab <- "Comp1"
    }
    p <- ggplot(merged_scores, aes(x = sample_name, y = comp1,
                                   fill = .data[["subclass1"]], text = sample_name)) +
      geom_bar(stat = "identity") +
      labs(title = "PLS‑DA Scores (Métabolomique)",
           x = "Sample", y = ylab) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggplotly(p, tooltip = "text")
  })
  
  output$plsdaMetaboLoadingsBar1 <- renderPlotly({
    req(plsdaMetaboResult())
    res <- plsdaMetaboResult()
    load_df <- as.data.frame(res$loadings$X)
    load_df$variable <- rownames(load_df)
    load_df <- load_df[order(load_df$comp1), ]
    load_df$rank <- seq_len(nrow(load_df))
    
    p <- ggplot(load_df, aes(x = rank, y = comp1, label = variable)) +
      geom_bar(stat = "identity") +
      geom_text(aes(angle = 90,
                    y = ifelse(comp1 > 0, comp1 + 0.002, comp1 - 0.002)),
                size = 2, hjust = ifelse(load_df$comp1 < 0, 1, 0)) +
      labs(title = "Loadings plot on Dim.1 (PLS‑DA Metabo)") +
      theme_light()
    ggplotly(p, tooltip = "label")
  })
  
  output$plsdaMetaboVIP <- renderPlotly({
    req(plsdaMetaboResult())
    res <- plsdaMetaboResult()
    vip_scores <- vip(res)[, 1]
    vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
    vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), ]
    vip_df$rank <- seq_len(nrow(vip_df))
    p <- ggplot(vip_df, aes(x = rank, y = VIP, label = variable)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "VIP Plot - PLS‑DA Métabolomique", x = "Rang", y = "VIP") +
      theme_light()
    ggplotly(p)
  })
  
  output$plsdaMetaboLoadingsVIP <- renderPlotly({
    req(plsdaMetaboResult())
    res <- plsdaMetaboResult()
    loadings <- data.frame(variable = rownames(res$loadings$X),
                           loading = res$loadings$X[, 1])
    vip_scores <- vip(res)[, 1]
    if (is.null(names(vip_scores)) || length(vip_scores) == 0) {
      names(vip_scores) <- rownames(res$loadings$X)
    }
    vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
    merged_df <- merge(loadings, vip_df, by = "variable")
    p <- ggplot(merged_df, aes(x = loading, y = VIP, label = variable)) +
      geom_point(color = "steelblue", size = 3) +
      geom_text_repel() +
      labs(title = "Loadings vs VIP - PLS‑DA Métabolomique",
           x = "Loadings (Comp1)", y = "VIP") +
      theme_light()
    ggplotly(p)
  })
  
  # Download PLSDA Metabo
  output$download_PLSDAMetabo <- downloadHandler(
    filename = function() { paste("PLSDA_Metabolomique_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- plsdaMetaboResult()
      
      # Scores
      scores <- as.data.frame(res$variates$X)
      scores$sample_name <- rownames(scores)
      scores$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      # VIP
      vip_scores <- vip(res)[, 1]
      vip_df <- data.frame(variable = names(vip_scores), VIP = as.vector(vip_scores))
      vip_df <- vip_df[order(vip_df$VIP, decreasing = TRUE), ]
      vip_df$rank <- seq_len(nrow(vip_df))
      addWorksheet(wb, "VIP")
      writeData(wb, "VIP", vip_df)
      
      # Loadings
      loadings <- data.frame(variable = rownames(res$loadings$X),
                             loading = res$loadings$X[, 1])
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      # VIP vs Loadings
      merged_df <- merge(loadings, vip_df, by = "variable")
      addWorksheet(wb, "VIP_vs_Loadings")
      writeData(wb, "VIP_vs_Loadings", merged_df)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  #  PLS‑DA Métabolomique via opls (ropls)
  
  # Définition du modèle OPLS
  oplsMetaboResult <- eventReactive(input$goOPLSMetabo, {
    req(metaboData(), metadataData())
    metabo_filtered <- cleanNumericData(metaboData(), idCol = 1)
    group <- as.factor(metadataData()$subclass1)
    # Calcul du modèle OPLS (nombre de composantes prédictives estimé automatiquement et 100 permutations)
    opls(x = metabo_filtered, y = group, predI = NA, permI = 100)
  })
  
  # Affichage Q²
  output$oplsMetaboPerfPlot <- renderPlotly({
    req(oplsMetaboResult())
    res <- oplsMetaboResult()
    q2_value <- as.numeric(res@summaryDF["Q2(cum)", ])
    q2_df <- data.frame(Metric = "Q2(cum)", Value = q2_value)
    
    p <- ggplot(q2_df, aes(x = Metric, y = Value)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Performance OPLS - Q2(cum) (Métabolomique)",
           x = "", y = "Q2") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$oplsMetaboScores <- renderPlotly({
    req(oplsMetaboResult())
    res <- oplsMetaboResult()
    scores <- as.data.frame(res@scoreMN)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    p <- ggplot(merged_scores, aes(x = p1, y = p2, color = .data[["subclass1"]], text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "Scores OPLS - Métabolomique", x = "p1", y = "p2") +
      theme_light()
    
    ggplotly(p, tooltip = "text")
  })
  
  output$oplsMetaboLoadings <- renderPlotly({
    req(oplsMetaboResult())
    res <- oplsMetaboResult()
    loadings <- as.data.frame(res@loadingMN)
    loadings$variable <- rownames(loadings)
    
    p <- ggplot(loadings, aes(x = p1, y = p2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "Loadings OPLS - Métabolomique", x = "p1", y = "p2") +
      theme_light()
    
    ggplotly(p)
  })
  
  
  # PLS Regression
  
  plsRegResult <- eventReactive(input$goPLSReg, {
    req(proteoData(), metaboData())
    proteo_numeric <- cleanNumericData(proteoData(), idCol = 1)
    metabo_numeric <- cleanNumericData(metaboData(), idCol = 1)
    pls(X = as.matrix(proteo_numeric), Y = as.matrix(metabo_numeric),
        ncomp = 2, scale = TRUE, mode = "regression")
  })
  
  output$plsRegSummary <- renderPrint({
    req(plsRegResult())
    capture.output(print(plsRegResult()))
  })
  
  output$plsRegScoresProteo <- renderPlotly({
    req(plsRegResult(), metadataData())
    scores <- as.data.frame(plsRegResult()$variates$X)
    scores$sample_name <- rownames(scores)
    
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    ev <- plsRegResult()$explained_variance$X
    percent1 <- round(ev[1] * 100, 1)
    percent2 <- round(ev[2] * 100, 1)
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PLS Regression Scores - Proteomique",
           x = paste0("comp1 (", percent1, "%)"),
           y = paste0("comp2 (", percent2, "%)")) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$plsRegScoresMetabo <- renderPlotly({
    req(plsRegResult(), metadataData())
    scores <- as.data.frame(plsRegResult()$variates$Y)
    scores$sample_name <- rownames(scores)
    
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    ev <- plsRegResult()$explained_variance$Y
    percent1 <- round(ev[1] * 100, 1)
    percent2 <- round(ev[2] * 100, 1)
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PLS Regression Scores - Métabolomique",
           x = paste0("comp1 (", percent1, "%)"),
           y = paste0("comp2 (", percent2, "%)")) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$plsRegLoadingsProteo <- renderPlotly({
    req(plsRegResult())
    loadings <- data.frame(plsRegResult()$loadings$X)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = comp1, y = comp2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PLS Regression Loadings - Proteomique", x = "comp1", y = "comp2") +
      theme_light()
    ggplotly(p)
  })
  
  output$plsRegLoadingsMetabo <- renderPlotly({
    req(plsRegResult())
    loadings <- data.frame(plsRegResult()$loadings$Y)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = comp1, y = comp2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PLS Regression Loadings - Métabolomique", x = "comp1", y = "comp2") +
      theme_light()
    ggplotly(p)
  })
  
  # Download PLS Regression
  output$download_PLSReg <- downloadHandler(
    filename = function() { paste("PLS_Regression_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- plsRegResult()
      
      scoresX <- as.data.frame(res$variates$X)
      scoresX$sample_name <- rownames(scoresX)
      scoresX$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores_Proteo")
      writeData(wb, "Scores_Proteo", scoresX)
      
      scoresY <- as.data.frame(res$variates$Y)
      scoresY$sample_name <- rownames(scoresY)
      scoresY$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores_Metabo")
      writeData(wb, "Scores_Metabo", scoresY)
      
      loadingsX <- data.frame(res$loadings$X)
      loadingsX$variable <- rownames(loadingsX)
      addWorksheet(wb, "Loadings_Proteo")
      writeData(wb, "Loadings_Proteo", loadingsX)
      
      loadingsY <- data.frame(res$loadings$Y)
      loadingsY$variable <- rownames(loadingsY)
      addWorksheet(wb, "Loadings_Metabo")
      writeData(wb, "Loadings_Metabo", loadingsY)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # PLS Canonical
  
  plsCanoResult <- eventReactive(input$goPLSCano, {
    req(proteoData(), metaboData())
    proteo_numeric <- cleanNumericData(proteoData(), idCol = 1)
    metabo_numeric <- cleanNumericData(metaboData(), idCol = 1)
    # Définition X et Y
    pls(X = as.matrix(proteo_numeric), Y = as.matrix(metabo_numeric),
        ncomp = 2, scale = TRUE, mode = "canonical")
  })
  
  output$plsCanoSummary <- renderPrint({
    req(plsCanoResult())
    capture.output(print(plsCanoResult()))
  })
  
  output$plsCanoScoresProteo <- renderPlotly({
    req(plsCanoResult(), metadataData())
    scores <- as.data.frame(plsCanoResult()$variates$X)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(plsCanoResult()$explained_variance)) {
      ev <- plsCanoResult()$explained_variance$X
      percent1 <- round(ev[1] * 100, 1)
      percent2 <- round(ev[2] * 100, 1)
      xlab <- paste0("comp1 (", percent1, "%)")
      ylab <- paste0("comp2 (", percent2, "%)")
    } else {
      xlab <- "comp1"
      ylab <- "comp2"
    }
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PLS Canonical Scores - Proteomique",
           x = xlab, y = ylab) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$plsCanoScoresMetabo <- renderPlotly({
    req(plsCanoResult(), metadataData())
    scores <- as.data.frame(plsCanoResult()$variates$Y)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(plsCanoResult()$explained_variance)) {
      ev <- plsCanoResult()$explained_variance$Y
      percent1 <- round(ev[1] * 100, 1)
      percent2 <- round(ev[2] * 100, 1)
      xlab <- paste0("comp1 (", percent1, "%)")
      ylab <- paste0("comp2 (", percent2, "%)")
    } else {
      xlab <- "comp1"
      ylab <- "comp2"
    }
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[["subclass1"]],
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "PLS Canonical Scores - Métabolomique",
           x = xlab, y = ylab) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$plsCanoLoadingsProteo <- renderPlotly({
    req(plsCanoResult())
    loadings <- data.frame(plsCanoResult()$loadings$X)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = comp1, y = comp2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PLS Canonical Loadings - Proteomique", x = "comp1", y = "comp2") +
      theme_light()
    ggplotly(p)
  })
  
  output$plsCanoLoadingsMetabo <- renderPlotly({
    req(plsCanoResult())
    loadings <- data.frame(plsCanoResult()$loadings$Y)
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = comp1, y = comp2, label = variable)) +
      geom_point(size = 3) +
      geom_text_repel() +
      labs(title = "PLS Canonical Loadings - Métabolomique", x = "comp1", y = "comp2") +
      theme_light()
    ggplotly(p)
  })
  
  # Download PLS Canonical
  output$download_PLSCano <- downloadHandler(
    filename = function() { paste("PLS_Canonical_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- plsCanoResult()
      
      # Scores X
      scoresX <- as.data.frame(res$variates$X)
      scoresX$sample_name <- rownames(scoresX)
      scoresX$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores_Proteo")
      writeData(wb, "Scores_Proteo", scoresX)
      
      # Scores Y
      scoresY <- as.data.frame(res$variates$Y)
      scoresY$sample_name <- rownames(scoresY)
      scoresY$subclass1 <- metadataData()$subclass1
      addWorksheet(wb, "Scores_Metabo")
      writeData(wb, "Scores_Metabo", scoresY)
      
      # Loadings X
      loadingsX <- data.frame(res$loadings$X)
      loadingsX$variable <- rownames(loadingsX)
      addWorksheet(wb, "Loadings_Proteo")
      writeData(wb, "Loadings_Proteo", loadingsX)
      
      # Loadings Y
      loadingsY <- data.frame(res$loadings$Y)
      loadingsY$variable <- rownames(loadingsY)
      addWorksheet(wb, "Loadings_Metabo")
      writeData(wb, "Loadings_Metabo", loadingsY)
      
      # Résumé
      res_text <- capture.output(print(res))
      addWorksheet(wb, "Résumé")
      writeData(wb, "Résumé", data.frame(Text = res_text))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  
  # ComDim
  
  comdimResult <- eventReactive(input$goComDim, {
    req(proteoData(), metaboData())
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    metabo_filtered <- cleanNumericData(metaboData(), idCol = 1)
    multiblock_data <- cbind(proteo_filtered, metabo_filtered)
    n_group <- c(ncol(proteo_filtered), ncol(metabo_filtered))
    ComDim(X = multiblock_data,
           block = n_group,
           name.block = c("Proteomics", "Metabolomics"),
           scale = TRUE,
           scale.block = TRUE)
  })
  
  output$comdimSaliences <- renderPlotly({
    req(comdimResult())
    res <- comdimResult()
    saliences <- res$saliences
    rownames(saliences) <- c("Proteomics", "Metabolomics")
    sal_df <- as.data.frame(t(saliences[, 1:4]))
    sal_df$Dim <- rownames(sal_df)
    saliences_melt <- melt(sal_df)
    p <- ggplot(saliences_melt, aes(x = Dim, y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_light() +
      labs(title = "Saliences - ComDim", x = "Global Components", y = "Weights")
    ggplotly(p)
  })
  
  output$comdimScores <- renderPlotly({
    req(comdimResult(), metadataData())
    res <- comdimResult()
    scores <- as.data.frame(res$Scor.g)
    scores$sample_name <- rownames(scores)
    
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    p <- ggplot(merged_scores, aes(x = Dim.1, y = Dim.2,
                                   color = .data[["subclass1"]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "ComDim Scores", x = "Dim.1", y = "Dim.2") +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$comdimLoadings <- renderPlotly({
    req(comdimResult())
    res <- comdimResult()
    loadings <- as.data.frame(res$Load.g)
    # Indiquer le block
    proteo_filtered <- cleanNumericData(proteoData(), idCol = 1)
    metabo_filtered <- cleanNumericData(metaboData(), idCol = 1)
    loadings$Block <- c(rep("Proteomics", ncol(proteo_filtered)),
                        rep("Metabolomics", ncol(metabo_filtered)))
    loadings$Variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = Dim.1, y = Dim.2, color = Block, text = Variable)) +
      geom_point(size = 3) +
      labs(title = "ComDim Loadings", x = "Dim.1", y = "Dim.2") +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$comdimSummary <- renderPrint({
    req(comdimResult())
    comdimResult() 
  })
  
  # Download ComDim
  output$download_ComDim <- downloadHandler(
    filename = function() { paste("ComDim_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- comdimResult()
      
      # Saliences
      saliences <- as.data.frame(res$saliences)
      addWorksheet(wb, "Saliences")
      writeData(wb, "Saliences", saliences)
      
      # Scores
      scores <- as.data.frame(res$Scor.g)
      scores$sample_name <- rownames(scores)
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      # Loadings
      loadings <- as.data.frame(res$Load.g)
      loadings$variable <- rownames(loadings)
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      # Résumé
      res_text <- capture.output(print(res))
      addWorksheet(wb, "Résumé")
      writeData(wb, "Résumé", data.frame(Text = res_text))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # Block PLS‑DA (mixOmics)
  
  blockPLSDA_res <- eventReactive(input$goBlockPLSDA, {
    req(proteoData(), metaboData(), metadataData(), input$groupVarBlock)
    proteo_mat <- as.matrix(cleanNumericData(proteoData(), idCol = 1))
    metabo_mat <- as.matrix(cleanNumericData(metaboData(), idCol = 1))
    
    # Intersection
    common_samples <- intersect(rownames(proteo_mat), rownames(metabo_mat))
    proteo_mat <- proteo_mat[common_samples, , drop = FALSE]
    metabo_mat <- metabo_mat[common_samples, , drop = FALSE]
    
    block_data <- list(proteo = proteo_mat, metabo = metabo_mat)
    group <- as.factor(metadataData()[[input$groupVarBlock]])
    block.plsda(X = block_data, Y = group, design = "full",
                ncomp = input$ncompBlock, all.outputs = TRUE, near.zero.var = TRUE)
  })
  
  output$blockPLS_perf <- renderPlotly({
    req(blockPLSDA_res())
    md <- metadataData()
    group <- as.factor(md[[input$groupVarBlock]])
    nfolds <- min(7, min(table(group)))
    perf_res <- tryCatch({
      perf(blockPLSDA_res(), validation = 'Mfold', folds = nfolds, nrepeat = 1, 
           auc = TRUE, progressBar = FALSE)
    }, error = function(e) NULL)
    if (is.null(perf_res) || is.null(perf_res$Q2)) return(plotly_empty())
    df <- as.data.frame(perf_res$Q2)
    df$comp <- as.numeric(rownames(df))
    p <- ggplot(df, aes(x = comp, y = Q2)) +
      geom_line() +
      geom_point() +
      labs(title = "Performance (Q2) par composante", x = "Composante", y = "Q2") +
      theme_light()
    ggplotly(p)
  })
  
  output$blockPLS_expl <- renderPlotly({
    req(blockPLSDA_res())
    ave_list <- blockPLSDA_res()$AVE$AVE_X[1:2]
    if (length(ave_list) < 2) return(plotly_empty())
    ave_df <- tryCatch(bind_rows(ave_list, .id = "block"), error = function(e) NULL)
    if (is.null(ave_df)) return(plotly_empty())
    outer <- blockPLSDA_res()$AVE[["AVE_outer"]]
    outer_df <- data.frame(block = "outer_model", t(outer))
    all_df <- merge(ave_df, outer_df, by = "block", all = TRUE)
    all_df <- melt(all_df, id.vars = "block", variable.name = "comp", value.name = "value")
    p <- ggplot(all_df, aes(x = comp, y = value, fill = block)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Variance expliquée par bloc", x = "Composante", y = "% de variance expliquée") +
      theme_light()
    ggplotly(p)
  })
  
  output$blockPLS_scores <- renderPlotly({
    req(blockPLSDA_res(), metadataData())
    scores <- as.data.frame(blockPLSDA_res()$variates$weighted.average)
    if (nrow(scores) == 0) {
      if (!is.null(blockPLSDA_res()$variates$proteo) && !is.null(blockPLSDA_res()$variates$metabo)) {
        scores <- (blockPLSDA_res()$variates$proteo + blockPLSDA_res()$variates$metabo) / 2
      } else return(plotly_empty())
    }
    scores <- as.data.frame(scores)
    scores$sample_name <- rownames(scores)
    
    md <- metadataData()
    merged_scores <- merge(scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(blockPLSDA_res()$explained_variance)) {
      ev <- blockPLSDA_res()$explained_variance$weighted.average
      percent1 <- round(ev[1] * 100, 1)
      percent2 <- round(ev[2] * 100, 1)
      xlab <- paste0("comp1 (", percent1, "%)")
      ylab <- paste0("comp2 (", percent2, "%)")
    } else {
      xlab <- "comp1"
      ylab <- "comp2"
    }
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[[input$groupVarBlock]], 
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "Scores - Block PLS‑DA (Dim.1 vs Dim.2)", x = xlab, y = ylab) +
      theme_light()
    ggplotly(p, tooltip = "text")
  })
  
  output$blockPLS_weighted <- renderPlotly({
    req(blockPLSDA_res(), metadataData())
    bp_res <- blockPLSDA_res()
    blockPLS_variates.weighted <- bp_res$variates[c("proteo", "metabo")]
    for (omic in c("proteo", "metabo")) {
      for (comp in c("comp1", "comp2")) {
        blockPLS_variates.weighted[[omic]][, comp] <-
          blockPLS_variates.weighted[[omic]][, comp] * bp_res$weights[omic, comp]
      }
    }
    weighted_scores <- abind::abind(blockPLS_variates.weighted[c("proteo", "metabo")], along = 3)
    weighted_scores <- apply(weighted_scores, c(1, 2), mean)
    weighted_scores <- as.data.frame(weighted_scores)
    weighted_scores$sample_name <- rownames(weighted_scores)
    
    md <- metadataData()
    merged_scores <- merge(weighted_scores, md, by.x = "sample_name", by.y = "sample_name", all.x = TRUE)
    
    if (!is.null(bp_res$explained_variance)) {
      ev <- bp_res$explained_variance$weighted.average
      percent1 <- round(ev[1] * 100, 1)
      percent2 <- round(ev[2] * 100, 1)
      xlab <- paste0("comp1 (", percent1, "%)")
      ylab <- paste0("comp2 (", percent2, "%)")
    } else {
      xlab <- "comp1"
      ylab <- "comp2"
    }
    p <- ggplot(merged_scores, aes(x = comp1, y = comp2,
                                   color = .data[[input$groupVarBlock]])) +
      geom_point(size = 4) +
      labs(title = "Scores pondérés (comp1 vs comp2)", x = xlab, y = ylab) +
      theme_light()
    ggplotly(p)
  })
  
  output$blockPLS_corrCircle <- renderPlotly({
    req(blockPLSDA_res())
    loadings_proteo <- as.data.frame(blockPLSDA_res()$loadings$proteo)
    loadings_metabo <- as.data.frame(blockPLSDA_res()$loadings$metabo)
    blockPLS_loadings <- rbind.data.frame(loadings_proteo, loadings_metabo)
    blockPLS_loadings$omic <- c(rep("proteo", nrow(loadings_proteo)),
                                rep("metabo", nrow(loadings_metabo)))
    blockPLS_loadings$variable <- rownames(blockPLS_loadings)
    p <- ggplot(blockPLS_loadings, aes(x = comp1, y = comp2,
                                       color = omic, label = variable)) +
      geom_point() +
      geom_text_repel() +
      labs(title = "Cercle de corrélation", x = "Comp1", y = "Comp2") +
      theme_light()
    ggplotly(p)
  })
  
  output$blockPLS_loadings <- renderPlotly({
    req(blockPLSDA_res())
    loadings_proteo <- as.data.frame(blockPLSDA_res()$loadings$proteo)
    loadings_proteo$variable <- rownames(loadings_proteo)
    loadings_proteo$block <- "proteo"
    loadings_metabo <- as.data.frame(blockPLSDA_res()$loadings$metabo)
    loadings_metabo$variable <- rownames(loadings_metabo)
    loadings_metabo$block <- "metabo"
    loadings_all <- rbind(loadings_proteo, loadings_metabo)
    p <- ggplot(loadings_all, aes(x = comp1, y = comp2,
                                  color = block, label = variable)) +
      geom_point(size = 3) +
      geom_text(aes(label = variable)) +
      labs(title = "Loadings - Block PLS‑DA", x = "Comp1", y = "Comp2") +
      theme_light()
    ggplotly(p)
  })
  
  output$blockPLS_permtest <- renderPrint({
    req(blockPLSDA_res())
    test_res <- DIABLO.test(blockPLSDA_res(), progress = FALSE)
    test_res
  })
  
  output$blockPLS_summary <- renderPrint({
    req(blockPLSDA_res())
    capture.output(print(blockPLSDA_res()))
  })
  
  # Download block PLS-DA
  output$download_BlockPLSDA <- downloadHandler(
    filename = function() { paste("Block_PLSDA_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      md <- metadataData()
      group <- as.factor(md[[input$groupVarBlock]])
      nfolds <- min(7, min(table(group)))
      perf_res <- tryCatch({
        perf(blockPLSDA_res(), validation = 'Mfold', folds = nfolds, nrepeat = 1, 
             auc = TRUE, progressBar = FALSE)
      }, error = function(e) NULL)
      if (!is.null(perf_res) && !is.null(perf_res$Q2)) {
        perf_df <- as.data.frame(perf_res$Q2)
        perf_df$comp <- rownames(perf_df)
        addWorksheet(wb, "Performance")
        writeData(wb, "Performance", perf_df)
      }
      
      ave_list <- blockPLSDA_res()$AVE$AVE_X[1:2]
      if (length(ave_list) >= 2) {
        ave_df <- tryCatch(bind_rows(ave_list, .id = "block"), error = function(e) NULL)
        if (!is.null(ave_df)) {
          outer <- blockPLSDA_res()$AVE[["AVE_outer"]]
          outer_df <- data.frame(block = "outer_model", t(outer))
          all_ave <- merge(ave_df, outer_df, by = "block", all = TRUE)
          addWorksheet(wb, "Explained Variance")
          writeData(wb, "Explained Variance", all_ave)
        }
      }
      
      scores <- as.data.frame(blockPLSDA_res()$variates$weighted.average)
      scores$sample_name <- rownames(scores)
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      loadings_proteo <- as.data.frame(blockPLSDA_res()$loadings$proteo)
      loadings_proteo$variable <- rownames(loadings_proteo)
      loadings_proteo$block <- "proteo"
      loadings_metabo <- as.data.frame(blockPLSDA_res()$loadings$metabo)
      loadings_metabo$variable <- rownames(loadings_metabo)
      loadings_metabo$block <- "metabo"
      loadings_all <- rbind(loadings_proteo, loadings_metabo)
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings_all)
      
      perm_res <- capture.output(DIABLO.test(blockPLSDA_res(), progress = FALSE))
      addWorksheet(wb, "Permutation Test")
      writeData(wb, "Permutation Test", data.frame(Text = perm_res))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  # Corr Network Interactif interactif avec visNetwork
  network_plot_res <- eventReactive(input$goNetwork, {
    req(blockPLSDA_res(), input$networkCutoff)
    
    net_obj <- tryCatch({
      mixOmics::network(blockPLSDA_res(), 
                        comp = list(proteo = 1, metabo = 1),
                        cutoff = input$networkCutoff)
    }, error = function(e) {
      showNotification(paste("Erreur dans mixOmics::network :", e$message), type = "error")
      return(NULL)
    })
    
    req(net_obj)
    
    if (!"gR" %in% names(net_obj)) {
      showNotification("L'objet renvoyé par mixOmics::network ne contient pas l'élément 'gR'.", type = "error")
      return(NULL)
    }
    
    net_obj$gR
  })
  
  output$networkPlot <- renderVisNetwork({
    req(network_plot_res())
    g <- network_plot_res() 
    visIgraph(g) %>% visOptions(highlightNearest = TRUE, selectedBy = "name")
  })
  output$download_Network <- downloadHandler(
    filename = function() { paste("Network_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      # Récupérer l'objet réseau (un objet igraph)
      g <- network_plot_res()
      
      # Extraire la liste des arêtes et des nœuds
      edges <- igraph::as_data_frame(g, what = "edges")
      nodes <- igraph::as_data_frame(g, what = "vertices")
      
      addWorksheet(wb, "Edges")
      writeData(wb, "Edges", edges)
      
      addWorksheet(wb, "Nodes")
      writeData(wb, "Nodes", nodes)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # Consensus OPLS 
  
  consOPLS_res <- eventReactive(input$goConsensusOPLS, {
    req(proteoData(), metaboData(), metadataData(), input$groupVarCons)
    # Nettoyage et intersection des données
    proteo_clean <- cleanNumericData(proteoData(), idCol = 1)
    metabo_clean <- cleanNumericData(metaboData(), idCol = 1)
    common_samples <- intersect(rownames(proteo_clean), rownames(metabo_clean))
    
    proteo_clean <- proteo_clean[common_samples, , drop = FALSE]
    metabo_clean <- metabo_clean[common_samples, , drop = FALSE]
    
    COPLS_data <- list(
      proteo = as.matrix(proteo_clean),
      metabo = as.matrix(metabo_clean)
    )
    COPLS_data <- lapply(COPLS_data, function(x) {
      x <- scale(x, center = TRUE, scale = TRUE)
      if (any(is.na(x)) || any(!is.finite(x))) {
        stop("Matrice avec NA ou Inf après scaling.")
      }
      x
    })
    
    md <- metadataData()
    md <- md[md$sample_name %in% common_samples, ]
    subclass1 <- factor(md[[input$groupVarCons]])
    if (nlevels(subclass1) < 2) {
      showNotification("Le groupe sélectionné doit contenir au moins deux classes pour ConsensusOPLS.", type = "error")
      return(NULL)
    }
    dummy_subclass1 <- tryCatch(model.matrix(~ subclass1 - 1), error = function(e) NULL)
    if (is.null(dummy_subclass1) || any(is.na(dummy_subclass1))) {
      stop("La matrice dummy pour la variable de regroupement contient des valeurs manquantes.")
    }
    
    res <- tryCatch({
      ConsensusOPLS(
        data       = COPLS_data,
        Y          = dummy_subclass1,
        maxPcomp   = 1,
        maxOcomp   = 1,
        modelType  = "da",
        cvType     = "nfold",
        nfold      = min(5, nrow(COPLS_data$proteo)),
        nperm      = 5,
        verbose    = TRUE
      )
    }, error = function(e) {
      showNotification(paste("ConsensusOPLS error:", e$message), type = "error")
      return(NULL)
    })
    res
  })
  
  output$consOPLS_Q2 <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    Q2perm <- data.frame(Q2perm = consOPLS_res()@permStats$Q2Y)
    p <- ggplot(Q2perm, aes(x = Q2perm)) +
      geom_histogram(color = "grey", fill = "grey") +
      geom_vline(aes(xintercept = consOPLS_res()@Q2["po1"]),
                 color = "blue", linetype = "dashed", size = 1) +
      theme_classic() +
      ggtitle("Q2 Permutation test")
    ggplotly(p)
  })
  
  output$consOPLS_DQ2 <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    DQ2perm <- data.frame(DQ2perm = consOPLS_res()@permStats$DQ2Y)
    p <- ggplot(DQ2perm, aes(x = DQ2perm)) +
      geom_histogram(color = "grey", fill = "grey") +
      geom_vline(aes(xintercept = consOPLS_res()@DQ2["po1"]),
                 color = "blue", linetype = "dashed", size = 1) +
      theme_classic() +
      ggtitle("DQ2 Permutation test")
    ggplotly(p)
  })
  
  output$consOPLS_R2Y <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    R2Yperm <- data.frame(R2Yperm = consOPLS_res()@permStats$R2Y)
    p <- ggplot(R2Yperm, aes(x = R2Yperm)) +
      geom_histogram(color = "grey", fill = "grey") +
      geom_vline(aes(xintercept = consOPLS_res()@R2Y["po1"]),
                 color = "blue", linetype = "dashed", size = 1) +
      theme_classic() +
      ggtitle("R2Y Permutation test")
    ggplotly(p)
  })
  
  output$consOPLS_contrib <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    contributions <- consOPLS_res()@blockContribution
    if (is.null(dim(contributions))) {
      contrib_df <- data.frame(dataset = names(contributions),
                               value   = as.numeric(contributions))
      p <- ggplot(contrib_df, aes(x = dataset, y = value, fill = dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Contributions par bloc (vecteur)", x = "Bloc", y = "Contribution") +
        theme_minimal()
    } else {
      df_mat <- as.data.frame(contributions)
      df_mat$dataset <- rownames(contributions)
      contrib_df <- melt(df_mat, id.vars = "dataset", variable.name = "comp", value.name = "value")
      p <- ggplot(contrib_df, aes(x = comp, y = value, fill = dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Contributions par bloc (matrice)", x = "Composante", y = "Contribution") +
        theme_minimal()
    }
    ggplotly(p)
  })
  
  output$consOPLS_scores <- renderPlotly({
    req(consOPLS_res(), metadataData())
    if (is.null(consOPLS_res())) return(plotly_empty())
    scores <- as.data.frame(consOPLS_res()@scores)
    scores$sample_name <- rownames(scores)
    md <- metadataData()
    merged_scores <- merge(scores, md, by = "sample_name", all.x = TRUE)
    xlab <- "Predictive"
    ylab <- "Orthogonal"
    if ("explainedVar" %in% slotNames(consOPLS_res())) {
      ev <- consOPLS_res()@explainedVar
      xlab <- sprintf("Predictive (%.1f%%)", ev["predictive"] * 100)
      ylab <- sprintf("Orthogonal (%.1f%%)", ev["orthogonal"] * 100)
    }
    p <- ggplot(merged_scores, aes(x = p_1, y = o_1,
                                   color = .data[[input$groupVarCons]],
                                   text = sample_name)) +
      geom_point(size = 4) +
      labs(title = "Scores: Predictive vs Orthogonal", x = xlab, y = ylab) +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })
  
  output$consOPLS_loadings <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    loadings_proteo <- consOPLS_res()@loadings$proteo
    loadings_metabo <- consOPLS_res()@loadings$metabo
    loadings <- rbind.data.frame(loadings_proteo, loadings_metabo)
    loadings$dataset <- c(rep("proteo", nrow(loadings_proteo)),
                          rep("metabo", nrow(loadings_metabo)))
    loadings$variable <- rownames(loadings)
    p <- ggplot(loadings, aes(x = p_1, y = o_1, color = dataset, label = variable)) +
      geom_point(size = 2) +
      geom_text(aes(label = variable)) +
      labs(title = "Loadings: Predictive vs Orthogonal", x = "Predictive", y = "Orthogonal") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })
  
  output$consOPLS_VIP <- renderPlotly({
    req(consOPLS_res())
    if (is.null(consOPLS_res())) return(plotly_empty())
    VIP_df <- data.frame(
      VIP = c(consOPLS_res()@VIP$proteo$p, consOPLS_res()@VIP$metabo$p),
      variable = c(rownames(consOPLS_res()@VIP$proteo), rownames(consOPLS_res()@VIP$metabo)),
      dataset  = c(rep("proteo", nrow(consOPLS_res()@VIP$proteo)),
                   rep("metabo", nrow(consOPLS_res()@VIP$metabo)))
    )
    p <- ggplot(VIP_df, aes(x = variable, y = VIP, fill = dataset)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "VIP", x = "Variable", y = "VIP") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggplotly(p)
  })
  
  output$consOPLS_summary <- renderPrint({
    req(consOPLS_res())
    if (!is.null(consOPLS_res())) {
      capture.output(print(consOPLS_res()))
    }
  })
  
  # Download Consensus OPLS
  output$download_ConsensusOPLS <- downloadHandler(
    filename = function() { paste("ConsensusOPLS_", Sys.Date(), ".xlsx", sep = "") },
    content = function(file) {
      wb <- createWorkbook()
      res <- consOPLS_res()
      
      perm_stats <- as.data.frame(res@permStats)
      addWorksheet(wb, "Permutations")
      writeData(wb, "Permutations", perm_stats)
      
      q2 <- data.frame(Q2 = res@Q2)
      addWorksheet(wb, "Q2")
      writeData(wb, "Q2", q2)
      
      dq2 <- data.frame(DQ2 = res@DQ2)
      addWorksheet(wb, "DQ2")
      writeData(wb, "DQ2", dq2)
      
      r2y <- data.frame(R2Y = res@R2Y)
      addWorksheet(wb, "R2Y")
      writeData(wb, "R2Y", r2y)
      
      contrib <- as.data.frame(res@blockContribution)
      addWorksheet(wb, "Contributions")
      writeData(wb, "Contributions", contrib)
      
      scores <- as.data.frame(res@scores)
      scores$sample_name <- rownames(scores)
      addWorksheet(wb, "Scores")
      writeData(wb, "Scores", scores)
      
      loadings <- rbind.data.frame(res@loadings$proteo, res@loadings$metabo)
      loadings$dataset <- c(rep("proteo", nrow(res@loadings$proteo)),
                            rep("metabo", nrow(res@loadings$metabo)))
      loadings$variable <- rownames(loadings)
      addWorksheet(wb, "Loadings")
      writeData(wb, "Loadings", loadings)
      
      VIP_df <- data.frame(
        VIP = c(res@VIP$proteo$p, res@VIP$metabo$p),
        variable = c(rownames(res@VIP$proteo), rownames(res@VIP$metabo)),
        dataset  = c(rep("proteo", nrow(res@VIP$proteo)),
                     rep("metabo", nrow(res@VIP$metabo)))
      )
      addWorksheet(wb, "VIP")
      writeData(wb, "VIP", VIP_df)
      
      res_text <- capture.output(print(res))
      addWorksheet(wb, "Résumé")
      writeData(wb, "Résumé", data.frame(Text = res_text))
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  output$consOPLS_VIP_Loadings <- renderPlotly({
    req(consOPLS_res())
    res <- consOPLS_res()
    
    loadings_proteo <- as.data.frame(res@loadings$proteo)
    loadings_metabo <- as.data.frame(res@loadings$metabo)
    loadings_proteo$variable <- rownames(loadings_proteo)
    loadings_metabo$variable <- rownames(loadings_metabo)
    loadings_proteo$dataset <- "proteo"
    loadings_metabo$dataset <- "metabo"
    loadings_all <- rbind(loadings_proteo, loadings_metabo)
    
    # Crée la table des VIP
    VIP_proteo <- data.frame(VIP = res@VIP$proteo$p,
                             variable = rownames(res@VIP$proteo),
                             dataset = "proteo")
    VIP_metabo <- data.frame(VIP = res@VIP$metabo$p,
                             variable = rownames(res@VIP$metabo),
                             dataset = "metabo")
    VIP_all <- rbind(VIP_proteo, VIP_metabo)
    
    # Fusionne loadings et VIP
    loadings_VIP <- merge(loadings_all, VIP_all, by = c("variable", "dataset"))
    
    # Ajoute un label si le VIP dépasse un seuil (ici 1)
    loadings_VIP$label <- ifelse(loadings_VIP$VIP > 1, loadings_VIP$variable, NA)
    
    #plot
    p <- ggplot(loadings_VIP, aes(x = p_1, y = VIP, color = dataset, label = label)) +
      geom_point(size = 2) +
      labs(x = "Predictive",
           y = "VIP",
           title = "VIP vs loadings on predictive latent variable") +
      geom_text_repel(size = 3, max.overlaps = 50, segment.size = 0.1) +
      theme_light()
    
    ggplotly(p, tooltip = c("label", "x", "y", "color"))
  })
  
}

shinyApp(ui, server)