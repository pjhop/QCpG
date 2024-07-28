## shiny app ----------------------------------------------------------------

run_app <- function(qcmetrics, samplesheet = NULL, background = NULL, show_all_relationships = TRUE, variables = "all") {

    ## Load and preprocess data ---------------------------
    detp_threshold <- qcmetrics@detp_threshold
    beadnr_threshold <- qcmetrics@beadnr_threshold
    array <- qcmetrics@metadata$array

    # Preprocess data
    data <- prepareData(qcmetrics, samplesheet = samplesheet, variables = "all")
    samplesheet <- data$samplesheet
    outlier_variables <- data$outlier_variables
    variables <- data$variables ## Update: .y variables
    data <- data$data

    # Check how many samples have missing sex info
    samples_missing_sex <- samplesheet %>% dplyr::filter(is.na(Sex)) %$% Sample_Name

    # Check how many samples did not have genotype data available
    if(!is.null(qcmetrics@ibs$ibs_geno_identical)) {
      samples_missing_geno <- samplesheet %>%
        dplyr::filter(!(Sample_Name %in% qcmetrics@ibs$ibs_geno_identical$Sample_Name)) %$% Sample_Name
    }

    # Tranformed PCs df -> for defining PC outliers
    if(!is.null(qcmetrics@PCs)) {
      PCs_melt <- qcmetrics@PCs %>% tidyr::gather(key = "PC", value = "value", dplyr::starts_with("PC"))
    }

    controlPCs_melt <- qcmetrics@controlPCs %>% tidyr::gather(key = "PC_control", value = "value", dplyr::ends_with("_control"))

    ## Zhou Masking
    zhou_choices <- c("MASK_general_AFR","MASK_snp5_AFR","MASK_general_EAS", "MASK_snp5_EAS",
     "MASK_general_EUR","MASK_snp5_EUR", "MASK_general_SAS","MASK_snp5_SAS","MASK_general_AMR",
     "MASK_snp5_AMR","MASK_general_GWD", "MASK_snp5_GWD",
     "MASK_general_YRI","MASK_snp5_YRI","MASK_general_TSI",
     "MASK_snp5_TSI","MASK_general_IBS","MASK_snp5_IBS",
     "MASK_general_CHS", "MASK_snp5_CHS", "MASK_general_PUR",
     "MASK_snp5_PUR", "MASK_general_JPT", "MASK_snp5_JPT",
     "MASK_general_GIH", "MASK_snp5_GIH","MASK_general_CH_B",
     "MASK_snp5_CH_B", "MASK_general_STU","MASK_snp5_STU",
     "MASK_general_ITU","MASK_snp5_ITU","MASK_general_LWK",
     "MASK_snp5_LWK","MASK_general_KHV","MASK_snp5_KHV",
     "MASK_general_FIN","MASK_snp5_FIN", "MASK_general_ESN",
     "MASK_snp5_ESN", "MASK_general_CEU", "MASK_snp5_CEU",
     "MASK_general_PJL","MASK_snp5_PJL","MASK_general_AC_B",
     "MASK_snp5_AC_B", "MASK_general_CLM", "MASK_snp5_CLM",
     "MASK_general_CDX","MASK_snp5_CDX", "MASK_general_GBR",
     "MASK_snp5_GBR","MASK_general_BE_B", "MASK_snp5_BE_B",
     "MASK_general_PEL","MASK_snp5_PEL", "MASK_general_MSL",
     "MASK_snp5_MSL", "MASK_general_MXL", "MASK_snp5_MXL",
     "MASK_general_ASW","MASK_snp5_ASW")

     ## Thresholds
     # if(qcmetrics@metadata$array == "450k") {
     #     thresholds <- list(MU = 1500, RG_ratio = 0.5, GR_ratio = 0.5, OP = 11.75, HC = 12, bscon = 80,
     #                            detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2)
     # } else if (qcmetrics@metadata$array == "EPIC") {
     #     thresholds <- list(MU = 2000, RG_ratio = 0.4, GR_ratio = 0.4, OP = 12, HC = 12.75, bscon = 80,
     #                            detP = 0.05, beadNr = 0.05, IBS_mean = 1.9, IBS_var = 0.1, XY_diff = -2)
     # } else {
     #     stop("Array type not known!")
     # }
     thresholds <- qcmetrics@sample_qc_outliers$thresholds
     probe_thresholds <- qcmetrics@probe_qc_outliers$thresholds
     relatedness_thresholds <- qcmetrics@relatedness_thresholds

     if(!show_all_relationships && !is.null(qcmetrics@ibs$ibs_dnam)) {
       qcmetrics@ibs <- .summarize_relationships(
         ibs = qcmetrics@ibs,
         relatedness_thresholds = qcmetrics@relatedness_thresholds
       )
       qcmetrics@ibs$save_all_relationships <- FALSE
     }

    ## Create background plots (if background is specified)
    if(!is.null(background)) {
      if(is(background, "QCmetrics")) {
        background <- makeBackground(background)
      }
      data_bg <- background@Sample_Metrics

      # Check if any background data is also in present in foreground data
      overlap <- sum(data_bg$Sample_Name %in% data$Sample_Name)
      if(overlap > 0) {
        print(sprintf("%s samples in the foreground data are present in the background data. They will be removed from the background.", overlap))
        data_bg <- data_bg %>% dplyr::filter(!(Sample_Name %in% data$Sample_Name))
      }
      if(nrow(data_bg) == 0) {
        # If there is total overlap between foreground and background, the background option is turned off
        background <- NULL
      } else {
        # Create hexabin backgrounds
        bg_xy <- ggplot(data_bg, aes(x = xMed, y = yMed)) + geom_hex(alpha = 0.6, show.legend=FALSE) + scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
        bg_mu <- ggplot(data_bg, aes(x = Unmethylated, y = Methylated)) + geom_hex(alpha = 0.6, show.legend=FALSE) + scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
        bg_rg <- ggplot(data_bg, aes(x = Red, y = Green)) + geom_hex(alpha = 0.6, show.legend=FALSE) + scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
        bg_op <- ggplot(data_bg, aes(x = OP_x, y = OP_y)) + geom_hex(alpha = 0.6, show.legend=FALSE) + scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
        bg_hc <- ggplot(data_bg, aes(x = HC_x, y = HC_y)) + geom_hex(alpha = 0.6, show.legend=FALSE) + scale_fill_gradient2(low ="white", mid = "gray95", high = "#0072B5")
      }
    } else {
      background <- NULL
    }

  ### App ---------------------------------------------------------
  shinyApp(

    ui <- fluidPage(
      navbarPage(HTML(paste0("QCpG",'<sup>', 'v', packageVersion("QCpG"))),

         ### Sample QC tab -------
         tabPanel("Sample QC",

          ## Sidebarpanel -------
          sidebarPanel(width=3,
                        h5(HTML(sprintf("<b>Name:</b> %s", qcmetrics@metadata$experiment_name))),
                        h5(HTML(sprintf("<b>Experimenter:</b> %s", qcmetrics@metadata$experimenter))),
                        h5(HTML(sprintf("<b>Sample size:</b> %s", nrow(samplesheet)))),
                        hr(),
                       selectInput(inputId = "selected_colorscheme",
                                   label = "Color scheme",
                                   choices = c("viridis", "Okabe-Ito", "ggplot"),
                                   selected = "viridis", multiple = FALSE, selectize = TRUE, width = "98%"),
                       selectInput(inputId = "scatterColorByColumn", label = "Color datapoints by", choices = c(variables, "predictedSex"),
                                    selected = "None", multiple = FALSE, selectize = TRUE, width = "98%"),

                       conditionalPanel(
                         condition = "input.define_outliers == 1",
                             selectInput(inputId = "colorbyQCmetric", label = "Color outliers per metric",
                                         choices = c(outlier_variables, paste0(outlier_variables[outlier_variables != "None"], "_unique")),
                                         selected = "None", multiple = FALSE, selectize = TRUE, width = "98%")
                       ),
                       checkboxInput(inputId = "webGL", label = "Use WebGL?", value = FALSE),
                        uiOutput("reset_selected_samples"),
                       conditionalPanel(
                         condition="input.tabselected == 0 && input.pca == 0.1",
                         selectInput(inputId = "scattercontrolPCcolumns", label = "PCs to plot", choices = paste0("PC", 1:10, "_control"),
                                     selected = c("PC1_control", "PC2_control"), multiple = TRUE, selectize = TRUE, width = "98%"),
                         checkboxInput(inputId = "pairsplot2", label = "Make pairs plot?", value = FALSE)
                       ),
                        conditionalPanel(
                          condition = "input.tabselected == 0 && input.pca == 0.2",
                          selectInput(inputId = "scatterPCcolumns", label = "PCs to plot",
                                      choices = paste0("PC", 1:10),selected = c("PC1", "PC2"),
                                      multiple = TRUE, selectize = TRUE, width = "98%"),
                          checkboxInput(inputId = "pairsplot", label = "Make pairs plot?", value = FALSE)
                        ),
                       uiOutput("show_background"),
                       hr(),
                       h4("Outliers"),
                       checkboxInput(inputId = "define_outliers", label = "Define Outliers?", value = FALSE),
                        conditionalPanel(
                          condition = "input.define_outliers == 1",

                          checkboxInput(inputId = "show_outliers", label = "Show all Outliers?", value = TRUE),

                          ## Thresholds
                          h4("Thresholds"),
                          checkboxInput(inputId = "show_all_thresholds", "Show all thresholds?",
                                        value = 0),
                          conditionalPanel(
                            condition="input.tabselected == 0 && input.pca == 0.1 && input.show_all_thresholds == 0",
                            checkboxInput(inputId = "controlPCA_outliers_tab",
                                          "Define control PCA outliers?",
                                          value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) TRUE else FALSE)
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 0 && input.pca == 0.1 && input.controlPCA_outliers_tab == 1 && input.show_all_thresholds == 0",
                            selectInput(inputId = "controlPCs_used_tab",
                                        label = "control PCs used to define outliers",
                                        choices = paste0("PC", 1:10, "_control"),
                                        selected = if(!is.null(thresholds$controlPCs_used) && !is.na(thresholds$controlPCs_used)) thresholds$controlPCs_used else "PC1_control", multiple = TRUE, selectize = TRUE, width = "98%"),
                            numericInput(inputId = "controlPCA_sd_tab", label = "Number of SDs",
                                         min = 0,
                                         max = 10,
                                         value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) thresholds$controlPCA_sd else 3, width = "98%")
                          ),
                          conditionalPanel(
                            condition="input.tabselected == 0 && input.pca == 0.2 && input.show_all_thresholds == 0",
                            checkboxInput(inputId = "PCA_outliers_tab",
                                          "Define PCA outliers?",
                                          value = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) TRUE else FALSE)
                          ),
                          conditionalPanel(
                            condition="input.tabselected == 0 && input.pca == 0.2 && input.PCA_outliers_tab == 1 && input.show_all_thresholds == 0",
                            selectInput(inputId = "PCs_used_tab",
                                        label = "PCs used to define outliers",
                                        choices = paste0("PC", 1:10),
                                        selected = if(!is.null(thresholds$PCs_used) && !is.na(thresholds$PCs_used)) thresholds$PCs_used else "PC1",
                                        multiple = TRUE,
                                        selectize = TRUE,
                                        width = "98%"),
                            numericInput(
                              inputId = "PCA_sd_tab",
                              label = "Number of SDs",
                              min = 0,
                              max = 10,
                              value = if (!is.null(thresholds$PCA_sd) &&
                                          !is.na(thresholds$PCA_sd))
                                thresholds$PCA_sd
                              else
                                3,
                              width = "98%"
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 1 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "XY_diff_tab",
                              "XY Diff",
                              value = thresholds[["XY_diff"]],
                              min = -5,
                              max = 0
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 2 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "MU_threshold_tab",
                              "MU",
                              value = thresholds[["MU"]],
                              min = 0,
                              max = 10000,
                              step = 500
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 3 && input.show_all_thresholds == 0",
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "RG_ratio_threshold_tab",
                                "R/G ratio",
                                value = thresholds[["RG_ratio"]],
                                min = 0,
                                max = 3,
                                step = 0.1
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "GR_ratio_threshold_tab",
                                "G/R ratio",
                                value = thresholds[["GR_ratio"]],
                                min = 0,
                                max = 3,
                                step = 0.1
                              )
                            ))
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 4 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "OP_threshold_tab",
                              "OP",
                              value = thresholds[["OP"]],
                              min = 0,
                              max = 100,
                              step = 0.5
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 5 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "HC_threshold_tab",
                              "Hyb",
                              value = thresholds[["HC"]],
                              min = 0,
                              max = 100,
                              step = 0.5
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 6 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "bscon_threshold_tab",
                              "bscon",
                              value = thresholds[["bscon"]],
                              min = 0,
                              max = 100,
                              step = 1
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 7 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "detP_threshold_tab",
                              "detP",
                              thresholds[["detP"]],
                              min = 0,
                              max = 1,
                              step = 0.01
                            )
                          ),
                          conditionalPanel(
                            condition = "input.tabselected == 8 && input.show_all_thresholds == 0",
                            numericInput(
                              inputId = "beadNr_threshold_tab",
                              "beadNr",
                              thresholds[["beadNr"]],
                              min = 0,
                              max = 1,
                              step = 0.01
                            )
                          ),
                          conditionalPanel(condition = "input.tabselected == 9 && input.show_all_thresholds == 0",
                                           fluidRow(column(
                                             6,
                                             numericInput(
                                               inputId = "IBS_mean_tab",
                                               "IBS mean",
                                               thresholds[["IBS_mean"]],
                                               min = 0,
                                               max = 2,
                                               step = 0.005
                                             )
                                           ),
                                           column(
                                             6,
                                             numericInput(
                                               inputId = "IBS_var_tab",
                                               "IBS var",
                                               thresholds[["IBS_var"]],
                                               min = 0,
                                               max = 1,
                                               step = 0.005
                                             )
                                           ))),
                          conditionalPanel(condition = "input.tabselected == 11 && input.show_all_thresholds == 0",
                                           fluidRow(column(
                                             10, uiOutput("custom_outliers_sidepanel_tab")
                                           ))),
                          conditionalPanel(
                            condition = "input.show_all_thresholds == 1",
                            checkboxInput(
                              inputId = "controlPCA_outliers_main",
                              "Define control PCA outliers?",
                              value = if (!is.null(thresholds$controlPCA_sd) &&
                                          !is.na(thresholds$controlPCA_sd))
                                TRUE
                              else
                                FALSE
                            ),
                            conditionalPanel(
                              condition = "input.controlPCA_outliers_main == 1",
                              selectInput(
                                inputId = "controlPCs_used_main",
                                label = "control PCs used to define outliers",
                                choices = paste0("PC", 1:10, "_control"),
                                selected = if (!is.null(thresholds$controlPCs_used) &&
                                               !is.na(thresholds$controlPCs_used))
                                  thresholds$controlPCs_used
                                else
                                  "PC1_control",
                                multiple = TRUE,
                                selectize = TRUE,
                                width = "98%"
                              ),
                              numericInput(
                                inputId = "controlPCA_sd_main",
                                label = "Number of SDs",
                                min = 0,
                                max = 10,
                                value = if (!is.null(thresholds$controlPCA_sd) &&
                                            !is.na(thresholds$controlPCA_sd))
                                  thresholds$controlPCA_sd
                                else
                                  3,
                                width = "98%"
                              )
                            ),
                            checkboxInput(
                              inputId = "PCA_outliers_main",
                              "Define PCA outliers?",
                              value = if (!is.null(thresholds$PCA_sd) &&
                                          !is.na(thresholds$PCA_sd))
                                TRUE
                              else
                                FALSE
                            ),
                            conditionalPanel(
                              condition = "input.PCA_outliers_main == 1",
                              selectInput(
                                inputId = "PCs_used_main",
                                label = "PCs used to define outliers",
                                choices = paste0("PC", 1:10),
                                selected = if (!is.null(thresholds$PCs_used) &&
                                               !is.na(thresholds$PCs_used))
                                  thresholds$PCs_used
                                else
                                  "PC1",
                                multiple = TRUE,
                                selectize = TRUE,
                                width = "98%"
                              ),
                              numericInput(
                                inputId = "PCA_sd_main",
                                label = "Number of SDs",
                                min = 0,
                                max = 10,
                                value = if (!is.null(thresholds$PCA_sd) &&
                                            !is.na(thresholds$PCA_sd))
                                  thresholds$PCA_sd
                                else
                                  3,
                                width = "98%"
                              )
                            ),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "XY_diff_main",
                                "XY Diff",
                                value = thresholds[["XY_diff"]],
                                min = -5,
                                max = 0
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "MU_threshold_main",
                                "MU",
                                value = thresholds[["MU"]],
                                min = 0,
                                max = 10000,
                                step = 500
                              )
                            )),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "RG_ratio_threshold_main",
                                "R/G ratio",
                                value = thresholds[["RG_ratio"]],
                                min = 0,
                                max = 3,
                                step = 0.1
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "GR_ratio_threshold_main",
                                "G/R ratio",
                                value = thresholds[["GR_ratio"]],
                                min = 0,
                                max = 3,
                                step = 0.1
                              )
                            )),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "OP_threshold_main",
                                "OP",
                                value = thresholds[["OP"]],
                                min = 0,
                                max = 100,
                                step = 0.5
                              )
                            )),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "HC_threshold_main",
                                "Hyb",
                                thresholds[["HC"]],
                                min = 0,
                                max = 100,
                                step = 0.5
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "bscon_threshold_main",
                                "bscon",
                                thresholds[["bscon"]],
                                min = 0,
                                max = 100,
                                step = 1
                              )
                            )),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "detP_threshold_main",
                                "detP",
                                thresholds[["detP"]],
                                min = 0,
                                max = 1,
                                step = 0.01
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "beadNr_threshold_main",
                                "beadNr",
                                thresholds[["beadNr"]],
                                min = 0,
                                max = 1,
                                step = 0.01
                              )
                            )),
                            fluidRow(column(
                              6,
                              numericInput(
                                inputId = "IBS_mean_main",
                                "IBS mean",
                                thresholds[["IBS_mean"]],
                                min = 0,
                                max = 2,
                                step = 0.005
                              )
                            ),
                            column(
                              6,
                              numericInput(
                                inputId = "IBS_var_main",
                                "IBS var",
                                thresholds[["IBS_var"]],
                                min = 0,
                                max = 1,
                                step = 0.005
                              )
                            )),
                            fluidRow(column(
                              10, uiOutput("custom_outliers_sidepanel_main")
                            )),
                          ),

                          br(),
                          actionButton("reset_defaults", "Reset thresholds to default"),
                          downloadButton("newqcmetrics", "Download QCmetrics file"),
                        )#,
                      #  conditionalPanel(
                       #   condition = "input.define_outliers == 1",
                       #   br(),
                         # downloadButton("downloadData", "Download outlier table")
                      #  )

          ),

          ## Main Panel -------------
          mainPanel(
            tabsetPanel(
              id = "tabselected",
              tabPanel("PCA", fluid = TRUE, value = 0,
                tabsetPanel(
                  id = "pca",
                  tabPanel("Control Probes", value = 0.1,
                           textOutput(outputId = "nr.outliers.total.controlpca"),
                           textOutput(outputId = "nr.outliers.controlpca"),
                           plotlyOutput(outputId = "control_PCA"),
                           plotOutput(outputId = "Heatmap_controlPCA"),
                           plotOutput(outputId = "PairsControlPCs")
                  ),
                  tabPanel("Array-Wide", value = 0.2,
                       textOutput(outputId = "nr.outliers.total.pca"),
                       textOutput(outputId = "nr.outliers.pca"),
                       textOutput(outputId = "nr.probes.used.pca"),
                       plotlyOutput(outputId = "PCA"),
                       plotOutput(outputId = "Heatmap"),
                       plotOutput(outputId = "Pairs")
                 )
                 )),
              tabPanel("Sex", value = 1, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.sex"),
                       textOutput(outputId = "sex.outliers"),
                       uiOutput("samples.missing.sex.text"),
                       plotlyOutput(outputId = "xyplot"),
                       plotlyOutput(outputId = "xydiff")

              ),
              tabPanel("MU", value = 2, fluid = TRUE,
                      textOutput(outputId = "nr.outliers.total.mu"),
                      textOutput(outputId = "MU.outliers"),
                      plotlyOutput(outputId = "MU")
              ),
              tabPanel("RG", value = 3, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.rg"),
                       textOutput(outputId = "RG.ratio.outliers"),
                       textOutput(outputId = "GR.ratio.outliers"),
                       plotlyOutput(outputId = "RG")
              ),
              tabPanel("OP", value = 4, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.op"),
                       textOutput(outputId = "OP.outliers"),
                       plotlyOutput(outputId = "OP")

              ),
              tabPanel("HC", value = 5, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.hc"),
                       textOutput(outputId = "HC.outliers"),
                       plotlyOutput(outputId = "HC")
              ),
              tabPanel("bscon",value = 6, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.bscon"),
                       textOutput(outputId = "bscon.outliers"),
                       plotlyOutput(outputId = "bscon")
              ),
               tabPanel("detP", value = 7, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.detp"),
                       textOutput(outputId = "detp.outliers"),
                       plotlyOutput(outputId = "detP")
              ),
               tabPanel("beadNr", value = 8, fluid = TRUE,
                       textOutput(outputId = "nr.outliers.total.beadnr"),
                       textOutput(outputId = "beadnr.outliers"),
                       plotlyOutput(outputId = "beadNr")
              ),
              tabPanel("IBS", value = 9, fluid = TRUE,
                       h5("Match between sample's (inferred) DNAm SNPs and Genotype Array/WGS SNPs"),
                       textOutput(outputId = "nr.outliers.total.ibs.intra"),
                       textOutput(outputId = "ibs.outliers"),
                       textOutput(outputId = "samples.missing.geno"),
                       textOutput(outputId = "nr.called.probes.geno"),
                       plotlyOutput(outputId = "intra_ibs")
                  ),
               tabPanel("PMS", value = 10, fluid=TRUE,
                        tabsetPanel(
                          id = "predpheno",
                          tabPanel("CTF",
                                   textOutput(outputId = "nr.outliers.total.ctf"),
                                   plotlyOutput(outputId = "CTF")

                          ),
                          tabPanel("Age",
                                   textOutput(outputId = "nr.outliers.total.age"),
                                   textOutput(outputId = "age.text"),
                                   textOutput(outputId = "age.text.missing"),
                                   selectInput(inputId = "select_age_predictor", label = "Which age predictor to use?",
                                               choices = c("Horvath", "Zhang_ElasticNet", "Zhang_BLUP")),
                                   plotlyOutput(outputId = "Age")

                          ),
                          tabPanel("Smoking",
                                   textOutput(outputId = "nr.outliers.total.smoking"),
                                   textOutput(outputId = "smoking.text"),
                                   textOutput(outputId = "smoking.text.missing"),
                                   plotlyOutput(outputId = "Smoking")
                          ),
                          tabPanel("Alcohol",
                                   textOutput(outputId = "nr.outliers.total.alcohol"),
                                   textOutput(outputId = "alcohol.text"),
                                   textOutput(outputId = "alcohol.text.missing"),
                                   plotlyOutput(outputId = "Alcohol")
                          ),
                          tabPanel("BMI",
                                   textOutput(outputId = "nr.outliers.total.bmi"),
                                   textOutput(outputId = "bmi.text"),
                                   textOutput(outputId = "bmi.text.missing"),
                                   plotlyOutput(outputId = "BMI")
                          )
                  )
                ),
               tabPanel("Custom", value = 11, fluid=TRUE,
                       textOutput(outputId = "nr.outliers.total.custom"),
                       textOutput(outputId = "custom.outliers"),
                       selectInput(inputId = "custom_columns", label = "Choose variables", choices = colnames(samplesheet)[!colnames(samplesheet) == "Sample_Name"],
                              selected = NULL, multiple = TRUE, selectize = TRUE, width = "50%"),
                      fluidRow(column(1,actionButton("add_custom_columns", "Add")),
                               column(1,actionButton("reset_custom_columns", "Reset"))),
                      textOutput(outputId = "test"),
                      h5("Note: columns must be logical (TRUE/FALSE)")

              ),
              tabPanel("OutlierTable", fluid = TRUE,
                       h3("Overview Sample QC outliers"),
                       uiOutput("download_overview_button"),
                       DT::dataTableOutput(outputId = "overview_outliers"),
                       h3("Sample QC outliers"),
                       DT::dataTableOutput(outputId = "outliers")

              )
            )
          )),

    ### Site QC tab -----------------

    tabPanel("Site QC",
        ## Sidebar panel
        sidebarPanel(width=3,
                        checkboxInput(inputId = "define_outliers_probes", label = "Define Outliers?", value = FALSE),

                        conditionalPanel(
                          condition = "input.define_outliers_probes == 1",
                          numericInput(inputId = "detP_probes", "detP percentage cut-off", probe_thresholds$detP, min = 0, max = 1, step = 0.01),
                          numericInput(inputId = "beadNr_probes", "beadNr percentage cut-off", probe_thresholds$beadNr, min = 0, max = 1, step = 0.01),
                          checkboxInput(inputId = "show_outliers_probes", label = "Show Outliers?", value = FALSE),
                          br(),
                          fluidRow(column(6, actionButton("reset_defaults_probes", "Reset thresholds to default")),
                                   column(6, downloadButton("newqcmetrics_probes", "Download QCmetrics file")))
          )),

        ## Main panel -----------
        mainPanel(
            tabsetPanel(
              id = "tabselected",
              tabPanel("detP", value = 1, fluid = TRUE,

                       textOutput(outputId = "nr.outliers.probes"),
                       textOutput(outputId = "nr.outliers.probes.detp"),
                       plotOutput(outputId = "detP_probes"),
                       p("Note that probes with proportion < 0.01 are not displayed, in order to decrease loading times."),
                       textOutput(outputId = "note_probeqc")
              ),
              tabPanel("beadNr", fluid = TRUE,
                       textOutput(outputId = "nr.outliers.probes2"),
                       textOutput(outputId = "nr.outliers.probes.beadnr"),
                       plotOutput(outputId = "beadNr_probes"),
                       p("Note that probes with proportion < 0.01 are not displayed, in order to decrease loading times."),
                       textOutput(outputId = "note_probeqc2")

              ),
              tabPanel("Masking", fluid = TRUE,
                       textOutput(outputId = "nr.outliers.probes3"),
                       textOutput(outputId = "nr_probes_zhou"),
                       textOutput(outputId = "zhou_population_used"),
                       selectInput(inputId = "zhou_column", "Choose a population*", choices = zhou_choices,
                                   selected = "None", multiple = FALSE, selectize = TRUE, width = "30%"),
                       actionButton("download_zhou", "Download"),
                       textOutput(outputId = "note_download_zhou"),
                       br(),
                       br(),
                       br(),
                       textOutput(outputId = "note_masking")

              )
            )
          )),
    tabPanel("Relatedness",
             textOutput(outputId = "ibs.summarized.text"),
             h5("DNAm-inferred SNPs x DNAm-inferred SNPs"),
             textOutput(outputId = "nr.called.probes.dnam"),
             fluidRow(column(9,plotlyOutput(outputId = "inter_ibs_dnam")),
                      column(3,
                             sliderInput(inputId = "ibs_relatedness_mean", "IBS Mean", relatedness_thresholds$IBS_mean,
                                          min = if(qcmetrics@ibs$save_all_relationships) 0 else relatedness_thresholds$IBS_mean,
                                          max = 2,
                                          step = 0.05),
                             sliderInput(inputId = "ibs_relatedness_var", "IBS Var", relatedness_thresholds$IBS_var,
                                          min = 0,
                                          max = if(qcmetrics@ibs$save_all_relationships) 1 else relatedness_thresholds$IBS_var,
                                          step = 0.05),
                             checkboxInput(inputId = "show_pairs", label = "Show related pairs?", value = FALSE),
                             downloadButton("download_related_samples", "Download related samples"))
             ),
             br(),

             h5("DNAm-inferred SNPs x Genotyped SNPs"),
             verbatimTextOutput(outputId = "nr.outliers.total.inter.ibs"),
             textOutput(outputId = "samples.missing.geno2"),
             fluidRow(column(9,plotlyOutput(outputId = "inter_ibs")))
             )
      )
    ),


    #### Server -------------------------------------------------------------------------

    server <- function(input, output, session) {

      melt_PCs <- reactive({
        dat <- PCs_melt %>% dplyr::mutate(PC = factor(PC, levels = paste0("PC", 1:10))) %>%
        dplyr::group_by(PC) %>%
          dplyr::mutate(mean = mean(value), sd = sqrt(var(value)),
                 sd_min = mean - as.numeric(values$PCA_sd) * sd, sd_plus = mean + as.numeric(values$PCA_sd) * sd,
                 outlier = ifelse(value > sd_plus | value < sd_min, TRUE, FALSE))

        PCs_melt_sel <- dat[dat$PC %in% values$PCs_used,]
        PCs_melt_sel
      })

      melt_controlPCs <- reactive({
        dat <- controlPCs_melt %>% dplyr::mutate(PC_control = factor(PC_control, levels = paste0("PC", 1:10, "_control"))) %>%
          dplyr::group_by(PC_control) %>%
          dplyr::mutate(mean = mean(value), sd = sqrt(var(value)),
                        sd_min = mean - as.numeric(values$controlPCA_sd) * sd, sd_plus = mean + as.numeric(values$controlPCA_sd) * sd,
                        outlier = ifelse(value > sd_plus | value < sd_min, TRUE, FALSE))

        PCs_melt_sel <- dat[dat$PC_control %in% values$controlPCs_used,]
        PCs_melt_sel
      })

      ### Reactive Values (outliers, selected samples, thresholds)
      values <- reactiveValues(outlier_vars = if(all(!is.na(thresholds[["custom_outliers"]]))) c(outlier_variables, thresholds[["custom_outliers"]]) else outlier_variables ,
                               selected_samples = if(!is.na(qcmetrics@selected_samples)) qcmetrics@selected_samples else c(),
                               webGL = FALSE,
                               selected_colorscheme = "viridis",
                               colorscheme = scale_color_viridis_d(),
                               XY_diff = thresholds[["XY_diff"]],
                               MU_threshold = thresholds[["MU"]],
                               RG_ratio_threshold = thresholds[["RG_ratio"]],
                               GR_ratio_threshold = thresholds[["GR_ratio"]],
                               OP_threshold = thresholds[["OP"]],
                               HC_threshold = thresholds[["HC"]],
                               bscon_threshold = thresholds[["bscon"]],
                               detP_threshold = thresholds[["detP"]],
                               beadNr_threshold = thresholds[["beadNr"]],
                               IBS_mean = thresholds[["IBS_mean"]],
                               IBS_var = thresholds[["IBS_var"]],
                               PCA_outliers = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) TRUE else FALSE,
                               PCs_used = if(!is.null(thresholds$PCs_used) && !is.na(thresholds$PCs_used)) thresholds$PCs_used else NULL,
                               PCA_sd = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) thresholds$PCA_sd else NULL,
                               controlPCA_outliers = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) TRUE else FALSE,
                               controlPCs_used = if(!is.null(thresholds$controlPCs_used) && !is.na(thresholds$controlPCs_used)) thresholds$controlPCs_used else NULL,
                               controlPCA_sd = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) thresholds$controlPCA_sd else NULL
                               )
      ### update colorscheme
      observeEvent(input$selected_colorscheme,{
        values$selected_colorscheme <- input$selected_colorscheme
      })

      ### Update reactiveValues when input changes
      observeEvent(input$XY_diff_main,{
        values$XY_diff <- input$XY_diff_main
      })
      observeEvent(input$XY_diff_tab,{
        values$XY_diff <- input$XY_diff_tab
      })
      observeEvent(input$MU_threshold_main,{
        values$MU_threshold <- input$MU_threshold_main
      })
      observeEvent(input$MU_threshold_tab,{
        values$MU_threshold <- input$MU_threshold_tab
      })
      observeEvent(input$RG_ratio_threshold_main,{
        values$RG_ratio_threshold <- input$RG_ratio_threshold_main
      })
      observeEvent(input$RG_ratio_threshold_tab,{
        values$RG_ratio_threshold <- input$RG_ratio_threshold_tab
      })
      observeEvent(input$GR_ratio_threshold_main,{
        values$GR_ratio_threshold <- input$GR_ratio_threshold_main
      })
      observeEvent(input$GR_ratio_threshold_tab,{
        values$GR_ratio_threshold <- input$GR_ratio_threshold_tab
      })
      observeEvent(input$OP_threshold_main,{
        values$OP_threshold <- input$OP_threshold_main
      })
      observeEvent(input$OP_threshold_tab,{
        values$OP_threshold <- input$OP_threshold_tab
      })
      observeEvent(input$HC_threshold_main,{
        values$HC_threshold <- input$HC_threshold_main
      })
      observeEvent(input$HC_threshold_tab,{
        values$HC_threshold <- input$HC_threshold_tab
      })
      observeEvent(input$bscon_threshold_main,{
        values$bscon_threshold <- input$bscon_threshold_main
      })
      observeEvent(input$bscon_threshold_tab,{
        values$bscon_threshold <- input$bscon_threshold_tab
      })
      observeEvent(input$detP_threshold_main,{
        values$detP_threshold <- input$detP_threshold_main
      })
      observeEvent(input$detP_threshold_tab,{
        values$detP_threshold <- input$detP_threshold_tab
      })
      observeEvent(input$beadNr_threshold_main,{
        values$beadNr_threshold <- input$beadNr_threshold_main
      })
      observeEvent(input$beadNr_threshold_tab,{
        values$beadNr_threshold <- input$beadNr_threshold_tab
      })
      observeEvent(input$IBS_mean_main,{
        values$IBS_mean <- input$IBS_mean_main
      })
      observeEvent(input$IBS_mean_tab,{
        values$IBS_mean <- input$IBS_mean_tab
      })
      observeEvent(input$IBS_var_main,{
        values$IBS_var <- input$IBS_var_main
      })
      observeEvent(input$IBS_var_tab,{
        values$IBS_var <- input$IBS_var_tab
      })
      observeEvent(input$PCA_outliers_main,{
        values$PCA_outliers <- input$PCA_outliers_main
      })
      observeEvent(input$PCA_outliers_tab,{
        values$PCA_outliers <- input$PCA_outliers_tab
      })
      observeEvent(input$PCs_used_main,{
        values$PCs_used <- input$PCs_used_main
      })
      observeEvent(input$PCA_sd_main,{
        values$PCA_sd <- input$PCA_sd_main
      })
      observeEvent(input$PCs_used_tab,{
        values$PCs_used <- input$PCs_used_tab
      })
      observeEvent(input$PCA_sd_tab,{
        values$PCA_sd <- input$PCA_sd_tab
      })
      observeEvent(input$controlPCA_outliers_main,{
        values$controlPCA_outliers <- input$controlPCA_outliers_main
      })
      observeEvent(input$controlPCA_outliers_tab,{
        values$controlPCA_outliers <- input$controlPCA_outliers_tab
      })
      observeEvent(input$controlPCs_used_main,{
        values$controlPCs_used <- input$controlPCs_used_main
      })
      observeEvent(input$controlPCA_sd_main,{
        values$controlPCA_sd <- input$controlPCA_sd_main
      })
      observeEvent(input$controlPCs_used_tab,{
        values$controlPCs_used <- input$controlPCs_used_tab
      })
      observeEvent(input$controlPCA_sd_tab,{
        values$controlPCA_sd <- input$controlPCA_sd_tab
      })

      ## Update relevant inputs
      observeEvent(input$show_all_thresholds, {
        if (input$show_all_thresholds == 1) {
            updateNumericInput(session, "XY_diff_main", value = values$XY_diff)
            updateNumericInput(session, "MU_threshold_main", value = values$MU_threshold)
            updateNumericInput(session, "RG_ratio_threshold_main", value = values$RG_ratio_threshold)
            updateNumericInput(session, "GR_ratio_threshold_main", value = values$GR_ratio_threshold)
            updateNumericInput(session, "OP_threshold_main", value = values$OP_threshold)
            updateNumericInput(session, "HC_threshold_main", value = values$HC_threshold)
            updateNumericInput(session, "bscon_threshold_main", value = values$bscon_threshold)
            updateNumericInput(session, "detP_threshold_main", value = values$detP_threshold)
            updateNumericInput(session, "beadNr_threshold_main", value = values$beadNr_threshold)
            updateNumericInput(session, "IBS_mean_main", value = values$IBS_mean)
            updateNumericInput(session, "IBS_var_main", value = values$IBS_var)
            updateCheckboxInput(session, "PCA_outliers_main", value = values$PCA_outliers)
            updateSelectInput(session, "PCs_used_main", selected = values$PCs_used)
            updateNumericInput(session, "PCA_sd_main", value = values$PCA_sd)
            updateCheckboxInput(session, "controlPCA_outliers_main", value = values$controlPCA_outliers)
            updateSelectInput(session, "controlPCs_used_main", selected = values$controlPCs_used)
            updateNumericInput(session, "controlPCA_sd_main", value = values$controlPCA_sd)
        }

        if (input$show_all_thresholds == 0) {
          updateNumericInput(session, "XY_diff_tab", value = values$XY_diff)
          updateNumericInput(session, "MU_threshold_tab", value = values$MU_threshold)
          updateNumericInput(session, "RG_ratio_threshold_tab", value = values$RG_ratio_threshold)
          updateNumericInput(session, "GR_ratio_threshold_tab", value = values$GR_ratio_threshold)
          updateNumericInput(session, "OP_threshold_tab", value = values$OP_threshold)
          updateNumericInput(session, "HC_threshold_tab", value = values$HC_threshold)
          updateNumericInput(session, "bscon_threshold_tab", value = values$bscon_threshold)
          updateNumericInput(session, "detP_threshold_tab", value = values$detP_threshold)
          updateNumericInput(session, "beadNr_threshold_tab", value = values$beadNr_threshold)
          updateNumericInput(session, "IBS_mean_tab", value = values$IBS_mean)
          updateNumericInput(session, "IBS_var_tab", value = values$IBS_var)
          updateCheckboxInput(session, "PCA_outliers_tab", value = values$PCA_outliers)
          updateSelectInput(session, "PCs_used_tab", selected = values$PCs_used)
          updateNumericInput(session, "PCA_sd_tab", value = values$PCA_sd)
          updateCheckboxInput(session, "controlPCA_outliers_tab", value = values$controlPCA_outliers)
          updateSelectInput(session, "controlPCs_used_tab", selected = values$controlPCs_used)
          updateNumericInput(session, "controlPCA_sd_tab", value = values$controlPCA_sd)
        }
      })

      ### Plotly click/select events

      # Set selected samples on click
      observeEvent(event_data("plotly_click"),{
        d <- event_data("plotly_click")
        values$selected_samples <- unique(c(values$selected_samples,d$key))
      })

      # Set selected samples on selection (box or lasso)
      observeEvent(event_data("plotly_selected"),{
        d <- event_data("plotly_selected")
        values$selected_samples <- unique(c(values$selected_samples,d$key))
      })

      # Reset selected samples
      observeEvent(input$reset_samples, {
        values$selected_samples <- c()
      }, ignoreInit = TRUE)

      # If samples are selected, the show_outliers button is set to 0
      observeEvent(event_data("plotly_click"), {
          updateCheckboxInput(session, inputId = "show_outliers", value = 0)
         })
      observeEvent(event_data("plotly_selected"), {
             updateCheckboxInput(session, inputId = "show_outliers", value = 0)
            })

      # If show_outliers button is clicked, selected samples are resetted
      observe({
        req(values$selected_samples)
        if(input$show_outliers == 1 && length(outliers()) > 0)
          values$selected_samples <- c()
      })
      # observeEvent(input$show_outliers, {
      #   if(input$show_outliers == 1)
      #     values$selected_samples <- c()
      #   }, ignoreInit = FALSE
      # )

      observe({
        x <- input$define_outliers
        updateCheckboxInput(session, "define_outliers_probes", value = x)
      })

      observe({
        y <- input$define_outliers_probes
        updateCheckboxInput(session, "define_outliers", value = y)
      })

      observe({
        if(input$webGL) {
          updateCheckboxInput(session, inputId = "show_background", value = FALSE)
          values$webGL <- TRUE
        } else {
          values$webGL <- FALSE
        }
      })

      ### Add custom columns
      observeEvent(input$add_custom_columns,{
        # Check if selected variables are present in the samplesheet, do not contain missings and the column is logical
        output$test <- renderText({
                validate(
                  need((all(isolate(input$custom_columns) %in% colnames(data))), "Variables should be in samplesheet!"),
                  need(all(purrr::map_lgl(data[,isolate(input$custom_columns)], .f = is.logical)) ,"Variables should be logical (TRUE/FALSE)!"),
                  need(sum(is.na(data[,isolate(input$custom_columns)])) == 0 ,"Variables should not containg missings!")
                )
                "Added the variables!"
            })
           validate(
                  need((all(isolate(input$custom_columns) %in% colnames(data))), "Variables should be in samplesheet!"),
                  need(all(purrr::map_lgl(data[,isolate(input$custom_columns)], .f = is.logical)) ,"Variables should be logical (TRUE/FALSE)!"),
                  need(sum(is.na(data[,isolate(input$custom_columns)])) == 0 ,"Variables should not containg missings!")
                )
          # Update outlier variables
          new <- c(outlier_variables, isolate(input$custom_columns))
          values$outlier_vars <- new
        })

      # Reset custom outlier variables
      observeEvent(input$reset_custom_columns,{
        values$outlier_vars <- outlier_variables
        updateSelectInput(session, inputId = "custom_columns", selected = "")
      })

      ### Color by variables

      # Variables that can be selected to color the plots,
      # When define_outliers is selected, extra options will become available as to color the plots by specific QC metric outliers
      # In addition, for every metric there is an 'unique' version that is TRUE if the sample only fails on that metric
      observe({
           updateSelectInput(session, inputId = "colorbyQCmetric", choices = c(values$outlier_vars, paste0(values$outlier_vars[values$outlier_vars != "None"], "_unique")))
         })

      # If the user chooses a metric to color the samples by, the 'scatterColorByColumn' input will be set to 'None'
      observe({
           if(input$colorbyQCmetric != "None") {
              updateSelectInput(session, inputId = "scatterColorByColumn", selected = "None")
           }
         })

      # If define_outliers == 0, the 'colorbyQCmetric' input will be set to 'None'
       observe({
           if(input$define_outliers == 0) {
              if(input$scatterColorByColumn != "None") {
                values$colorscheme = set_color_scheme(plotdata = plotdata(),
                                                      column = input$scatterColorByColumn,
                                                      colorscheme = values[["selected_colorscheme"]])
              }
              updateSelectInput(session, inputId = "colorbyQCmetric", selected = "None")
           } else if(input$define_outliers == 1) {
             if(input$scatterColorByColumn != "None") {
               values$colorscheme = set_color_scheme(plotdata = plotdata(),
                                                     column = input$scatterColorByColumn,
                                                     colorscheme = values[["selected_colorscheme"]])
             } else if(input$colorbyQCmetric != "None") {
               values$colorscheme = set_color_scheme(plotdata = plotdata(),
                                                     column = input$colorbyQCmetric,
                                                     colorscheme = values[["selected_colorscheme"]])
             }
           }
         })

       observeEvent(input$scatterColorByColumn,{
         if(input$scatterColorByColumn != "None" && input$colorbyQCmetric == "None")
         values$colorscheme = set_color_scheme(plotdata = plotdata(),
                                               column = input$scatterColorByColumn,
                                               colorscheme = values[["selected_colorscheme"]])

       })

       observeEvent(input$colorbyQCmetric,{
         if(input$colorbyQCmetric != "None") {
           values$colorscheme = set_color_scheme(plotdata = plotdata(),
                                                 column = input$colorbyQCmetric,
                                                 colorscheme = values[["selected_colorscheme"]])
         }
         # } else if(input$colorbyQCmetric == "None" && input$scatterColorByColumn != "None")
         #   values$colorscheme = set_color_scheme(plotdata = plotdata(),
         #                                         column = input$scatterColorByColumn,
         #                                         colorscheme = "viridis")

       })

      ### Sample QC
      # Reset sample qc thresholds to default
      observeEvent(input$reset_defaults, {
          updateNumericInput(session, "XY_diff_main", value = thresholds$XY_diff)
          updateNumericInput(session, "XY_diff_tab", value = thresholds$XY_diff)
          values$XY_diff <- thresholds$XY_diff
          updateNumericInput(session, "MU_threshold_main", value = thresholds$MU)
          updateNumericInput(session, "MU_threshold_tab", value = thresholds$MU)
          values$MU_threshold <- thresholds$MU
          updateNumericInput(session, "RG_ratio_threshold_main", value = thresholds$RG_ratio)
          updateNumericInput(session, "RG_ratio_threshold_tab", value = thresholds$RG_ratio)
          values$RG_ratio_threshold <- thresholds$RG_ratio
          updateNumericInput(session, "GR_ratio_threshold_main", value = thresholds$GR_ratio)
          updateNumericInput(session, "GR_ratio_threshold_tab", value = thresholds$GR_ratio)
          values$GR_ratio_threshold <- thresholds$GR_ratio
          updateNumericInput(session, "OP_threshold_main", value = thresholds$OP)
          updateNumericInput(session, "OP_threshold_tab", value = thresholds$OP)
          values$OP_threshold <- thresholds$OP
          updateNumericInput(session, "HC_threshold_main", value = thresholds$HC)
          updateNumericInput(session, "HC_threshold_tab", value = thresholds$HC)
          values$HC_threshold <- thresholds$HC
          updateNumericInput(session, "bscon_threshold_main", value = thresholds$bscon)
          updateNumericInput(session, "bscon_threshold_tab", value = thresholds$bscon)
          values$bscon_threshold <- thresholds$bscon
          updateNumericInput(session, "detP_threshold_main", value = thresholds$detP)
          updateNumericInput(session, "detP_threshold_tab", value = thresholds$detP)
          values$detP_threshold <- thresholds$detP
          updateNumericInput(session, "beadNr_threshold_main", value = thresholds$beadNr)
          updateNumericInput(session, "beadNr_threshold_tab", value = thresholds$beadNr)
          values$beadNr_threshold <- thresholds$beadNr
          updateNumericInput(session, "IBS_mean_main", value = thresholds$IBS_mean)
          updateNumericInput(session, "IBS_mean_tab", value = thresholds$IBS_mean)
          values$IBS_mean <- thresholds$IBS_mean
          updateNumericInput(session, "IBS_var_main", value = thresholds$IBS_var)
          updateNumericInput(session, "IBS_var_tab", value = thresholds$IBS_var)
          values$IBS_var <- thresholds$IBS_var

          values$controlPCA_outliers <- if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) TRUE else FALSE
          updateCheckboxInput(session, "controlPCA_outliers_main", value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) TRUE else FALSE)
          updateCheckboxInput(session, "controlPCA_outliers_tab", value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) TRUE else FALSE)

          values$controlPCs_used <- if(!is.null(thresholds$controlPCs_used) && !is.na(thresholds$controlPCs_used)) thresholds$controlPCs_used else "PC1"
          updateSelectInput(session, "controlPCs_used_main", selected = if(!is.null(thresholds$controlPCs_used) && !is.na(thresholds$controlPCs_used)) thresholds$controlPCs_used else "PC1")
          updateSelectInput(session, "controlPCs_used_tab", selected = if(!is.null(thresholds$controlPCs_used) && !is.na(thresholds$controlPCs_used)) thresholds$controlPCs_used else "PC1")

          values$controlPCA_sd <- if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) thresholds$controlPCA_sd else 3
          updateNumericInput(session, "controlPCA_sd_main", value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) thresholds$controlPCA_sd else 3)
          updateNumericInput(session, "controlPCA_sd_tab", value = if(!is.null(thresholds$controlPCA_sd) && !is.na(thresholds$controlPCA_sd)) thresholds$controlPCA_sd else 3)

          values$PCA_outliers <- if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) TRUE else FALSE
          updateCheckboxInput(session, "PCA_outliers_main", value = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) TRUE else FALSE)
          updateCheckboxInput(session, "PCA_outliers_tab", value = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) TRUE else FALSE)

          values$PCs_used <- if(!is.null(thresholds$PCs_used) && !is.na(thresholds$PCs_used)) thresholds$PCs_used else "PC1"
          updateSelectInput(session, "PCs_used_main", selected = if(!is.null(thresholds$PCs_used) && !is.na(thresholds$PCs_used)) thresholds$PCs_used else "PC1")
          updateSelectInput(session, "PCs_used_tab", selected = if(!is.null(thresholds$PCs_used) && !is.na(thresholds$PCs_used)) thresholds$PCs_used else "PC1")

          values$PCA_sd <- if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) thresholds$PCA_sd else 3
          updateNumericInput(session, "PCA_sd_main", value = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) thresholds$PCA_sd else 3)
          updateNumericInput(session, "PCA_sd_tab", value = if(!is.null(thresholds$PCA_sd) && !is.na(thresholds$PCA_sd)) thresholds$PCA_sd else 3)
        })

        # Reset probe qc thresholds to default
        observeEvent(input$reset_defaults_probes, {
          updateNumericInput(session, "detP_probes", value = 0.05)
          updateNumericInput(session, "beadNr_probes", value = 0.05)
        })

      # Vector of sample outliers
      outliers <- reactive({
        if(input$define_outliers == 1) {
          outliers <- outliers_table()
          outliers <- outliers %>% dplyr::filter(Total) %$% Sample_Name
          outliers
        } else {
          outliers <- c()
          outliers
        }
      })

      # Table of custom outliers
      custom_outliers <- reactive({
        diff <- setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier"))
        if(length(diff) > 0) {
          custom_outlier_table <- data[,c("Sample_Name", diff)]
          custom_outlier_table_tmp <- custom_outlier_table %>% dplyr::select(-Sample_Name)
          total <- apply(custom_outlier_table_tmp, 1, function(x) any(x))
          custom_outlier_table_tmp$Total <- total
          custom_outlier_table_tmp
        } else {
          custom_outlier_table <- tibble::tibble(Sample_Name = character(), Total = logical())
        }
        })


      # Creates a table of all samples and on which QC Metric they failed
      outliers_table <- reactive({
        if(!is.null(qcmetrics@PCs)) {
          PCs_melt_sel <- melt_PCs()
        }

        controlPCs_melt_sel <- melt_controlPCs()

        if(input$define_outliers == 1) {
          MU_outliers <- data %>% dplyr::filter(Methylated < values$MU_threshold | Unmethylated < values$MU_threshold) %$% Sample_Name
          RG_ratio_outliers <- data %>% dplyr::filter(RG_ratio < values$RG_ratio_threshold) %$% Sample_Name
          GR_ratio_outliers <- data %>% dplyr::filter(GR_ratio < values$GR_ratio_threshold) %$% Sample_Name
          OP_outliers <- data %>% dplyr::filter(OP_x < values$OP_threshold) %$% Sample_Name
          HC_outliers <- data %>% dplyr::filter(HC_x < values$HC_threshold) %$% Sample_Name
          bscon_outliers <- data %>% dplyr::filter(bscon < values$bscon_threshold) %$% Sample_Name
          detP_outliers <- data %>% dplyr::filter(detP_Percentage > values$detP_threshold) %$% Sample_Name
          beadNr_outliers <- data %>% dplyr::filter(BeadNr_Percentage > values$beadNr_threshold) %$% Sample_Name
          sex_fails <- data %>% dplyr::filter(!is.na(XYdiff) & (XYdiff < values$XY_diff & Sex == "M") | (XYdiff > values$XY_diff & Sex == "F")) %$% Sample_Name
          sex_missings <- data %>% dplyr::filter(is.na(Sex)) %$% Sample_Name
          outlier_table <- tibble::tibble(Sample_Name = data$Sample_Name)
          outlier_table <- outlier_table %>%
                              dplyr::mutate(MU_outlier = ifelse(Sample_Name %in% MU_outliers, TRUE, FALSE),
                                     RG_ratio_outlier = ifelse(Sample_Name %in% RG_ratio_outliers, TRUE, FALSE),
                                     GR_ratio_outlier = ifelse(Sample_Name %in% GR_ratio_outliers, TRUE, FALSE),
                                     OP_outlier = ifelse(Sample_Name %in% OP_outliers, TRUE, FALSE),
                                     HC_outlier = ifelse(Sample_Name %in% HC_outliers, TRUE, FALSE),
                                     bscon_outlier = ifelse(Sample_Name %in% bscon_outliers, TRUE, FALSE),
                                     detectionP_outlier = ifelse(Sample_Name %in% detP_outliers, TRUE, FALSE),
                                     beadNr_outlier = ifelse(Sample_Name %in% beadNr_outliers, TRUE, FALSE),
                                     sex_outlier = ifelse(Sample_Name %in% sex_fails, TRUE, FALSE),
                                     sex_outlier = ifelse(Sample_Name %in% sex_missings, TRUE, sex_outlier))

          # If selected: add PCA outliers
          if(values$PCA_outliers == 1) {
            PCA_outliers <- unique(PCs_melt_sel[PCs_melt_sel$outlier,]$Sample_Name)
            outlier_table <- outlier_table %>% dplyr::mutate(PCA_outlier = ifelse(Sample_Name %in% PCA_outliers, TRUE, FALSE))
          }

          # If selected: add controlPCA outliers
          if(values$controlPCA_outliers == 1) {
            controlPCA_outliers <- unique(controlPCs_melt_sel[controlPCs_melt_sel$outlier,]$Sample_Name)
            outlier_table <- outlier_table %>% dplyr::mutate(controlPCA_outlier = ifelse(Sample_Name %in% controlPCA_outliers, TRUE, FALSE))
          }

          # Optional: genotype outliers (only when genotype data is available)
          if(!is.null(qcmetrics@ibs$ibs_geno_identical)) {
            IBS_outliers <- data %>% dplyr::filter(IBS_mean < values$IBS_mean | IBS_var > values$IBS_var) %$% Sample_Name
            outlier_table <- outlier_table %>% dplyr::mutate(IBS_outlier = ifelse(Sample_Name %in% IBS_outliers, TRUE, FALSE))
          }
          # Optional: custom outliers (if specified)
          diff <- setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))
          if(length(diff) > 0) {
            outlier_table <- outlier_table %>% dplyr::left_join(data[,c("Sample_Name",diff)], by = "Sample_Name")
          }

          # Create 'Total' column
          outlier_table_tmp <- outlier_table %>% dplyr::select(-Sample_Name)
          total <- apply(outlier_table_tmp, 1, function(x) any(x))
          outlier_table$Total <- total
          outlier_table

        } else {
          outlier_table = tibble::tibble(Sample_Name = character(), MU_outlier = logical(),
                RG_ratio_outlier = logical(), GR_ratio_outlier = logical(),
                OP_outlier = logical(),  HC_outlier = logical(), bscon_outlier = logical(), detectionP_outlier = logical(), beadNr_outlier = logical(),
                sex_outlier = logical(), PCA_outlier = logical(), controlPCA_outlier = logical(), IBS_outlier = logical(), Total = logical())
        }
      })

      # Table used for plotting, includes 'unique' columns: samples that fail only on specified metric
      outliers_table_plot <- reactive({
        if(input$define_outliers == 1) {
            outlier_vars <- colnames(outliers_table())[!colnames(outliers_table()) %in% c("Sample_Name", "Total")]
            outlier_table_temp <- outliers_table()[,outlier_vars]

            add_unique_outliers <- function(var) {
              outliers_vars_current <- outlier_table_temp[,colnames(outlier_table_temp) == var]
              outliers_vars_other <- outlier_table_temp[,colnames(outlier_table_temp) != var]
              other_outliers <- apply(outliers_vars_other, 1, function(x) any(x))
              unique <- outliers_vars_current & !other_outliers
              unique
            }
            uniques <- lapply(outlier_vars, FUN = add_unique_outliers)
            uniques <- do.call(cbind, uniques)
            colnames(uniques) <- paste0(outlier_vars, "_unique")
            outliers_table_plot <- cbind(outliers_table(), uniques)
            outliers_table_plot
        } else {
          outlier_table = tibble::tibble()
        }
      })

      ## PlotData
      plotdata <- reactive({
        if(input$define_outliers == 1) {
        inter <- setdiff(colnames(outliers_table_plot()), colnames(data))
        inter <- c("Sample_Name", inter)
        plotdata <- data %>%
          dplyr::mutate(predictedSex = ifelse(XYdiff < values$XY_diff, "F", "M")) %>%
          dplyr::left_join(outliers_table_plot()[,inter], by = "Sample_Name")
        plotdata
        } else {
          plotdata <- data %>% dplyr::mutate(predictedSex = ifelse(XYdiff < values$XY_diff, "F", "M"))
          plotdata
        }
      })

      ## Relatedness

      # Table of related samples (based on DNAm x Genotype data)
      related_samples_geno <- reactive({
          relateds <- qcmetrics@ibs$ibs_geno_relatedness %>%
            dplyr::filter(IBS_mean > input$ibs_relatedness_mean, IBS_var < input$ibs_relatedness_var)
          relateds
      })

      # Table of related samples (based on DNAm x DNAm data)
      related_samples_dnam <- reactive({
          related_nogeno <- qcmetrics@ibs$ibs_dnam %>%
            dplyr::filter(IBS_mean > input$ibs_relatedness_mean, IBS_var < input$ibs_relatedness_var)
          related_nogeno
      })

      # Download handlers
      output$downloadData <- downloadHandler(
        filename = function() {
          if(is.null(qcmetrics@metadata$name)) {
            "sample_outlier_table.rds"
          } else {
            sprintf("sample_outlier_table_%s.rds", qcmetrics@metadata$name)
          }
        },
        content = function(file) {
          dat <- list()
          dat$sample_outliers <- outliers_table()
          dat$thresholds <- list(MU = values$MU_threshold,
                                 RG_ratio = values$RG_ratio_threshold,
                                 GR_ratio = values$GR_ratio_threshold,
                                 OP = values$OP_threshold,
                                 HC = values$HC_threshold,
                                 bscon = values$bscon_threshold,
                                 detP = values$detP_threshold,
                                 beadNr = values$beadNr_threshold,
                                 XY_diff = values$XY_diff,
                                 IBS_mean = values$IBS_mean,
                                 IBS_var = values$IBS_var,
                                 PCA_sd = if(values$PCA_outliers == 1) thresholds$PCA_sd else NA,
                                 PCs_used = if (values$PCA_outliers == 1) values$PCs_used else NA,
                                 controlPCA_sd = if(values$controlPCA_outliers == 1) values$controlPCA_sd else NA,
                                 controlPCs_used = if (values$controlPCA_outliers == 1) values$controlPCs_used else NA
                                 )
          dat$package_version <- packageVersion("QCpG")
          saveRDS(dat, file = file)
        }
      )

      output$newqcmetrics <- downloadHandler(
        filename = function() {
          "qcmetrics_updated.rds"
        },
        content = function(file) {
          metrics <- qcmetrics # copy of original metrics-file

          metrics@Sample_Metrics <- metrics@Sample_Metrics %>%
            dplyr::mutate(predictedSex = ifelse(XYdiff < values$XY_diff, "F", "M"))

          metrics@samplesheet <- samplesheet

          # Updates based on potentially updated thresholds
          metrics@relatedness_thresholds <- list(IBS_mean = input$ibs_relatedness_mean,
                                                 IBS_var = input$ibs_relatedness_var)

          diff <- setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))
          metrics@sample_qc_outliers <- list(
            sample_outliers = outliers_table(),
            thresholds = list(MU = values$MU_threshold,
                              RG_ratio = values$RG_ratio_threshold,
                              GR_ratio = values$GR_ratio_threshold,
                              OP = values$OP_threshold,
                              HC = values$HC_threshold,
                              bscon = values$bscon_threshold,
                              detP = values$detP_threshold,
                              beadNr = values$beadNr_threshold,
                              XY_diff = values$XY_diff,
                              IBS_mean = values$IBS_mean,
                              IBS_var = values$IBS_var,
                              PCA_sd = if(values$PCA_outliers == 1) values$PCA_sd else NA,
                              PCs_used = if (values$PCA_outliers == 1) values$PCs_used else NA,
                              controlPCA_sd = if(values$controlPCA_outliers == 1) values$controlPCA_sd else NA,
                              controlPCs_used = if (values$controlPCA_outliers == 1) values$controlPCs_used else NA,
                              custom_outliers = if(length(diff) > 0) diff else NA
                              )

          )
          thresholds = list(
            detP = input$detP_probes,
            beadNr = input$beadNr_probes
          )
          if(check$default == 1) {
            thresholds$zhou_population = zhou_pop_used()
          } else if (!is.null(probe_thresholds[["zhou_population"]])) {
            thresholds$zhou_population = probe_thresholds[["zhou_population"]]
          } else {
            zhou_pop_used <- NULL
          }
          metrics@probe_qc_outliers <- list(
            probe_outliers = probe_outlier_table(),
            thresholds = thresholds
          )

          if(length(values$selected_samples) > 0) {
            metrics@selected_samples <- values$selected_samples
          } else {
            metrics@selected_samples <- NA_character_
          }

          saveRDS(metrics, file = file)
        }
      )

      output$newqcmetrics_probes <- downloadHandler(
        filename = function() {
          "qcmetrics_updated.rds"
        },
        content = function(file) {
          metrics <- qcmetrics # copy of original metrics-file

          metrics@Sample_Metrics <- metrics@Sample_Metrics %>%
            dplyr::mutate(predictedSex = ifelse(XYdiff < values$XY_diff, "F", "M"))

          metrics@samplesheet <- samplesheet

          # Updates based on potentially updated thresholds
          metrics@relatedness_thresholds <- list(IBS_mean = input$ibs_relatedness_mean,
                                                 IBS_var = input$ibs_relatedness_var)

          diff <- setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))
          metrics@sample_qc_outliers <- list(
            sample_outliers = outliers_table(),
            thresholds = list(MU = values$MU_threshold,
                              RG_ratio = values$RG_ratio_threshold,
                              GR_ratio = values$GR_ratio_threshold,
                              OP = values$OP_threshold,
                              HC = values$HC_threshold,
                              bscon = values$bscon_threshold,
                              detP = values$detP_threshold,
                              beadNr = values$beadNr_threshold,
                              XY_diff = values$XY_diff,
                              IBS_mean = values$IBS_mean,
                              IBS_var = values$IBS_var,
                              PCA_sd = if(values$PCA_outliers == 1) values$PCA_sd else NA,
                              PCs_used = if (values$PCA_outliers == 1) values$PCs_used else NA,
                              controlPCA_sd = if(values$controlPCA_outliers == 1) values$controlPCA_sd else NA,
                              controlPCs_used = if (values$controlPCA_outliers == 1) values$controlPCs_used else NA,
                              custom_outliers = if(length(diff) > 0) diff else NA
            )

          )
          thresholds = list(
            detP = input$detP_probes,
            beadNr = input$beadNr_probes
          )
          if(check$default == 1) {
            thresholds$zhou_population = zhou_pop_used()
          } else if (!is.null(probe_thresholds[["zhou_population"]])) {
            thresholds$zhou_population = probe_thresholds[["zhou_population"]]
          } else {
            zhou_pop_used <- NULL
          }
          metrics@probe_qc_outliers <- list(
            probe_outliers = probe_outlier_table(),
            thresholds = thresholds
          )

          if(length(values$selected_samples) > 0) {
            metrics@selected_samples <- values$selected_samples
          } else {
            metrics@selected_samples <- NA_character_
          }

          saveRDS(metrics, file = file)
        }
      )

      # output$newqcmetrics <- downloadHandler(
      #   filename = function() {
      #     "qcmetrics_updated.rds"
      #   },
      #   content = .new_qcmetrics
      # )

      # Selected samples
      output$download_selected_samples <- downloadHandler(
        filename = function() {
          if(is.null(qcmetrics@metadata$name)) {
            "selected_samples.txt"
          } else {
            sprintf("selected_samples_%s.txt", qcmetrics@metadata$name)
          }

        },
        content = function(file) {
          readr::write_lines(values$selected_samples, path = file)
        }
      )

      # Download related samples, both DNAmxGeno and and DNAmxDNAm will be included
      output$download_related_samples <- downloadHandler(
        filename = function() {
          if(is.null(qcmetrics@metadata$name)) {
            "related_samples.rds"
          } else {
            sprintf("related_samples_%s.rds", qcmetrics@metadata$name)
          }

        },
        content = function(file) {
          if(!is.null(qcmetrics@ibs$ibs_geno_identical) && !is.null(qcmetrics@ibs$ibs_dnam)) {
            related_samples_geno <- related_samples_geno() %>% dplyr::select(Sample_Name, Sample_Name_geno, IBS_mean, IBS_var)
            related_samples_dnam <- related_samples_dnam() %>% dplyr::select(Sample_Name.x, Sample_Name.y, IBS_mean, IBS_var)
            related_samples <- list(related_samples_geno = related_samples_geno, related_samples_dnam = related_samples_dnam)
          } else if(!is.null(qcmetrics@ibs$ibs_geno_identical) && is.null(qcmetrics@ibs$ibs_dnam)) {
            related_samples_geno <- related_samples_geno() %>% dplyr::select(Sample_Name, Sample_Name_geno, IBS_mean, IBS_var)
            related_samples <- list(related_samples_geno = related_samples_geno, related_samples_dnam = NULL)
          } else {
            related_samples_dnam <- related_samples_dnam() %>% dplyr::select(Sample_Name.x, Sample_Name.y, IBS_mean, IBS_var)
            related_samples <- list(related_samples_geno = NULL, related_samples_dnam = related_samples_dnam)
          }

          saveRDS(related_samples, file = file)
        }
      )


      ### Probe QC

      ## Cross-reactive/SNP-probes by Zhou et al.

      # Little workaround to set a default value when 'download' button has not yet been pressed
      # Adapted from: https://stackoverflow.com/questions/33662033/shiny-how-to-make-reactive-value-initialize-with-default-value
      check <- reactiveValues(default = 0)

      observeEvent(input$download_zhou,{
        check$default <- input$download_zhou
      })

      zhou_masking <- reactive({
        if(check$default == 0) {
          if(is.null(probe_thresholds[["zhou_population"]])) {
            c()
          } else {
            qcmetrics@probe_qc_outliers$probe_outliers %>% dplyr::filter(zhou_outlier) %$% Probe
          }
        } else {
          zhou_masking_temp()
        }
      })

      # Download annotation and select specified column
      zhou_masking_temp <- eventReactive(
        input$download_zhou, {
          if(array == "450k") {
            mask <- data.frame(readr::read_tsv("https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.pop.tsv.gz"))
            mask <- mask[mask[[input$zhou_column]],][["probeID"]]
            mask
          } else {
            mask <- data.frame(readr::read_tsv("https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.pop.tsv.gz"))
            mask <- mask[mask[[input$zhou_column]],][["probeID"]]
            mask
          }
        })

      # Tag which population was used (for output)
      zhou_pop_used <- eventReactive(
        input$download_zhou, {
          pop <- input$zhou_column
          pop
        }
      )

      ## Creates a table of all probes and on which QC Metric they failed
      probe_outlier_table <- reactive({
        if(input$define_outliers_probes == 1) {
            detp_probes_outliers <- qcmetrics@Probe_Metrics %>% dplyr::filter(detP_Percentage > input$detP_probes) %$% Probe
            beadnr_probes_outliers <- qcmetrics@Probe_Metrics %>% dplyr::filter(BeadNr_Percentage > input$beadNr_probes) %$% Probe

            if(check$default == 1 || !is.null(probe_thresholds[["zhou_population"]])) {
              masking_outliers <- zhou_masking()
              probe_outlier_table <- tibble::tibble(Probe = qcmetrics@Probe_Metrics$Probe,
                                                    detP_outlier = ifelse(Probe %in% detp_probes_outliers, TRUE, FALSE),
                                                    beadNr_outlier = ifelse(Probe %in% beadnr_probes_outliers, TRUE, FALSE),
                                                    zhou_outlier = ifelse(Probe %in% masking_outliers, TRUE, FALSE),
                                                    Total = detP_outlier | beadNr_outlier | zhou_outlier)
              probe_outlier_table %>% dplyr::filter(Total)
            } else {
              probe_outlier_table <- tibble::tibble(Probe = qcmetrics@Probe_Metrics$Probe,
                                                    detP_outlier = ifelse(Probe %in% detp_probes_outliers, TRUE, FALSE),
                                                    beadNr_outlier = ifelse(Probe %in% beadnr_probes_outliers, TRUE, FALSE),
                                                    Total = detP_outlier | beadNr_outlier)
              probe_outlier_table %>% dplyr::filter(Total)
            }

        } else {
            probe_outlier_table <- tibble::tibble(Probe = character(),
                                          detP_outlier = logical(),
                                          beadNr_outlier = logical(),
                                          zhou_outlier = logical(),
                                          Total = logical())
            probe_outlier_table
        }
      })

      ## Vector of probe outliers
       probe_outliers <- reactive({
        if(input$define_outliers_probes == 1) {
          outliers <- probe_outlier_table() %$% Probe
          outliers
        } else {
          outliers <- c()
          outliers
        }
      })

      ## Download handlers
      output$download_probe_table <- downloadHandler(
         filename = function() {
           if(is.null(qcmetrics@metadata$name)) {
             "probe_outlier_table.rds"
           } else {
             sprintf("probe_outlier_table_%s.rds", qcmetrics@metadata$name)
           }
         },
         content = function(file) {
           dat <- list()
           dat$probe_outliers <- probe_outlier_table()
           dat$thresholds <- list(detP = input$detP_probes, beadNr = input$beadNr_probes)
           if(check$default == 1) {
             dat$thresholds$zhou_population = zhou_pop_used()
           } else if (!is.null(probe_thresholds[["zhou_population"]])) {
             dat$thresholds$zhou_population = probe_thresholds[["zhou_population"]]
           } else {
             zhou_pop_used <- NULL
           }
           dat$probe_qc_method <- qcmetrics@probe_qc_method
           dat$package_version <- packageVersion("QCpG")
           saveRDS(dat, file = file)
         }
       )


      ### Interactive UI elements

      # Number of sex missings
      output$samples.missing.sex.text <- renderUI({
        if(length(samples_missing_sex) > 1) {
              textOutput("samples.missing.sex.multiple")
         } else if(length(samples_missing_sex) == 1) {
              textOutput("samples.missing.sex.1")
         }
      })
      # Reset selected samples (appears only when samples are selected)
      output$reset_selected_samples <- isolate(renderUI({
        if(length(values$selected_samples) > 0)
          tagList(
            h5(sprintf("%s samples selected", length(values$selected_samples))),
            actionButton("reset_samples", "Reset selected samples"),
            br(),
            br(),
            downloadButton("download_selected_samples", "Download selected samples")
            )
      }))

      # Show background (appears only if background is specified)
      output$show_background <- renderUI({
        if(!is.null(background))
          checkboxInput(inputId = "show_background", label = "Show background?", value = TRUE)
      })

      # Custom outliers
      output$custom_outliers_sidepanel_main <- renderUI({
       if(length(setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))) > 0) {
         HTML(paste(c("<b>Custom outlier columns:</b>", setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))), collapse = " <br> "))
        }
      })
      output$custom_outliers_sidepanel_tab <- renderUI({
        if(length(setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))) > 0) {
          HTML(paste(c("<b>Custom outlier columns:</b>", setdiff(values$outlier_vars, c(outlier_variables, "PCA_outlier", "controlPCA_outlier"))), collapse = " <br> "))
        }
      })

      ### Reactive Text
      output$note_masking <- renderText({
          "*Note that for every population there is a 'General' column and a 'snp5' column. The 'snp5' column contains probes with a SNP with MAF > 0.01
           within 5 bp of the 3' end of the probe for the chosen population. The 'General' column merged the 'snp5' column with all other masking recommendations,
           such as cross-reactive probes and probes with low mapping quality"
        })

      output$note_download_zhou <- renderText({
          "Download masking info, takes a few seconds.."
        })

      output$note_probeqc <- renderText({
          if(qcmetrics@probe_qc_method == "post_sampleqc")
            "Note that probe QC metrics were calculated after removing sample outliers based on thresholds specified in the getQCmetrics function"
          })
      output$note_probeqc2 <- renderText({
              if(qcmetrics@probe_qc_method == "post_sampleqc")
                "Note that probe QC metrics were calculated after removing sample outliers based on thresholds specified in the getQCmetrics function."
              })
      output$nr_probes_zhou <- renderText({
          nr <- length(zhou_masking())
          sprintf("%s probes are masked based on Zhou's masking recommendations.", nr)
        })
      output$zhou_population_used <- renderText({
        if(check$default == 1) {
          zhou_pop_used <- zhou_pop_used()
        } else if(!is.null(probe_thresholds[["zhou_population"]])) {
          zhou_pop_used <- probe_thresholds[["zhou_population"]]
        } else {
          zhou_pop_used <- "None"
        }
        sprintf("Population used for masking: %s", zhou_pop_used)
      })

      ## Total number of outliers (identical for each tab)
      output$nr.outliers.total.pca <- output$nr.outliers.total.controlpca <- output$nr.outliers.total.sex  <-
        output$nr.outliers.total.mu <- output$nr.outliers.total.rg <- output$nr.outliers.total.op <-
        output$nr.outliers.total.hc <- output$nr.outliers.total.bscon <-
        output$nr.outliers.total.detp <- output$nr.outliers.total.beadnr <- output$nr.outliers.total.ibs.intra <-
        output$nr.outliers.total.ibs.inter <- output$nr.outliers.total.ctf <- output$nr.outliers.total.age <-
        output$nr.outliers.total.smoking <- output$nr.outliers.total.custom <-  renderText({
          if(input$define_outliers) {
            nr <- length(outliers())
            sprintf("In total, there are %s outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
          }
        })

      output$nr.outliers.pca <- renderText({
        if(values$PCA_outliers == 1 && input$define_outliers == 1) {
          nr <- nrow(outliers_table() %>% dplyr::filter(PCA_outlier))
          sprintf("There are %s PCA outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })

      output$nr.probes.used.pca <- renderText({
        if(!is.null(qcmetrics@PCs)) {
          sprintf("%s probes were used to estimate the PCs", qcmetrics@nprobes_PCA)
        }
      })

      output$nr.outliers.controlpca <- renderText({
        if(values$controlPCA_outliers == 1 && input$define_outliers == 1) {
          nr <- nrow(outliers_table() %>% dplyr::filter(controlPCA_outlier))
          sprintf("There are %s control PCA outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })

      output$sex.outliers <- renderText({
        if(input$define_outliers == 1) {
          nr <- nrow(outliers_table() %>% dplyr::filter(sex_outlier))
          sprintf("There are %s sex outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$samples.missing.sex.1 <- renderText({
        "Note: There is 1 sample for which reported Sex is unavailable.
                   This sample will be flagged as an outlier."

      })
      output$samples.missing.sex.multiple <- renderText({
        sprintf("Note: There are %s samples for which reported Sex is unavailable.
                   These samples will be flagged as outliers.", length(samples_missing_sex))

      })
      output$MU.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(MU_outlier))
          sprintf("There are %s MU outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }

      })
      output$RG.ratio.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(RG_ratio_outlier))
          sprintf("There are %s RG-ratio outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$GR.ratio.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(GR_ratio_outlier))
          sprintf("There are %s GR-ratio outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$OP.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(OP_outlier))
          sprintf("There are %s OP outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })

      output$HC.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(HC_outlier))
          sprintf("There are %s HC outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$bscon.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(bscon_outlier))
          sprintf("There are %s bscon outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$detp.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(detectionP_outlier))
          sprintf("There are %s detP outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
       })
      output$samples.missing.geno <- renderText({
        if(!is.null(qcmetrics@ibs$ibs_geno_identical)) {
          sprintf("Note: There are %s samples for which genotype is unavailable.", length(samples_missing_geno))
        }
       })
      output$ibs.outliers <- renderText({
        if(!is.null(qcmetrics@ibs$ibs_geno_identical) & input$define_outliers == 1) {
          nr <- nrow(outliers_table() %>% dplyr::filter(IBS_outlier))
          sprintf("There are %s IBS outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
       })
      output$nr.called.probes.geno <- renderText({
        if(!is.null(qcmetrics@ibs$called_probes_geno)) {
          sprintf("IBS was calculated using %s SNP-containing probes", length(qcmetrics@ibs$called_probes_geno))
        }
      })
      output$nr.called.probes.dnam <- renderText({
        if(!is.null(qcmetrics@ibs$called_probes)) {
          sprintf("Relatedness was calculated using %s SNP-containing probes", length(qcmetrics@ibs$called_probes))
        }
      })
      output$ibs.summarized.text <- renderText({
        if(!qcmetrics@ibs$save_all_relationships) {
          sprintf("Note: all unrelated samples have been presummarized.\nHence the minimum value for IBS_mean=%s and the maximum value for IBS_var=%s",
                  qcmetrics@relatedness_thresholds[["IBS_mean"]], qcmetrics@relatedness_thresholds[["IBS_var"]])
        }
      })

      output$age.text <- renderText({
        if(!("Age" %in% colnames(samplesheet))) {
          "No 'Age' column was found in the samplesheet. Therefore, a histogram of predicted age will be displayed instead of a scatterplot"
        }
      })
      output$age.text.missing <- renderText({
        if(("Age" %in% colnames(samplesheet)) && sum(is.na(samplesheet$Age)) > 0) {
          sprintf("Chronological age is missing for %s samples", sum(is.na(samplesheet$Age)))
        }
      })
      output$smoking.text <- renderText({
        if(!("Smoking_Status" %in% colnames(samplesheet))) {
          "No 'Smoking_Status' column was found in the samplesheet. Therefore, a histogram of the smoking score will be displayed instead of boxplots"
        }
      })
      output$smoking.text.missing <- renderText({
        if(("Smoking_Status" %in% colnames(samplesheet)) && sum(is.na(samplesheet$Smoking_Status)) > 0) {
        sprintf("Smoking status is missing for %s samples",  sum(is.na(samplesheet$Smoking_Status)))
        }
      })
      output$alcohol.text <- renderText({
        if(!("Alcohol" %in% colnames(samplesheet))) {
          "No 'Alcohol' column was found in the samplesheet. Therefore, a histogram of the alcohol PMS will be displayed instead of a scatterplot"
        }
      })
      output$alcohol.text.missing <- renderText({
        if(("Alcohol" %in% colnames(samplesheet)) && sum(is.na(samplesheet$Alcohol)) > 0) {
          sprintf("Alcohol info is missing for %s samples",  sum(is.na(samplesheet$Alcohol)))
        }
      })
      output$bmi.text <- renderText({
        if(!("BMI" %in% colnames(samplesheet))) {
          "No 'BMI' column was found in the samplesheet. Therefore, a histogram of the BMI PMS will be displayed instead of a scatterplot"
        }
      })
      output$bmi.text.missing <- renderText({
        if(("BMI" %in% colnames(samplesheet)) && sum(is.na(samplesheet$BMI)) > 0) {
          sprintf("BMI info is missing for %s samples",  sum(is.na(samplesheet$BMI)))
        }
      })
      output$beadnr.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(outliers_table() %>% dplyr::filter(beadNr_outlier))
          sprintf("There are %s beadnr outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$custom.outliers <- renderText({
        if(input$define_outliers) {
          nr <- nrow(custom_outliers() %>% dplyr::filter(Total))
          sprintf("There are %s custom outliers (%s%%)", nr, signif((nr/nrow(samplesheet)) * 100,2))
        }
      })
      output$nr.outliers.probes <- renderText({
        if(input$define_outliers_probes) {
           nr <- length(probe_outliers())
           sprintf("In total, there are %s outliers (%s%%)", nr, signif((nr/nrow(qcmetrics@Probe_Metrics)) * 100,2))
        }
      })
      output$nr.outliers.probes.detp <- renderText({
        if(input$define_outliers_probes) {
          nr <- sum(probe_outlier_table()$detP_outlier)
          sprintf("There are %s detP outliers (%s%%)", nr, signif((nr/nrow(qcmetrics@Probe_Metrics)) * 100,2))
        }
      })
      output$nr.outliers.probes2 <- renderText({
        if(input$define_outliers_probes) {
          nr <- length(probe_outliers())
          sprintf("In total, there are %s outliers (%s%%)", nr, signif((nr/nrow(qcmetrics@Probe_Metrics)) * 100,2))
        }
      })
      output$nr.outliers.probes.beadnr <- renderText({
        if(input$define_outliers_probes) {
          nr <- sum(probe_outlier_table()$beadNr_outlier)
          sprintf("There are %s beadNr outliers (%s%%)", nr, signif((nr/nrow(qcmetrics@Probe_Metrics)) * 100,2))
        }
      })
      output$nr.outliers.probes3 <- renderText({
        if(input$define_outliers_probes) {
          nr <- length(probe_outliers())
          sprintf("In total, there are %s outliers (%s%%)", nr, signif((nr/nrow(qcmetrics@Probe_Metrics)) * 100,2))
        }
      })


      ###### PLOTS ########

      ## Scatter plot of PCs
      output$PCA <- renderPlotly({

        validate(
          need(!is.null(qcmetrics@PCs), "The metrics-file does not contain PCs!")
        )

        shiny_PCA_plot(scatterPCcolumns = input$scatterPCcolumns,
                 melt_PCs = melt_PCs(),
                 define_outliers = input$define_outliers,
                 webGL = values$webGL,
                 plotdata = plotdata(),
                 colorbyQCmetric = input$colorbyQCmetric,
                 scatterColorByColumn = input$scatterColorByColumn,
                 selected_samples = values$selected_samples,
                 PCA_outliers = values$PCA_outliers,
                 outliers = outliers(),
                 show_outliers = input$show_outliers,
                 colorscheme=values$colorscheme)


      })

      ## Scatter plot of PCs
      output$Heatmap <- renderPlot({

          validate(
            need(!is.null(qcmetrics@heatmap_pvals), message = FALSE)
          )
          heatmap <- qcmetrics@heatmap_pvals %>% #dplyr::mutate(Variable = factor(Variable, levels = variables)) %>%
                      ggplot(aes(x = Variable, y = PC, fill = value)) +
                      geom_tile()  +
                      scale_fill_gradient(low = "white", high = "Red", name = expression(-log[10](italic(p)))) +
                      theme_minimal() +
                      theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                            text = element_text(size = 15),
                                           axis.title.x = element_blank(),
                                           axis.title.y = element_blank())

          heatmap

      })

      output$Pairs <- renderPlot({
        if(input$define_outliers) {
        }
        if(input$pairsplot & input$scatterColorByColumn == "None" & input$colorbyQCmetric == "None") {
          cols <- grep("PC[1-5]{1}$", colnames(plotdata()))
          GGally::ggpairs(plotdata(),diag = list(continuous = "densityDiag" ), upper = "blank", columns = cols)
        } else if(input$pairsplot){
          cols <- grep("PC[1-5]{1}$", colnames(plotdata()))
          GGally::ggpairs(plotdata(),diag = list(continuous = "densityDiag" ), upper = "blank", columns = cols,
                  mapping = aes_string(colour = if(input$colorbyQCmetric != "None") input$colorbyQCmetric else input$scatterColorByColumn))
        }
      })

      output$control_PCA <- renderPlotly({
        validate(
          need(!is.null(qcmetrics@controlPCs), "The metrics-file does not contain control PCs!")
        )
        shiny_controlPCA_plot(scattercontrolPCcolumns = input$scattercontrolPCcolumns,
                     melt_controlPCs = melt_controlPCs(),
                     define_outliers = input$define_outliers,
                     webGL = values$webGL,
                     plotdata = plotdata(),
                     colorbyQCmetric = input$colorbyQCmetric,
                     scatterColorByColumn = input$scatterColorByColumn,
                     selected_samples = values$selected_samples,
                     controlPCA_outliers = values$controlPCA_outliers,
                     outliers = outliers(),
                     show_outliers = input$show_outliers,
                     colorscheme = values$colorscheme
                     )

    })

      output$Heatmap_controlPCA <- renderPlot({

        validate(
          need(!is.null(qcmetrics@heatmap_pvals_controlPCA), message = FALSE)
        )
        heatmap <- qcmetrics@heatmap_pvals_controlPCA %>% #dplyr::mutate(Variable = factor(Variable, levels = variables)) %>%
          ggplot(aes(x = Variable, y = PC, fill = value)) +
          geom_tile()  +
          scale_fill_gradient(low = "white", high = "Red", name = expression(-log[10](italic(p)))) +
          theme_minimal() +
          theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                text = element_text(size = 15),
                axis.title.x = element_blank(),
                axis.title.y = element_blank())

        heatmap

      })

      output$PairsControlPCs <- renderPlot({
        if(input$define_outliers) {
        }
        if(input$pairsplot2 & input$scatterColorByColumn == "None" & input$colorbyQCmetric == "None") {
          cols <- grep("PC[1-5]{1}_control$", colnames(plotdata()))
          GGally::ggpairs(plotdata(),diag = list(continuous = "densityDiag" ), upper = "blank", columns = cols)
        } else if(input$pairsplot2){
          cols <- grep("PC[1-5]{1}_control$", colnames(plotdata()))
          GGally::ggpairs(plotdata(),diag = list(continuous = "densityDiag" ), upper = "blank", columns = cols,
                          mapping = aes_string(colour = if(input$colorbyQCmetric != "None") input$colorbyQCmetric else input$scatterColorByColumn))
        }
      })

      output$xyplot <- renderPlotly({
        shiny_XYplot(define_outliers = input$define_outliers,
               background = background,
               show_background = input$show_background,
               webGL = values$webGL,
               plotdata = plotdata(),
               colorbyQCmetric = input$colorbyQCmetric,
               scatterColorByColumn = input$scatterColorByColumn,
               selected_samples = values$selected_samples,
               outliers = outliers(),
               show_outliers = input$show_outliers,
               bg = if(!is.null(background)) bg_xy else NULL,
               colorscheme = values$colorscheme)
        })

      output$xydiff <- renderPlotly({
        shiny_XYdiff_plot(
               XY_diff = values$XY_diff,
               define_outliers = input$define_outliers,
               webGL = values$webGL,
               plotdata = plotdata(),
               colorbyQCmetric = input$colorbyQCmetric,
               scatterColorByColumn = input$scatterColorByColumn,
               selected_samples = values$selected_samples,
               outliers = outliers(),
               show_outliers = input$show_outliers,
               colorscheme = values$colorscheme
               )
      })

      output$MU <- renderPlotly({
        shiny_plot_general(type = "MU",
                     define_outliers = input$define_outliers,
                     background = background,
                     show_background = input$show_background,
                     webGL = values$webGL,
                     plotdata = plotdata(),
                     colorbyQCmetric = input$colorbyQCmetric,
                     scatterColorByColumn = input$scatterColorByColumn,
                     thresholds = c(values$MU_threshold),
                     selected_samples = values$selected_samples,
                     outliers = outliers(),
                     show_outliers = input$show_outliers,
                     bg = if(!is.null(background)) bg_mu else NULL,
                     colorscheme = values$colorscheme
                     )
        })

      output$RG <- renderPlotly({
        shiny_plot_general(type = "RG",
                     define_outliers = input$define_outliers,
                     background = background,
                     show_background = input$show_background,
                     webGL = values$webGL,
                     plotdata = plotdata(),
                     colorbyQCmetric = input$colorbyQCmetric,
                     scatterColorByColumn = input$scatterColorByColumn,
                     thresholds = c(values$RG_ratio_threshold, values$GR_ratio_threshold),
                     selected_samples = values$selected_samples,
                     outliers = outliers(),
                     show_outliers = input$show_outliers,
                     bg = if(!is.null(background)) bg_rg else NULL,
                     colorscheme = values$colorscheme
                     )
        })

      output$OP <- renderPlotly({
        shiny_plot_general(type = "OP",
                     define_outliers = input$define_outliers,
                     background = background,
                     show_background = input$show_background,
                     webGL = values$webGL,
                     plotdata = plotdata(),
                     colorbyQCmetric = input$colorbyQCmetric,
                     scatterColorByColumn = input$scatterColorByColumn,
                     thresholds = c(values$OP_threshold),
                     selected_samples = values$selected_samples,
                     outliers = outliers(),
                     show_outliers = input$show_outliers,
                     bg = if(!is.null(background)) bg_op else NULL,
                     colorscheme = values$colorscheme
                     )
        })

      output$HC <- renderPlotly({
        shiny_plot_general(type = "HC",
                     define_outliers = input$define_outliers,
                     background = background,
                     show_background = input$show_background,
                     webGL = values$webGL,
                     plotdata = plotdata(),
                     colorbyQCmetric = input$colorbyQCmetric,
                     scatterColorByColumn = input$scatterColorByColumn,
                     thresholds = c(values$HC_threshold),
                     selected_samples = values$selected_samples,
                     outliers = outliers(),
                     show_outliers = input$show_outliers,
                     bg = if(!is.null(background)) bg_hc else NULL,
                     colorscheme = values$colorscheme
                     )
        })

      output$bscon <- renderPlotly({
        shiny_bscon_plot(bscon_threshold = values$bscon_threshold,
                   define_outliers = input$define_outliers,
                   webGL = values$webGL,
                   plotdata = plotdata(),
                   colorbyQCmetric = input$colorbyQCmetric,
                   scatterColorByColumn = input$scatterColorByColumn,
                   selected_samples = values$selected_samples,
                   outliers = outliers(),
                   show_outliers = input$show_outliers,
                   colorscheme = values$colorscheme
                   )
        })

      output$detP <- renderPlotly({
        shiny_detp_beadnr_plot(
          type = "detP",
          threshold = values$detP_threshold,
          threshold2 = detp_threshold,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )
        })

      output$beadNr <- renderPlotly({
        shiny_detp_beadnr_plot(
          type = "beadNr",
          threshold = values$beadNr_threshold,
          threshold2 = beadnr_threshold,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )
        })

       output$intra_ibs <- renderPlotly({

         validate(
           need(!is.null(qcmetrics@ibs$ibs_geno_identical), "No genotype data available!")
         )

         shiny_IBS_plot(IBS_mean = values$IBS_mean,
                  IBS_var = values$IBS_var,
                  define_outliers = input$define_outliers,
                  webGL = values$webGL,
                  plotdata = plotdata(),
                  colorbyQCmetric = input$colorbyQCmetric,
                  scatterColorByColumn = input$scatterColorByColumn,
                  selected_samples = values$selected_samples,
                  outliers = outliers(),
                  show_outliers = input$show_outliers,
                  colorscheme = values$colorscheme
                  )
      })


      output$inter_ibs <- renderPlotly({
        validate(
          need(!is.null(qcmetrics@ibs$ibs_geno_relatedness), "No genotype data available!")
        )
        shiny_relatedness_geno_plot(
                              qcmetrics = qcmetrics,
                              ibs_relatedness_mean = input$ibs_relatedness_mean,
                              ibs_relatedness_var = input$ibs_relatedness_var,
                              define_outliers = input$define_outliers,
                              webGL = values$webGL,
                              plotdata = plotdata(),
                              colorbyQCmetric = input$colorbyQCmetric,
                              scatterColorByColumn = input$scatterColorByColumn,
                              selected_samples = values$selected_samples,
                              outliers = outliers(),
                              show_outliers = input$show_outliers,
                              show_pairs = input$show_pairs,
                              related_samples_geno = related_samples_geno()
                              )
      })

      output$inter_ibs_dnam <- renderPlotly({
        validate(
          need(!is.null(qcmetrics@ibs$ibs_dnam) | !is.null(qcmetrics@ibs$hexagons_relatedness_dnam), "No relatedness data available!")
        )
        shiny_relatedness_dnam_plot(
          qcmetrics = qcmetrics,
          ibs_relatedness_mean = input$ibs_relatedness_mean,
          ibs_relatedness_var = input$ibs_relatedness_var,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers
        )
      })

      output$CTF <- renderPlotly({
        shiny_CTF_plot(
          qcmetrics = qcmetrics,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )

      })

      output$Age <- renderPlotly({
        shiny_age_plot(
          qcmetrics = qcmetrics,
          samplesheet = samplesheet,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          select_age_predictor = input$select_age_predictor,
          colorscheme = values$colorscheme
        )
        })

      output$Smoking <- renderPlotly({
        shiny_smoking_plot(
          qcmetrics = qcmetrics,
          samplesheet = samplesheet,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )
      })

      output$Alcohol <- renderPlotly({
        shiny_plot_PMS(
          type = "Alcohol",
          qcmetrics = qcmetrics,
          samplesheet = samplesheet,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )
        })

      output$BMI <- renderPlotly({
        shiny_plot_PMS(
          type = "BMI",
          qcmetrics = qcmetrics,
          samplesheet = samplesheet,
          define_outliers = input$define_outliers,
          webGL = values$webGL,
          plotdata = plotdata(),
          colorbyQCmetric = input$colorbyQCmetric,
          scatterColorByColumn = input$scatterColorByColumn,
          selected_samples = values$selected_samples,
          outliers = outliers(),
          show_outliers = input$show_outliers,
          colorscheme = values$colorscheme
        )
      })

      output$overview_outliers <- DT::renderDataTable({
        validate(
          need(input$define_outliers == 1, "The 'define_outliers' option is not selected.")
        )
        outliers_total <- outliers_table() %>%
          tidyr::gather(key = "Metric", value = "Outliers", -c("Sample_Name", "Total")) %>%
          dplyr::group_by(Metric) %>%
          dplyr::summarize(Nr_Outliers = sum(Outliers)) %>%
          dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier$", ""))
        outliers_unique <- outliers_table_plot() %>%
          tidyr::gather(key = "Metric", value = "Outliers", dplyr::ends_with("unique")) %>%
          dplyr::group_by(Metric) %>%
          dplyr::summarize(Nr_Unique_Outliers = sum(Outliers)) %>%
          dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier_unique$|_unique$", ""))
        outliers_total <- outliers_total %>%
          dplyr::left_join(outliers_unique[,c("Metric", "Nr_Unique_Outliers")], by = "Metric")
        outliers_total

      })

      output$outliers <- DT::renderDataTable({
        validate(
          need(input$define_outliers == 1, "The 'define_outliers' option is not selected.")
        )
        outliers_table() %>% dplyr::filter(Total)
      })

      output$download_overview <- downloadHandler(
        filename = function() {
          if(is.null(qcmetrics@metadata$name)) {
            "sample_outlier_overview.rds"
          } else {
            sprintf("sample_outlier_overview_%s.rds", qcmetrics@metadata$name)
          }
        },
        content = function(file) {

          outliers_total <- outliers_table() %>%
            tidyr::gather(key = "Metric", value = "Outliers", -c("Sample_Name", "Total")) %>%
            dplyr::group_by(Metric) %>%
            dplyr::summarize(Nr_Outliers = sum(Outliers)) %>%
            dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier$", ""))
          outliers_unique <- outliers_table_plot() %>%
            tidyr::gather(key = "Metric", value = "Outliers", dplyr::ends_with("unique")) %>%
            dplyr::group_by(Metric) %>%
            dplyr::summarize(Nr_Unique_Outliers = sum(Outliers)) %>%
            dplyr::mutate(Metric = stringr::str_replace(Metric, "_outlier_unique$|_unique$", ""))
          outliers_total <- outliers_total %>%
            dplyr::left_join(outliers_unique[,c("Metric", "Nr_Unique_Outliers")], by = "Metric")

          saveRDS(outliers_total, file = file)
        }
      )

      output$download_overview_button <- renderUI({
        if (input$define_outliers == 1) {
          downloadButton("download_overview", "Download overview")
        }
      })

      ### Probe QC
      output$detP_probes <- renderPlot({
        if(input$define_outliers_probes == 1) {

            g <- qcmetrics@Probe_Metrics %>% dplyr::filter(detP_Percentage > 0.01) %>% ggplot(aes(x = Probe, y = detP_Percentage)) +
              geom_point(alpha=0.6) + theme(
                axis.text.x = element_blank(),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
              geom_hline(yintercept = input$detP_probes, linetype = "dashed", color = "red") +
              ylab(sprintf("Proportion of samples with detP > %s", detp_threshold)) +
              xlab("Probes")
            if(input$show_outliers_probes == 1) {
               g <- g + geom_point(data = qcmetrics@Probe_Metrics %>% dplyr::filter(Probe %in% probe_outliers()), aes(x = Probe, y = detP_Percentage), shape = 8, color = "red")
              }
            g

      } else {
             g <- qcmetrics@Probe_Metrics %>% dplyr::filter(detP_Percentage > 0.01) %>% ggplot(aes(x = Probe, y = detP_Percentage)) +
              geom_point(alpha=0.6) + theme(
                axis.text.x = element_blank(),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
                ylab(sprintf("Proportion of samples with detP > %s", detp_threshold)) +
              xlab("Probes")
            g
         }
      })

      output$beadNr_probes <- renderPlot({
        if(input$define_outliers_probes == 1) {
            g <- qcmetrics@Probe_Metrics %>% dplyr::filter(BeadNr_Percentage > 0.01) %>% ggplot(aes(x = Probe, y = BeadNr_Percentage)) +
              geom_point(alpha=0.6) + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
              geom_hline(yintercept = input$beadNr_probes, linetype = "dashed", color = "red") +
              ylab(sprintf("Proportion of samples with beadNr < %s", beadnr_threshold)) +
              xlab("Probes")
            if(input$show_outliers_probes == 1) {
               g <- g + geom_point(data = qcmetrics@Probe_Metrics %>% dplyr::filter(Probe %in% probe_outliers()), aes(x = Probe, y = BeadNr_Percentage), shape = 8, color = "red")
              }
            g

      } else {
             g <- qcmetrics@Probe_Metrics %>% dplyr::filter(BeadNr_Percentage > 0.01) %>% ggplot(aes(x = Probe, y = BeadNr_Percentage)) +
              geom_point(alpha=0.6) + theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12)) +
                ylab(sprintf("Proportion of samples with beadNr < %s", beadnr_threshold)) +
              xlab("Probes")
            g
         }
      })

    },
    options = list(launch.browser = TRUE)
    )
}

## Additional theme (see issue: https://github.com/tidyverse/ggplot2/issues/1555)
theme_minimal2 <- function(base_size = 12, base_family = "") {
  # Starts with theme_bw and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      strip.background  = element_blank(),
      plot.background   = element_blank(),
      # panel.grid.major.y = element_blank(),
      axis.ticks        = element_line(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.ticks.length = unit(0, "lines")
    )
}

## Color schemes ----------------------------------------------------------------
# adapted from: https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
scale_color_OkabeIto <- function() {

  values <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  n <- length(values)
  pal <- function(n) {
    if (n > length(values)) {
      warning("Insufficient values in manual scale. ", n, " needed but only ",
              length(values), " provided.", call. = FALSE)
    }
    values
  }
  ggplot2::discrete_scale("colour", "manual", pal)
}

set_color_scheme <- function(plotdata,column, colorscheme) {
  if(is.factor(plotdata[[column]]) || is.character(plotdata[[column]]) || is.logical(plotdata[[column]])) {
    type <- "d"
  } else {
    type <- "c"
  }
  if(colorscheme == "viridis") {
    if(type == "d") return(scale_color_viridis_d()) else return(scale_color_viridis_c())
  }
  if (colorscheme == "Okabe-Ito") {
    if(type == "d") return(scale_color_OkabeIto()) else return(scale_color_viridis_c())
  }
  if(colorscheme == "ggplot") {
    if(type == "d") return(scale_color_discrete()) else return(scale_color_continuous())
  }
}
