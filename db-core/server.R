

server <- function(input, output,session = session) {
    #HOME-------
    output$HOME_output_text <- renderUI({
        div(p("In DRCTdb, we collected and processed overall 2.6 million cells with transcriptome and epigenetics information from 16 studies. We also integrated GWAS data of 48 genetic diseases with single-cell multiomics data to identify disease related cell types. We will continuously enhance the DRCTdb with the advancements in single-cell multiomics. As new single-cell multiomics data becomes available, we will regularly update and upgrade the DRCTdb to ensure its relevance and comprehensiveness."),
            style ='display: inline;'
        )
    })
    
    #Contact-------
    output$Contact_text <- renderUI({
        div(
            h3('About'),
            p('We hope you find this data resource useful. Please contact us with your experiences and suggestions'),
            br(),
            
            h3('Contact'),
            p("If your have any question, please don't hesitate to contact us."),
            p("Yunhui: kongyunhui1@gmail.com"),
            p("Junyao:jyjiang@link.cuhk.edu.hk"),
            br(),
            h3('Citation'),
            p("xxx"),
            style ='display: inline;'
        )
    })
    #search-----
    output$coretable = renderDataTable(
        coretable[,c(1,4:6,9)],
        selection = 'single',
        rownames =FALSE,
        server = TRUE,
        options = list(
            autoWidth = FALSE,
            pageLength = 10,
            searchHighlight = TRUE,
            lengthChange = FALSE)
    )
    #unique(coretable$Dataset)
    #core_table_rows_selected
    output$download_table = renderDataTable(
        download_table[, c("Study name","dataset", "Sample", "Processed scRNA", "Processed scATAC",'PMID')],
        rownames =FALSE,
        selection = 'none',
        server = TRUE,
        escape = FALSE,
        options = list(
            autoWidth = FALSE,
            pageLength = 16,
            searching = FALSE,
            lengthChange = FALSE)
    )
    
    
    output$seleted_row = renderPrint({
        s = input$coretable_rows_selected
        if (length(s)) {
            cat('These rows were selected:\n\n')
            cat(coretable[s,]$dataset, sep = ', ')
        }
    })
    observeEvent(input$go_to_panel, {
        updateTabsetPanel(session,"inTabset", selected =  "Results")
    })

    #Results-----------
    ##left side-----
    Select_dataseted <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$dataset
        }else{
            Select_dataset = 'dataset1'
        }
    })

    # output$test <- renderText({
    #     Select_dataseted()
    # })
    output$study_name <- renderText({
        s = input$coretable_rows_selected
        if (length(s)) {
            study_name =  coretable[s,]$`Study name`
        }else{
            study_name = 'Hocker et al. (Sci Adv, 2021)'
        }
        return(study_name)
    })
    
    
    output$show_umap <- renderImage({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$dataset
        }else{
            Select_dataset = 'dataset1'
        }
        list(
            src = file.path("SVG/", paste0(Select_dataset, ".svg")),
            width = 600
        )
    }, deleteFile = FALSE)
    
    
    drcttable <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        #file_dir <- 'downstream_result/sample1/sample1_related_disease.xls'
        file_dir <- list.files(paste0('downstream_result/',Select_dataset),pattern = 'related_disease.xls',full.names = T)
        df <- data.table::fread(file_dir)
        df$`p value` <- signif(df$`p value`,digits =3)
        # df$`p value` <-  sapply(df$`p value`, function(x) {
        #     if (nchar(strsplit(format(x), split = "\\.")[[1]][2]) > 4) {
        #         format(x, scientific = TRUE)
        #     } else {
        #         format(x, nsmall = 4)
        #     }
        # })
        return(df)
    })
    output$show_DRCT_table <- renderDT({
        drcttable()
    }
    )
    der_table <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        file_dir <- list.files(paste0('downstream_result/',Select_dataset),pattern = 'DERs.txt.gz',full.names = T)
        df <- data.table::fread(file_dir)
        df <- df[,c(1:3,6,10,11)]
        colnames(df) <- c('chr','start','end','Cell type','Closed gene','region')
        return(df)
    })
    output$show_DER_table <- renderDT({
        der_table()
    }
    )
    
    tf_activity <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        file_dir <- list.files(paste0('downstream_result/',Select_dataset),pattern = 'tf_activity',full.names = T)
        df <- data.table::fread(file_dir)
        colnames(df) <- c('Cell type','TF')
        return(df)
    })
    output$show_tf_activity<- renderDT({
        tf_activity()
    }
    )
    
    
    ##right side-----
    output$dropdown <- renderUI({
        fluidRow(
            selectInput(inputId = 'ct',
                        label = 'Choose Cell type',
                        choices = unique(drcttable()[[2]]),
                        selected = unique(drcttable()[[2]])[1],
                        width = '300px')
        )
    })
    
   
    output$diseaseDropdown <- renderUI({
        if (is.null(input$ct)) return(NULL)  
        selectInput(inputId = 'disease',
                    label = 'Choose disease',
                    choices = unique(drcttable()[drcttable()[[2]] == input$ct, ][[1]]),
                    selected = unique(drcttable()[drcttable()[[2]] == input$ct, ][[1]])[1],
                    width = '300px')
    })
    
    ###ccc
    output$show_ccc_plot <- renderImage({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        ccc_path = list.files(path = paste0('downstream_result/',Select_dataset,'/ccc/'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset('svg')
        if (length(ccc_path) == 0) {
          ccc_path = 'www/CCC_error.png'
        }
        list(
            src = ccc_path,
            width = 600
        )
    }, deleteFile = FALSE)
    
    ccc_table <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        ccc_tabe_path = list.files(path = paste0('downstream_result/',Select_dataset,'/ccc/'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset('txt')
        if (length(ccc_tabe_path) == 0) {
          ccc_tabe_path = 'www/ccc_error.txt'
        }
        df <- data.table::fread(ccc_tabe_path)
        df <- df[,c(12,7,10,11,13)]
        df$prob.original <- signif(df$prob.original,digits =3)
        return(df)
    })
    output$show_ccc_table <- renderDT({
        #ccc_table()
        datatable(
            ccc_table(),
            extensions = 'Buttons',
            options = list(
                dom = 'Bfrtip',
                buttons = list(
                    list(
                        extend = 'excel',
                        text = 'Download to Excel',
                        filename = 'CCC_table',
                        title = '', 
                        messageTop = '',  
                        messageBottom = '',
                        exportOptions = list(
                            modifier = list(page = "all")
                        )
                    )
                )
            )
        )
    },server = FALSE
    )
    output$show_ccc <-  renderUI({
        if (input$switchccc) {
        imageOutput("show_ccc_plot")
        }else{
            DT::dataTableOutput('show_ccc_table',width = "100%")
        }
        
    })
    
    ###atac
    
    atac_table <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        atac_path = list.files(path = paste0('downstream_result/',Select_dataset,'/atac_snp'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset(gsub(' ','_',input$ct)) %>%
            stringr::str_subset('txt')
        
        df <- data.table::fread(atac_path)
        df <- df[,c(1,6,7,9,10,11)]
        return(df)
    })
    
    output$atac_table <- renderDT({
        #atac_table()
        datatable(
            atac_table(),
            extensions = 'Buttons',
            options = list(
                dom = 'Bfrtip',
                buttons = list(
                    list(
                        extend = 'excel',
                        text = 'Download to Excel',
                        filename = 'ATAC_table',
                        title = '', 
                        messageTop = '',  
                        messageBottom = '',
                        exportOptions = list(
                            modifier = list(page = "all")
                        )
                    )
                )
            )
        )
    },server = FALSE)
    
    output$show_atac_plot <- renderImage({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        atac_png_path = list.files(path = paste0('downstream_result/',Select_dataset,'/atac_snp'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset(gsub(' ','_',input$ct)) %>%
            stringr::str_subset('png')
        
        list(
            src = atac_png_path,
            width = 800
        )
    }, deleteFile = FALSE)
    
    output$show_atac_enrich <-  renderUI({
        if (input$switchatac) {
            imageOutput("show_atac_plot")
        }else{
            DT::dataTableOutput('atac_table',width = "100%")
        }
        
    })
    
    
    
    ###rna
    rna_table <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        rna_path = list.files(path = paste0('downstream_result/',Select_dataset,'/rna_snp'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset(gsub(' ','_',input$ct)) %>%
            stringr::str_subset('txt')
        
        df <- data.table::fread(rna_path)
        df <- df[,c(1,2,3,8,9,10)]
        return(df)
    })
    output$rna_table <- renderDT({
        #rna_table()
        
        datatable(
            rna_table(),
            extensions = 'Buttons',
            options = list(
                dom = 'Bfrtip',
                buttons = list(
                    list(
                        extend = 'excel',
                        text = 'Download to Excel',
                        filename = 'rna_table',
                        title = '', 
                        messageTop = '',  
                        messageBottom = '',
                        exportOptions = list(
                            modifier = list(page = "all")
                        )
                    )
                )
            )
        )
    },server=FALSE)
    
    output$show_rna_plot <- renderImage({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        rna_png_path = list.files(path = paste0('downstream_result/',Select_dataset,'/rna_snp'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset(gsub(' ','_',input$ct)) %>%
            stringr::str_subset('png')
        
        list(
            src = rna_png_path,
            width = 800
        )
    }, deleteFile = FALSE)
    
    output$show_rna_enrich <-  renderUI({
        if (input$switchrna) {
            imageOutput("show_rna_plot")
        }else{
            DT::dataTableOutput('rna_table',width = "100%")
        }
        
        
        
    })
    
    
    
    
    
    
    ##grn
    grn_file <- reactive({
        s = input$coretable_rows_selected
        if (length(s)) {
            Select_dataset =  coretable[s,]$Sample
        }else{
            Select_dataset = 'sample1'
        }
        
        grn_path = list.files(path = paste0('downstream_result/',Select_dataset,'/grn_cor02'),pattern = gsub(' ','_',input$disease),full.names = TRUE) %>% 
            stringr::str_subset(gsub(' ','_',input$ct)) %>%
            stringr::str_subset('Rds')
        grn <- readRDS(grn_path)
        return(grn)
    })
    output$grn_output <-  networkD3::renderForceNetwork({
        grn1 <- networkD3::forceNetwork(Links = grn_file()$link, Nodes = grn_file()$nodes,
                     Source = "source", Target = "target",fontSize = 10,
                     Value = "value", NodeID = "gene",zoom = TRUE,arrows =TRUE,
                     Group = "group2", opacity = 0.8)
        htmlwidgets::onRender(grn1, "
          function(el, x) {
            // 选择所有节点并为每个节点添加文本标签
            d3.select(el).selectAll('.node').append('text')
              .attr('dx', 12)
              .attr('dy', '.35em')
              .text(function(d) { return d.name; }) // 确保这里使用正确的属性来显示名称
              .style('font-size', '10px')
              .attr('class', 'node-label');
          }
        ")
    })
    
    
    ###tools part 
    enrich_data <- reactive({
        req(input$enrich_upload)
        ext <- tools::file_ext(input$enrich_upload$name)
        #print(input$enrich_upload$datapath)
        enriched_data <- switch(ext,
               xls = readxl::read_xls(input$enrich_upload$datapath),
               xlsx = readxl::read_xlsx(input$enrich_upload$datapath)
               #validate("Invalid file; Please upload a .xls or .xlsx file")
        )
        colnames(enriched_data) <- toupper(colnames(enriched_data))
        return(enriched_data)
    })
    
    enrich_res <- reactive({
        req(input$enrich_upload)
        #enrichres()
        if ("SYMBOL"  %in% toupper(colnames(enrich_data()))){
            enrich_GO(enrich_data())
        }else{
            stop('Must have SYMBOL columns')
        }
    })
    output$enrichres_table <- renderDataTable({
        req(input$enrich_upload)
        datatable(enrich_res()[[1]])
    })
    
    
    output$download_enrich_data <- downloadHandler(
        filename = function() {
            paste("enriched_reasuts", Sys.Date(), ".xlsx", sep = "")
        },
        content = function(file) {
            writexl::write_xlsx(enrich_res(), file)
        }
    )

    
}
