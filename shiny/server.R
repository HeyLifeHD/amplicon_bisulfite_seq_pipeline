source("utilities.R")


server <- function(input,output,session){
  status <- reactiveVal()
  

    #############MassArray###################


# Load sample data --------------------------------------------------------

    observeEvent(input$MA_sampleData,{
      updateTextInput(session,"textRegion",value="chr19:43203328-43203689")
      updateSelectInput(session,"selectGenome",selected = "hg19")
      #shinyjs::click("SubmitBTNRegion")
    })


    #####general ouput and variables########


    #repeating the input
    inReg <- eventReactive(input$SubmitBTNRegion,{
      status("region")
      input$textRegion
    })



    #determine chromosome, begin and end of region on chromosome from the input region
    region <- eventReactive(input$SubmitBTNRegion,{
      req(inReg())
      gsub("[,]","",inReg())
    })

    chromosome <- reactive({
      req(inReg())
      stringr::str_split_fixed(inReg(),":",2)[1]
    })
    beginFirst <-reactive({
      req(inReg())
      as.numeric(substr(region(),(regexpr(':',region())[1]+1),(regexpr("-",region())[1]-1)))
    })

    endFirst <- reactive({
      req(inReg())
      as.numeric(substr(region(),(regexpr('-',region())[1]+1),nchar(region())))

    })
    observeEvent(input$SliderPlot,{
      status("slider")
    })


# Slider ------------------------------------------------------------------
    # show info text for slider
    output$sliderInfo <- renderText({
      req(inReg())
      "Slider to change the size of the input region."
    })

    #show sliderInput for zooming options
    output$coordinateSlider <- renderUI({
      sliderInput("SliderPlot","Choose range for zoom: ",min = max(1,beginFirst()-10000), max = endFirst()+10000, value= c(beginFirst(),endFirst()),width="1800px",step=1)
    })


# Specify input -----------------------------------------------------------

    begin <- reactive({
      st <- status()
      if(st == "slider" && !is.null(input$SliderPlot[1]) && length(input$SliderPlot)>0){
        input$SliderPlot[1]
      }else{
        beginFirst()
      }
    })

    end <- reactive({
      st <- status()
      if(st == "slider" && !is.null(input$SliderPlot[2]) && length(input$SliderPlot)>0){
        input$SliderPlot[2]
      }else{
        endFirst()
      }
    })

    output$inputMA <- renderTable({
      req(inReg())
      data.frame(
        "Genomic_Location"=c("Selected Genome","Chromosome","Start","End"),
        "Input"=c(input$selectGenome,chromosome(),begin(),end())
      )}
    )

    output$inputRegion <- renderText({
      req(inReg())
      if(end()-begin()+1 >2000){
        "Chosen region is bigger than 2kb"
      }
    })


# Download Button ---------------------------------------------------------

    #downloads all the information that is also shown in the shiny app
    output$downloadB <- renderUI({
      req(region())
      downloadButton("downloadData", "Download information")
    })



# Export Sequence ---------------------------------------------------------

    #upon the region entered look for the related sequence in UCSC
    fwSeq <- reactive({

      if(begin() <= 0 || begin() > end() || !(chromosome() %in% paste0("chr",c(as.character(c(1:22)),"X","Y")))){
        showNotification("Enter region like: chr19:43203328-43203389.",type="default")
        shiny::validate(
          need(begin() > 0, "The begin position of the entered region must be greater than zero"),
          need(begin() <= end(),"The end position of the entered region has to be after the begin position"),
          need(chromosome() %in% paste0("chr",c(as.character(c(1:22)),"X","Y")), "The chromosome must be one of 1-22 or X, Y. Check that you have entered the region in the correct format.")
        )
      }else{
        discoverSequence(chromosome(),begin(),end(),input$selectGenome)
      }
    })



    #textInput for the user to decide the letters per line
    output$selectlinebreak <- renderUI({
      req(fwSeq())
      textInput("numChar","Letters per line: (1-100)",value="80",width = "150px",placeholder ="80")
    })
    #build reverse complement of forward sequence
    revCompl <- reactive({
      req(fwSeq())
      makeReverseComplement(fwSeq())
    })

    #delete linebreaks from forward sequence
    FWseqAsString <- reactive({
      gsub("[\r\n]","",as.character(fwSeq()))
    })


    #delete linebreaks from reverse complement
    RCseqAsString <- reactive({
      gsub("[\r\n]","",as.character(revCompl()))
    })

    #print texts for the sequences
    output$textfw <- renderText({
      req(fwSeq())
      "The forward sequence of the region you entered is:"
    })

    output$textColoring <- renderText({
      req(fwSeq())
      HTML("Coloring: &emsp; <b><span style='color:green'>Repeats</span> &emsp; <span style='color:red'>CpG</span> &emsp; <span style='color:violet'>SNP</span></b>")
    })

    output$textRC <- renderText({
      req(fwSeq())
      "The reverse complement of the sequence is:"
    })

    #count number of CGs in FW- and RC-sequences and determine their positions
    CGsinFWseq <- reactive({
      x <- findCG(FWseqAsString())[[1]]
      x <- apply(x,c(1,2),function(y) as.integer(y-1+begin()))
      x
    })
    CGsinRC <- reactive({
      x <- findCG(RCseqAsString())[[1]]
      x <- apply(x,c(1,2),function(y) as.integer(y+begin()))
      x
    })

    output$CGFWtext <- renderText({
      paste0("There are ", nrow(CGsinFWseq()), " CGs in the forward sequence at the following positions: ")
    })
    output$CGinFWseq <- renderTable(
      CGsinFWseq(),rownames=TRUE
    )
    output$CGinRCtext <- renderText({
      paste0("There are ", nrow(CGsinRC())," CGs in the reverse complement at the following positions: ")
    })
    output$CGinRC <- renderTable(
      CGsinRC(),rownames =TRUE
    )




# Create the plot ---------------------------------------------------------

    output$PlotOutputs <- renderPlot({
      gviz_primer(input$selectGenome,chromosome(),begin(),end(),strand="+",primer=F)
    })


# HTML output of sequence -------------------------------------------------

    # Map SNPs
    SNPpositions <- reactive({
      dbQuery(chromosome(),begin(),end(),input$selectGenome,S=T)
    })

    # Map repeats
    repeats <- reactive({
      dbQuery(chromosome(),begin(),end(),input$selectGenome,S=F,R=T)
    })


    #change the lineLength according to user input
    observeEvent(input$numChar,{
      req(input$numChar)
      if(is.null(input$numChar)|| as.numeric(input$numChar) < 1 || as.numeric(input$numChar) > 100)
        showNotification("Choose a number between 1 and 100 (default is 80)",type="warning")
    })

    lineLength <- reactive({
      req(input$numChar)
      if(is.null(input$numChar) || as.numeric(input$numChar)<1 || as.numeric(input$numChar) > 100){
        80
      }else{
        as.numeric(input$numChar)
      }}
    )

    #save the start positions of the SNPs into a variable
    vecSNPstart <- reactive({
      x <- c()
      #if there were  SNPs loaded, add the highlighting
      if(!is.null(SNPpositions())){
        #if there are no SNPs do not color
        if (length(start(SNPpositions()))>0){
          for(i in 1:length(start(SNPpositions()))){
            x <- append(x,start(SNPpositions())[i])
          }
        }else{NULL}
        x
      }else{NULL}
    })

    #save the widths of the SNPs into variables
    vecSNPwidth <- reactive({
      x <- c()
      #if there were  SNPs loaded, add the highlighting
      if(!is.null(SNPpositions())){
        #if there are no SNPs do not color
        if (length(width(SNPpositions()))>0){
          for(i in 1:length(width(SNPpositions()))){
            x <- append(x,width(SNPpositions())[i])
          }
        }else{NULL}
        x
      }else{NULL}

    })

    #calculate the SNP positons for the backward sequence on the reverse strand
    vecSNPstartBack <- reactive({
      y <- c()
      if(!is.null(vecSNPstart()) && !is.null(vecSNPwidth()) && length(vecSNPwidth() > 0) && length(vecSNPstart() >0)){
        for(i in 1:length(vecSNPstart())){
          #calculate the SNP position on the reverse strand
          y <- append(y,(end()-vecSNPstart()[i]+begin()))
        }
        y
      }else{
        NULL
      }
    })
    #
    #save the repeat starting positions into vector
    #
    vecRepStart <- reactive({
      rbegin <- c()
      if(!is.null(repeats()) && length(start(repeats()))>0){
        for(x in 1:length(start(repeats()))){
          rbegin <- append(rbegin,start(repeats()[x]))
        }
        rbegin
      }else{NULL}
    })


    #
    #save the widths of the repeats into vector
    #
    vecRepWidth <- reactive({
      rwidth <- c()
      if(!is.null(repeats()) && length(width(repeats()))>0){
        for(x in 1:length(width(repeats()))){
          rwidth <- append(rwidth,width(repeats()[x]))
        }
        rwidth
      }else{NULL}
    })

    #calculate the repeat starting positions on the reverse complement strand
    vecRepStartBack <- reactive({
      y <- c()
      if(!is.null(vecRepStart()) && !is.null(vecRepWidth()) && length(vecRepStart() > 0) && length(vecRepWidth() >0)){
        for(i in 1:length(vecRepStart())){
          #calculate the SNP position on the reverse strand
          y <- append(y,(end()-vecRepStart()[i]+begin()-vecRepWidth()[i]+1))
        }
        y
      }else{
        NULL
      }
    })

    #first color the SNPs,CpGs and repeats of the fw sequence and then print it
    fwseqHTML <- reactive({
      HTML(paste0("<pre style='word-wrap: break-word;width:100%;margin-left:5%'><font size ='4'><p align='justify'>",doColoring(begin(),FWseqAsString(),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,NULL,NULL),"</p><font size></pre>"))
      print(HTML(paste0("<pre style='word-wrap: break-word;width:100%;margin-left:5%'><font size ='4'><p align='justify'>",doColoring(begin(),FWseqAsString(),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,NULL,NULL),"</p><font size></pre>")))
    })

    #generate a string that shows the numbering of lines for forw and backw sequence
    stringNumbering <- reactive({
      req(FWseqAsString())
      req(input$numChar)

      str <-"1\n"
      x <- 1
      if(!is.null(FWseqAsString()) && length(FWseqAsString) > 0 && !is.null(input$numChar)&& length(input$numChar)>0 ){
        while(x+as.numeric(input$numChar) <= nchar(FWseqAsString())){
          str <- paste0(str,x+as.numeric(input$numChar))
          str <- paste0(str,"\n")
          x <- x + as.numeric(input$numChar)
        }
        str
      }
    })

    #print the numbering of the shown forward and backward sequence
    output$letterNumberForw <- renderText({
      HTML(paste0("<pre style='word-wrap: break-word;width:100%;margin-left:5%'><font size ='4'>",stringNumbering(),"<font size></pre>"))
    })
    output$letterNumberBack <- renderText({
      HTML(paste0("<pre style='word-wrap: break-word;width:100%;margin-left:5%'><font size ='4'>",stringNumbering(),"<font size></pre>"))
    })
    #print the sequence with all the colored features(SNPs,CGs and repeats)
    output$forwardSeq <- renderText(
      fwseqHTML()
    )

    #color only CpGs of rc sequence and print this
    rcseqHTML <- reactive({
      HTML(paste0("<pre style='word-wrap: break-word;width:100%;margin-left:5%'><font size='4'>", doColoring(begin(),RCseqAsString(),vecSNPstartBack(),vecSNPwidth(),vecRepStartBack(),vecRepWidth(),80,NULL,NULL),"<font size></pre>"))

    })
    output$revC <- renderText({
      rcseqHTML()
    })



# Primer Design -----------------------------------------------------------

    #show slider input for choosing the primer sequences for forward and backward primer
    output$slider1 <- renderUI({
      if(!is.null(input$SliderPlot))
        sliderInput("SliderForw","Choose range for forward primer:",min = round(as.numeric(input$SliderPlot[1]),digits =0),max = round(as.numeric(input$SliderPlot[2]),digits =0),step =1,value = c(round(as.numeric(input$SliderPlot[1]),digits =0),round(as.numeric(input$SliderPlot[2])),digits =0),width = "700px")
    })
    output$slider2 <- renderUI({
      if(!is.null(input$SliderPlot))
        sliderInput("SliderBack","Choose range for backward primer:",min = round(as.numeric(input$SliderPlot[1]),digits =0),max = round(as.numeric(input$SliderPlot[2]),digits =0),step =1,value = c(round(as.numeric(input$SliderPlot[1]),digits =0),round(as.numeric(input$SliderPlot[2])),digits =0),width = "700px")
    })

    strandSel <- reactive({
      if(input$standSelector=="Watson"){
        FWseqAsString()
      }else{
        makeComplement(FWseqAsString())
      }
    })



    #mark the chosen primer in the zoomed in sequence
    htmlOutBackward <- reactive({
      tmp <- ""
      if(input$standSelector=="Watson"){
          if(!is.null(input$SliderBack) && length(input$SliderBack) > 0){
            tmp <- HTML(paste0("The bisulfit converted complement of the chosen zoom range (",input$SliderPlot[1],"-",input$SliderPlot[2],") is \n <br> <pre style='word-wrap: break-word;width:80%;margin-left:5%'>3' ",doColoring_rev(begin(),makeComplement(bisulfitConvert(strandSel())),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)," 5'</pre>"))
          }}else{
            if(!is.null(input$SliderForw) && length(input$SliderForw) > 0){#TODO mit repeat gut??
              tmp <- HTML(paste0("The bisulfit converted Crick strand sequence (3'-5') of the chosen zoom range (",input$SliderPlot[1],"-",input$SliderPlot[2],")is\n <br> <pre> 3'",doColoring_rev(begin(),bisulfitConvert_rev(strandSel()),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)," 5'</pre>"))
          }
        }
      tmp
    })
    output$seqTextBackward <- renderUI({
      htmlOutBackward()
    })

    #mark the chosen primer in the zoomed in forward sequence
    htmlOutForward <- reactive({
      if(input$standSelector=="Watson"){
      if(!is.null(input$SliderForw) && length(input$SliderForw) > 0){#TODO mit repeat gut??
        HTML(paste0("The bisulfit converted forward sequence of the chosen zoom range (",input$SliderPlot[1],"-",input$SliderPlot[2],")is\n <br> <pre> 5'",doColoring(begin(),bisulfitConvert(strandSel()),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1)," 3'</pre>"))
      }}else{
        if(!is.null(input$SliderBack) && length(input$SliderBack) > 0){
          HTML(paste0("The bisulfit converted Crick strand complement (5'-3') of the chosen zoom range (",input$SliderPlot[1],"-",input$SliderPlot[2],") is \n <br> <pre style='word-wrap: break-word;width:80%;margin-left:5%'>5' ",doColoring(begin(),makeComplement(bisulfitConvert_rev(strandSel())),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),80,input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1)," 3'</pre>"))
        }
      }
    })
    output$seqTextForward <- renderUI({
      htmlOutForward()
    })
    output$suggestionPrimer <- renderText({
      req(fwSeq())
      HTML(paste0("<pre><font size ='3'><b>Instruction for MassArray primers:</b>\n\n","- Tm >= 52°C, optimal ~ 60°C\n\n","- should not contain CpG or SNP; if \n  unavoidable, then only within first\n  quarter(5’ to 3’) of primer\n\n",
                  "- after primer design, forward primer must\n  be preceded at 5’-end by AGGAAGAGAG and \n  reverse primer at 5’-end by CAGTAATACGACTCACTATAGGGAGAAGCT before ordering</font size></pre>"))
    })

    output$suggestionAmplicon <- renderText({
      req(fwSeq())
      HTML(paste0("<pre><font size ='3'><b>Instruction for MassArray amplicons:</b>\n\n","- must contain at least 1 CpG\n\n- if DNA from fresh or frozen \n  specimens: size between 100-500 bp\n\n- if DNA from FFPE specimens: size \n  between 100-200 bp</font size></pre>"))
    })
    #calculate the melting tmeperature of the chosen forward primer
    #print the number of letters per line
    output$letterNumberBP <- renderText({
      HTML(paste0("<pre><font size ='1'>\n",stringNumbering(),"<font size></pre>"))
    })
    output$letterNumberFP <- renderText({
      HTML(paste0("<pre><font size ='1'>\n",stringNumbering(),"<font size></pre>"))

    })


    #calculate the melting tmeperature of the chosen forward primer
    calcTempForw <- reactive({
      if(input$standSelector=="Watson"){
      req(input$SliderForw)
      primerForw <- bisulfitConvert(substr(FWseqAsString(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1))
      paste0("The melting temperature of the above primer sequence is:   ",calcMeltTemp(primerForw))}else{
        req(input$SliderForw)
        primerForw <-  makeComplement(bisulfitConvert(substr(strandSel(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1)))
        paste0("The melting temperature of the above primer sequence is:   ",calcMeltTemp(primerForw))
      }
    })
    #calculate the melting temp of the chosen backwards primer
    calcTempBack <- reactive({
      if(input$standSelector=="Watson"){
      req(input$SliderBack)
          p <- bisulfitConvert(makeComplement(substr(FWseqAsString(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)))
          paste0("The melting temperature of the above primer sequence is:   ",calcMeltTemp(p))}else{
            req(input$SliderBack)
            p <- Biostrings::reverse(bisulfitConvert(substr(strandSel(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)))
            paste0("The melting temperature of the above primer sequence is:   ",calcMeltTemp(p))
          }
    })

    #print the melting temperatures
    output$temperatureForw <- renderText({
      calcTempForw()
    })
    output$temperatureBack <- renderText({
      calcTempBack()
    })


    # create a table to show the primer output
    primer.df <- reactive({
      req(FWseqAsString())
      req(input$SliderForw)
      req(input$SliderBack)
      req(inReg())
      if(input$standSelector=="Watson"){
      data.frame("Forward_Primer"=c(substr(bisulfitConvert(substr(FWseqAsString(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1)),1,35),
                                    paste0(stringr::str_split_fixed(inReg(),":",2)[,1],":",input$SliderForw[1],"-",input$SliderForw[2]),
                                    as.numeric(input$SliderForw[2]-input$SliderForw[1]),
                                    stringr::str_remove(calcTempForw(),"The melting temperature of the above primer sequence is:")),
                 "RevComp_Primer"=c(Biostrings::reverse(substr(makeComplement(bisulfitConvert(substr(FWseqAsString(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1))),1,35)),
                                        paste0(stringr::str_split_fixed(inReg(),":",2)[,1],":",input$SliderBack[1],"-",input$SliderBack[2]),
                                        as.numeric(input$SliderBack[2]-input$SliderBack[1]),
                                        stringr::str_remove(calcTempBack(),"The melting temperature of the above primer sequence is:"))
      ,row.names = c("Sequence[max.35 Bases shown]","Coordinates","Primer_Size","Melting_Temperature"))}else{
        data.frame(
          "Primer_Forward"=c(substr(makeComplement(bisulfitConvert(substr(strandSel(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1))),1,35),
                             paste0(stringr::str_split_fixed(inReg(),":",2)[,1],":",input$SliderForw[1],"-",input$SliderForw[2]),
                             as.numeric(input$SliderForw[2]-input$SliderForw[1]),
                             stringr::str_remove(calcTempForw(),"The melting temperature of the above primer sequence is:")),
          "Primer_Reverse"=c(Biostrings::reverse(substr(bisulfitConvert(substr(strandSel(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)),1,35)),
                                      paste0(stringr::str_split_fixed(inReg(),":",2)[,1],":",input$SliderBack[1],"-",input$SliderBack[2]),
                                      as.numeric(input$SliderBack[2]-input$SliderBack[1]),
                                      stringr::str_remove(calcTempBack(),"The melting temperature of the above primer sequence is:"))
                   ,row.names = c("Sequence[max.35 Bases shown]","Coordinates","Primer_Size","Melting_Temperature"))
      }
    })

    output$ampsize <- renderText({
      req(input$SliderForw)
      req(input$SliderBack)
      paste0("The choosen amplicon has a size of ", input$SliderBack[2]-input$SliderForw[1]," bp.")
    })

    output$primerInfo <- renderTable({
      primer.df()
    },rownames=T)



# MassArray Amplicon prediction Plots -------------------------------------


    #only do amplikonPrediction if the sequence contais CGs
    #the sliderInput gives the possibility to zoom into the ampliconPrediction plot
    prediction <- reactive({
      req(CGsinFWseq())
      if(nrow(CGsinFWseq()) > 0){
        #doAmpliconPrediction(newSeqForw())
        doAmpliconPrediction((fwSeq()))
      }
    })

    #print a legend description for the above plot
    output$desciptionPlot <- renderText({
      req(CGsinFWseq())
      if(nrow(CGsinFWseq()) > 0)
        "The inSilico Assay prediction of the chosen zoom range:"
    })
    #plot the ampliconPrediction with the updated range
    output$amlikonPrediction <- renderPlot({
      prediction()
    })
    output$numberCG <- renderText({
      paste0("The chosen genomic regions contains ", nrow(prediction()$summary)," CpGs.")
    })
    output$assayCG <- renderTable({
      cbind("CpG"=1:nrow(prediction()$summary),prediction()$summary[,-1])
    })

    #legend for the inSilico amplicon Prediction plot
    output$legend <- renderPlot({
      if(nrow(CGsinFWseq()>0)){
        legend("center",legend= c("CG dinucleotides\n(numbered)\n","fragment molecular\nweight overlapping\nwith another fragment\n",
                                  "fragment molecular\nweight outside the\ntestable mass window\n" ,
                                  "fragment containing\na potential\nconversion control",
                                  "fragment uniquely\nassayable but\ncontaining no CGs\n",
                                  "linked arrowheads\ndenote molecular\nweight overlaps\nbetween multiple\nCG-containing fragments",
                                  "yellow highlights\nrepresent tagged or\nprimer sequences"),
               col =c("blue","red","gray","green","black","gray","yellow"),cex = 0.8,lwd =c(5,5,5,5,5,5,5), bty="n",pch = c(16,NA,NA,NA,NA,NA,NA),lty =1,ncol =4)
      }
    })



# Automated Primer Design -------------------------------------------------

    #bisulfit convert the forward sequence
    bs_seq <- reactive({
      bs_convert(FWseqAsString())
    })

    #count number of cg in sequence
    output$cgincluded <- renderUI({
      var <- stringr::str_count(as.character(fwSeq()),"CG")
      checkboxGroupInput("cgincluded","Include CpGs in region", choices = 1:var)})


    # pr_design <- eventReactive(input$MA_primerDesign,{
    #   print(paste0("chr",chromosome(),":",begin(),"-",end()))
    #   my_amplicons <- generate_amplicons_coord(as.data.frame(GenomicRanges::GRanges(
    #     paste0("chr",chromosome(),":",begin(),"-",end()))),
    #     input$selectGenome)
    #
    #
    #   primer_design(my_amplicons,names(my_amplicons),
    #                 primer_size = input$MaPrimS,
    #                 primMaxSize = input$MaMaxPrimS,
    #                 primMinSize = input$MaMinPrimS,
    #                 primTm = input$MaPrimT,
    #                 primMinTm = input$MaMinPrimT,
    #                 primMaxTm = input$MaMaxPrimT,
    #                 amp_min = input$MaAmpMin,
    #                 amp_max = input$MaAmpMax,
    #                 exclude_cpg = input$MaExcludeCG,
    #                 CGtarget = round((as.numeric(end())-as.numeric(begin()))/2))
    # })

    output$textPrimerDesign <- renderText({
      req(pr_design())
      "Primer3 generated the following primers for the input sequence. Select a primer pair for more information."
    })

    #design primer with primer3
    pr_design <- reactive({
        print(paste0(chromosome(),":",begin(),"-",end()))
        print(input$cgincluded)
        my_amplicons <- generate_amplicons_coord(as.data.frame(paste0(chromosome(),":",begin(),"-",end())),
          input$selectGenome)

        #specify the target CpGs
        var <- stringr::str_locate_all(as.character(fwSeq()),"CG")
        targ <- data.frame(var)[as.numeric(input$cgincluded),1]

        if(input$MA_primerDesign){
          primer_design(my_amplicons,names(my_amplicons),
                        primer_size = input$MaPrimS,
                        primMaxSize = input$MaMaxPrimS,
                        primMinSize = input$MaMinPrimS,
                        primTm = input$MaPrimT,
                        primMinTm = input$MaMinPrimT,
                        primMaxTm = input$MaMaxPrimT,
                        amp_min = input$MaAmpMin,
                        amp_max = input$MaAmpMax,
                        exclude_cpg = input$MaExcludeCG,
                        CGtarget = targ
                        )
        }else{
          primer_design(my_amplicons,names(my_amplicons),

                      CGtarget = round((as.numeric(end())-as.numeric(begin()))/2))}
    })




    outputDesignedPrimers <- reactive({
      req(pr_design())
        MA_table <<- cutTable(d=0,NULL,NULL)
    })
    # output$primerDesign <- renderTable({
    #   outputDesignedPrimers()
    # })

    output$primerDesign <- DT::renderDataTable({

      outputDesignedPrimers()

    },selection='single')


# Create the Pop-Up with Primer Infos -------------------------------------

    observeEvent(input$primerDesign_rows_selected, {
      shinyBS::toggleModal(session, "modalMassArray", "open")
    })

    MA_amplicon <- reactive({
      MA_table[input$primerDesign_rows_selected,18]
    })

    MA_output_table <- reactive({
      mass_output(MA_amplicon(),input$selectGenome)
    })

    output$MassArrayPopup <- renderUI({
      shinyBS::bsModal("modalMassArray", h3(strong("Amplicon Report"),align="center"), "", size = "large",
                       # div(conditionalPanel(id='loader',
                       #                      condition="($('html').hasClass('shiny-busy'))",
                       #                      img(src="http://thinkfuture.com/wp-content/uploads/2013/10/loading_spinner.gif")),align="center"),
                       conditionalPanel(id='spin',
                                        condition="($('html').hasClass('shiny-busy'))",
                                        add_busy_spinner(spin = "fading-circle",position="full-page",timeout = 1000)),
                       conditionalPanel(id='MA_user_output',
                                        condition="!($('htmlYes').hasClass('shiny-busy'))",
                                        column(12,align="center",
                                               renderText(paste0("The chose amplicon contains ",MA_output_table()$CpGs," CpGs." )),
                                               br(),br(),
                                               renderText("The following table predicts the MassArray fragments:"),
                                               br(),
                                               renderTable(MA_output_table()$Fragments)))
                       )

    })




    ############ make download of all the data available#################

    #optional download of data into a .csv file
    #the region,DNA sequence,information about Gviz tracks and ampliconPrediction information is printed
    output$downloadData <- downloadHandler(
      filename = "download.html",
      content = function(file) {

        #just produce a nicer output for these parameters
        p <- renderPrint( {
          tryCatch({
            prediction()
          },error= function(e){
            HTML(paste0("~sequence contains no CG nucleotides~"))
          })
        })
        rep <- renderPrint({
            range(repeats())
        })
        snps <- renderPrint({
            HTML(paste0(start(SNPpositions())," - ",start(SNPpositions())+width(SNPpositions()),"\n"))
        })

        #print the chosen forward and backward primer in the file
        primerF <- renderPrint({
          if(input$standSelector=="Watson"){
            bisulfitConvert(substr(FWseqAsString(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1))}else{
              Biostrings::reverse(makeComplement(bisulfitConvert(substr(FWseqAsString(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1))))
            }
        })
        primerB <- renderPrint({
          if(input$standSelector=="Watson"){
            makeComplement(bisulfitConvert(substr(strandSel(),input$SliderBack[1]-begin()+1,input$SliderBack[2]-begin()+1)))}else{
              Biostrings::reverse(bisulfitConvert(substr(strandSel(),input$SliderForw[1]-begin()+1,input$SliderForw[2]-begin()+1)))
            }
        })
        designedPrimers <- renderTable({
            outputDesignedPrimers()
        })
        #print region and FW sequence and CGs in FW sequence
        write.table(HTML(paste0("<pre><font size='4'>The region you entered was ","chr",chromosome(),":",begin()," - ",end()," on the genome ",input$selectGenome,
                                ".\n\nThe corresponding DNA sequence is:\n<font size></pre>", fwseqHTML(),
                                "\n\n<font size='4'>The corresponding reverse complement is:</font size> ",rcseqHTML(),
                                "\n\n<font size='4'>There are <font size>", nrow(CGsinFWseq()),"<font size ='4'>  CGs in the forward sequence at the following positions:</font size><br>"
        )),
        file,row.names=FALSE,col.names = FALSE,append=FALSE,quote =FALSE)
        write.table(CGsinFWseq()[,1],file,row.names=FALSE,col.names = FALSE,append=TRUE,quote = FALSE,sep="<br>",eol = "<br>")
        #print CGs in RC sequence
        #write.table(paste0("\n<font size='4'><br>There are </font size>",nrow(CGsinRC())," <font size='4'>CGs in the reverse complement at the following positions:</font size><br>"),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
        #write.table(CGsinRC()[,1],file,row.names=FALSE,col.names = FALSE,append=TRUE,quote = FALSE, sep="<br>",eol = "<br>")
        #print Gviz tracks into file
        write.table(paste0("\n\n<b><font size ='4'><br><br>Selected information about the entered region:</font size></b>\n\n"),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
        #print additional information if selected

          write.table(HTML(paste0("<b><br>Repeat positions:<br></b>\n\n",rep(),"\n\n<b><br><br>SNP positions:<br></b>\n\n",snps())),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE,sep="<br>",eol="<br>")
        #print information about the chosen forward and backward primer
        #here the new begin position for coloring needs to be calculated: for the backward sequence the amount of base pairs cut off at the end (the sliderPlot) need to be added
        #and the amount of non marked basepairs at the beginning of the backward zoom sequence need to be added to the beginning of the genomic search coordinates
        #write.table(HTML(paste0("<pre>\n\n<b><font size='4'>Primer information:</font size></b>\n\nThe selected forward primer sequence:\n\n",doColoring(input$SliderForw[1],primerF(),vecSNPstart(),vecSNPwidth(),vecRepStart(),vecRepWidth(),60,NULL,NULL), "\n\n",calcTempForw(),"\n</pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
          #write.table(HTML(paste0("<pre>\n\n<b><font size='4'>Primer information:</font size></b>\n\nThe selected forward primer sequence:\n\n",primerF(),"\n</pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)

        #calculate new begin and end position of the backward primer here as well
        #write.table(paste0("\n\nThe selected backward primer sequence:",HTML(paste0("\n<pre>",doColoring(input$SliderBack[1],IRanges::reverse(primerB()),vecSNPstartBack(),vecSNPwidth(),vecRepStartBack(),vecRepWidth(),60,NULL,NULL)), "\n\n",calcTempBack(),"\n</pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
          #write.table(HTML(paste0("<pre>\n\n<b><font size='4'>Primer information:</font size></b>\n\nThe selected reverse primer sequence:\n\n",primerB(),"\n</pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
          write.table(HTML(paste0("<pre>\n\n<b><font size='4'>Primer information:</font size></b></pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
          write.table(as.data.frame(primer.df()),file,quote = FALSE,row.names = T,col.names = T,append=T,sep="&emsp;",eol="<br>")

        #print information about the ampliconPrediction

        write.table(HTML(paste0("<pre>\n\n\n<b><font size='4'>AmpliconPrediction information:</font size></b>\n\n","range: ",input$SliderPlot[1]," - ",input$SliderPlot[2],"\n\n",
                                "$summary: TRUE/FALSE means the specific CG is/is not assayable by the given combination of cleavage reaction and DNA strand.\n\n",
                                "$counts: amount of CGs that are assayable by each assay\n\n",p(),"</pre>")),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
        #write table of designed primers in file
        write.table(paste0("\n\nPrimer3 designed the following primers for your input sequence:",designedPrimers()),file,quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)

      })







###################################################################
#                                                                 #
#                                                                 #
#                                                                 #
#                 BS-AMP-seq Functions                            #
#                                                                 #
#                                                                 #
#                                                                 #
###################################################################


  output$uiBisulfit <- renderUI({


# Choose Input option -----------------------------------------------------

    # choose which file format
    if(input$uiInputChoiceBisulfit == "fasta"){
      fileInput("fastaFile","Upload fasta file")
    }else if(input$uiInputChoiceBisulfit == "coord"){
      textAreaInput("genomicCoord", "Enter genomic coordinate",height="500px",placeholder="chr14:54,422,085-54,424,958/\nchr2:50,000,000-50,000,958/")
    }else{
      textAreaInput("cgID","Enter CpG ID",height="500px",placeholder=" cg12445832/\ncg16005420/")
    }
  })

    # choose which genome and then the respective file format (CpG-ID just for hg19)
  output$uiInputSelection <- renderUI({
    if(input$selectGenomeBisulfit=="hg19"){
      input_choices <- c("Genome coordinate" = "coord", "Fasta file" = "fasta","CpG ID" = "cpgID")
    }
    else{
      input_choices <- c("Fasta file" = "fasta", "Genome coordinate" = "coord")
    }
    selectInput("uiInputChoiceBisulfit","Select type of input",choices=input_choices,width="100%")
  })


  # enable download of the data
  tmp <- eventReactive(input$SubmitBTNRegionBisulfit,{

    output$downloadBut <- renderUI({
        downloadButton("downloadData_amp", "Download Table")
    })

    output$downloadData_amp <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(),"_BS_primer", ".csv")
      },
      content = function(file) {
        write.csv(tmp(), file, row.names = FALSE)
      }
    )


# Amplicons for CpG ID ----------------------------------------------------

    if(input$uiInputChoiceBisulfit == "cpgID"){

      # isolate single cpg ids
      req(input$cgID)
      comma <- stringr::str_split(input$cgID,"/")[[1]]
      print(comma)
      print(input$selectGenomeBisulfit)

      #advanced extenstion method
      if(!is.null(input$regionExtend)){
        if(input$regionExtend==T){
          ext <- extend_input()
          }else{
          ext <- 250
          }
      }
        else{
        ext <- 250
      }

      # generate amplicon from input cgID


      if(input$uiEPIC=="EPIC"){
       input_cgid <- cpg_coordinates(comma,epic=T)
      }else{
       input_cgid <- cpg_coordinates(comma)
       }
      my_amplicons <- generate_amplicons_coord(input_cgid,input$selectGenomeBisulfit,extend=ext,bs_convert = T)
      print(input$uiHandles)

      #advanced option implementation
      if(!is.null(input$incTarget)){
        if(input$incTarget==T)
          {cg_target <- paste0(targ.pos_input())}else{
            cg_target <- paste0(nchar(my_amplicons[[1]][1])/2)
        }
      }else{
        cg_target <- paste0(nchar(my_amplicons[[1]][1])/2)
      }


      # design primer for the amplicon
      print(extend_input)
      print(ext)
      print(input$incTarget)
      print(cg_target)
      primer_IDs <- primer_design(my_amplicons,names(my_amplicons),
                                  CGtarget = cg_target,
                                  primer_size = input$uiPrimS,
                                  primMaxSize = input$uiMaxPrimS,
                                  primMinSize = input$uiMinPrimS,
                                  primTm = input$uiPrimT,
                                  primMinTm = input$uiMinPrimT,
                                  primMaxTm = input$uiMaxPrimT,
                                  amp_min = input$uiAmpMin,
                                  amp_max = input$uiAmpMax,
                                  exclude_cpg = input$uiExcludeCG,
                                  primerNumber = input$nPrimer)

      # print table that contains all possible primers
      cutTable(input$uiHandles,adapter_fwd(),adapter_rev())
      primer.table <<- cutTable(input$uiHandles,adapter_fwd(),adapter_rev())


# Design primer based on genomic coordinates ------------------------------

    }else if(input$uiInputChoiceBisulfit == "coord"){

      req(input$genomicCoord)
      #single cpg coordinates
      coord.vec <- GenomicRanges::GRanges(stringr::str_split(input$genomicCoord,"/")[[1]])
      coord.vec_amp <- as.data.frame(stringr::str_split(input$genomicCoord,"/")[[1]])

      #advanced extenstion method
      if(!is.null(input$regionExtend)){
        if(input$regionExtend==T){
          ext <- extend_input()
        }else{
          ext <- 0
        }
      }
      else{
        ext <- 0
      }

      my_amplicons <- generate_amplicons_coord(as.data.frame(coord.vec_amp),input$selectGenomeBisulfit,extend=ext,bs_convert = T)

      #advanced option implementation
      if(!is.null(input$incTarget)){
        if(input$incTarget==T)
        {
          advanced <- T
          in_region <- GenomicRanges::GRanges(targ.coord_input())
          if(GenomicRanges::findOverlaps(coord.vec,in_region)@to>0){
            cg_target <- in_region@ranges@start-coord.vec@ranges@start
          }else{
            advanced <- F
            cg_target <- 0
          }

        }else{
          advanced <- F
          cg_target <- 0
        }
      }else{
        advanced <- F
        cg_target <- 0
      }

      # design primer for the amplicon
      primer_IDs <- primer_design(my_amplicons,names(my_amplicons),
                                  primer_size = input$uiPrimS,
                                  primMaxSize = input$uiMaxPrimS,
                                  primMinSize = input$uiMinPrimS,
                                  primTm = input$uiPrimT,
                                  primMinTm = input$uiMinPrimT,
                                  primMaxTm = input$uiMaxPrimT,
                                  amp_min = input$uiAmpMin,
                                  amp_max = input$uiAmpMax,
                                  exclude_cpg = input$uiExcludeCG,
                                  advSet = advanced,
                                  CGtarget = cg_target,
                                  primerNumber = input$nPrimer)

      # print table that contains all possible primers
      cutTable(input$uiHandles,#input$uiAdaptF,
               adapter_fwd(),adapter_rev())
      primer.table <<- cutTable(input$uiHandles,#input$uiAdaptF,
                                adapter_fwd(),adapter_rev())


# Input as FASTA ----------------------------------------------------------

    }else{
      req(input$fastaFile)
      #write sequences from fasta in list
      inFile <- input$fastaFile
      inFile <- inFile$datapath
      amplicons <- seqinr::read.fasta(inFile,as.string = T,set.attributes = F)
      names_amp <- names(amplicons)
      amplicons <- toupper(amplicons)
      names(amplicons) <- names_amp
      amp_minus <- secStrand(amplicons)
      amplicons <- c(amplicons,amp_minus)
      if (input$uiConvert == TRUE){
      amplicons <- bs_convert(amplicons)
      }
      #design primer for that list use names in fasta file for sequence ID
      attr(amplicons,"Amplicons") <- stringr::str_remove_all(names(amplicons),"_Minus")

      #advanced option implementation
      if(!is.null(input$incTarget)){
        if(input$incTarget==T)
        {
          advanced <- T
          cg_target <- paste0(targ.pos_input())}else{
          advanced <- F
          cg_target <- 0
        }
      }else{
        advanced <- F
        cg_target <- 0
      }

      primerFile <- primer_design(amplicons,names(amplicons),
                                  CGtarget = cg_target,
                                  primer_size = input$uiPrimS,
                                  primMaxSize = input$uiMaxPrimS,
                                  primMinSize = input$uiMinPrimS,
                                  primTm = input$uiPrimT,
                                  primMinTm = input$uiMinPrimT,
                                  primMaxTm = input$uiMaxPrimT,
                                  amp_min = input$uiAmpMin,
                                  amp_max = input$uiAmpMax,
                                  exclude_cpg = input$uiExcludeCG,
                                  advSet = advanced,
                                  primerNumber = input$nPrimer)

      # print table that contains all possible primers
      cutTable(input$uiHandles,adapter_fwd(),adapter_rev())
      primer.table <<- cutTable(input$uiHandles,adapter_fwd(),adapter_rev())
    }
  })



# Select the plot parameters ----------------------------------------------

  output$test <- DT::renderDataTable({

    tmp()

  },selection='single')

  observeEvent(input$test_rows_selected, {
    shinyBS::toggleModal(session, "modalPrimer", "open")
  })

  plot_chr <- reactive({
    if(as.character(primer.table[input$test_rows_selected,3]=="+")){
      as.character(GenomicRanges::GRanges(primer.table[input$test_rows_selected,15])@seqnames@values)}else{
        stringr::str_split_fixed(primer.table[input$test_rows_selected,15],":",2)[,1]
      }
    })

  plot_start <- reactive({
      GenomicRanges::GRanges((primer.table[input$test_rows_selected,15]))@ranges@start
      })

  plot_end <- reactive({
      GenomicRanges::GRanges(primer.table[input$test_rows_selected,15])@ranges@start+GenomicRanges::GRanges(primer.table[input$test_rows_selected,15])@ranges@width-1
    })

  fwd_start <- reactive({
    as.numeric(GenomicRanges::GRanges(primer.table[input$test_rows_selected,16])@ranges@start)-as.numeric(GenomicRanges::GRanges((primer.table[input$test_rows_selected,15]))@ranges@start)
  })

  fwd_end <- reactive({
    GenomicRanges::GRanges(primer.table[input$test_rows_selected,16])@ranges@width
  })

  rev_start <- reactive({
    as.numeric(GenomicRanges::GRanges(primer.table[input$test_rows_selected,17])@ranges@start)+as.numeric(GenomicRanges::GRanges(primer.table[input$test_rows_selected,17])@ranges@width)-as.numeric(GenomicRanges::GRanges((primer.table[input$test_rows_selected,15]))@ranges@start)-1
  })

  rev_end <- reactive({
    GenomicRanges::GRanges(primer.table[input$test_rows_selected,17])@ranges@width
  })

  strand <- reactive({
    #as.character(primer.table[input$test_rows_selected,3])
    as.character("+")
  })

  primerF <- reactive({
    if(!is.null(input$uiHandles)){
      stringr::str_remove_all(primer.table[input$test_rows_selected,4],input$uiAdaptF)
    }else{
      primer.table[input$test_rows_selected,4]
    }
  })

  primerR <- reactive({
    if(!is.null(input$uiHandles)){
      stringr::str_remove_all(primer.table[input$test_rows_selected,5],input$uiAdaptR)
    }else{
      primer.table[input$test_rows_selected,5]
    }
  })




# Create the plot ---------------------------------------------------------

  output$popup <- renderUI({
    shinyBS::bsModal("modalPrimer", "Chosen Primer", "", size = "large",
                       # div(
                       # conditionalPanel(id='loader',
                       #                condition="($('html').hasClass('shiny-busy'))",
                       #                img(src="loading_spinner.gif")),align="center"),
                     conditionalPanel(id='spin',
                                      condition="($('html').hasClass('shiny-busy'))",
                                      add_busy_spinner(spin = "fading-circle",position="full-page",
                                                       timeout = 1000)),
                     conditionalPanel(id='plot',
                                      condition="!($('html').hasClass('shiny-busy'))",
                       renderPlot(gviz_primer(input$selectGenomeBisulfit,
                                   plot_chr(),
                                   plot_start(),
                                   plot_end(),
                                   fwd_start(),
                                   fwd_end(),
                                   rev_start(),
                                   rev_end(),
                                   strand()
                                   )),
                       br(),
                       actionButton("BS_blast_primer","Blast Primer against BS genome"),
                       br(),
                       conditionalPanel(id='bs_primer_blast_cond',
                                        condition = "input.BS_blast_primer%2==1",
                                        tabsetPanel(
                                        tabPanel("Forward",
                                                 DT::renderDataTable(primer_bs_blast(primerF(),path_to_genome=paste0("/var/ressources/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa")))
                                                 ),
                                        tabPanel("Reverse",
                                                 DT::renderDataTable(primer_bs_blast(primerR(),path_to_genome=paste0("/var/ressources/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa")))
                                                 )
                                        ))
                       )
                       # h3("Primer blast against bisulfite converted reference genome"),
                       # h5("Forward"),
                       # DT::renderDataTable(primer_bs_blast(primerF(),path_to_genome=paste0("/var/ressources/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa"))),
                       # #MAC#DT::renderDataTable(primer_bs_blast(primerF(),path_to_bowtie="/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/bowtie",path_to_genome=paste0("/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/genomes/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa"))),
                       # h5("Reverse"),
                       # DT::renderDataTable(primer_bs_blast(primerR(),path_to_genome=paste0("/var/ressources/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa")))
                       # #MAC#DT::renderDataTable(primer_bs_blast(primerR(),path_to_bowtie="/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/bowtie",path_to_genome=paste0("/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/genomes/",input$selectGenomeBisulfit,"_bowtie_ct_ga_indexed.fa")))
                       )

  })


# Advanced Settings -------------------------------------------------------

  ## save the settings
  df.adv_settings <- reactive({
    data.frame("regEx"=if(is.null(input$regionExtend)){F}else{input$regionExtend},
               "inHandles"=if(is.null(input$uiHandles)){F}else{input$uiHandles},
               "inTarget"=if(is.null(input$incTarget)){F}else{input$incTarget},
               "inAdaptF"=if(is.null(input$uiHandles)){" "}else{input$uiAdaptF},
               "inAdaptR"=if(is.null(input$uiHandles)){" "}else{input$uiAdaptR},
               "inExtend"=if(is.null(input$regionExtend)){0}else{input$uiExtend},
               "inTarg.coord"=if(is.null(input$incTarget)){""}else{input$uiCG.coord},
               "inTarg.pos"=if(is.null(input$incTarget)){0}else{input$uiCG.pos})
  })

  ## reactive assignment to reset the input if checkboxes are unchecked
  adapter_fwd <- reactive({
    if(input$uiHandles==T){
      input$uiAdaptF
    }else{
      " "
    }
  })

  adapter_rev <- reactive({
    if(input$uiHandles==T){
      input$uiAdaptR
    }else{
     " "
    }
  })

  extend_input <- reactive({
    if(input$regionExtend==T){
      as.numeric(input$uiExtend)
    }else{
      0
    }
  })

  targ.coord_input <- reactive({
    if(input$incTarget==T){
      print(input$uiCG.coord)
    }else{
      ""
    }
  })

  targ.pos_input <- reactive({
    if(input$incTarget==T){
      as.numeric(input$uiCG.pos)
    }else{
      NULL
    }
  })


  observeEvent(input$advSet, {
    showModal(modalDialog(
      title = h1("Advanced Primer Design Settings"),
      h3("Choose which advanced settings should be applied"),
      br(),
      br(),
      checkboxInput("regionExtend","Extend the input genomic region",value=df.adv_settings()$regEx,width="100%"),
      checkboxInput("uiHandles","Add sequencing adapter",value=df.adv_settings()$inHandles,width="100%"),
      checkboxInput("incTarget","Specify a target region (not for batch input)",value=df.adv_settings()$inTarget,width="100%"),
      br(),
      br(),
      conditionalPanel(
        condition = "input.uiHandles == 1 || input.regionExtend ==1 || input.incTarget ==1",
        h3("Advanced Settings - Input Options")
      ),

      conditionalPanel(
        condition = "input.regionExtend == 1 && input.uiInputChoiceBisulfit != 'fasta'",
        numericInput("uiExtend","Extend the region in both directions by [input] base pairs",value=df.adv_settings()$inExtend,width="100%")
      ),
      conditionalPanel(
        condition = "input.uiHandles == 1",
        textInput("uiAdaptF","Adapter Sequence Fwd (5' - 3')",placeholder="CGTAGTCAGTCAGC",width="100%",value = df.adv_settings()$inAdaptF),
        textInput("uiAdaptR","Adapter Sequence Rev (5' - 3')",placeholder="CGTAGTCAGTCAGC",width="100%",value = df.adv_settings()$inAdaptR)
      ),
      conditionalPanel(
        condition = "input.incTarget == 1 && input.uiInputChoiceBisulfit == 'coord'",
        textInput("uiCG.coord","Target coordinates which should be included (separated by ´/´)",placeholder = "chr1:124324-124412",width="100%",value = df.adv_settings()$inTarg.coord)
      ),
      conditionalPanel(
        condition = "input.incTarget == 1 && input.uiInputChoiceBisulfit != 'coord'",
        numericInput("uiCG.pos","Target position of amplicon which should be included (separated by ´/´)",width="100%",value = df.adv_settings()$inTarg.pos)
      ),




      #Pop-Up Settings
      easyClose = TRUE,
      footer = modalButton("Confirm"),
      size="l"
    ))
  })


# Table Output of the used Settings ---------------------------------------

output$settings <- renderTable(
  data.frame(
    "Settings"=c("Selected Genome","Optimal Primer Size","Minimal Primer Size","Maximal Primer Size",
      "Optimal TM","Minimal TM","Maximal TM","Minimal Amplicon Size","Maximal Amplicon Size"),
    "Input"=c(input$selectGenomeBisulfit,input$uiPrimS,input$uiMaxPrimS,input$uiMinPrimS,input$uiPrimT,
              input$uiMinPrimT,input$uiMaxPrimT,input$uiAmpMin,input$uiAmpMax)
  )
)



# Switch to primer design tab ---------------------------------------------


  observeEvent(input$SubmitBTNRegionBisulfit, {
    updateTabsetPanel(session, "BSampTabset",
                      selected = "designedPrimer")
  })

  # observeEvent(input$SubmitBTNRegion, {
  #   updateTabsetPanel(session, "MassArrayTabset",
  #                     selected = "inputRegion_tabset")
  # })

  observeEvent(input$submit_blast, {
    updateTabsetPanel(session, "bisblast",
                      selected = "blast_query")
  })

  # observeEvent(input$submit_analysis, {
  #   updateTabsetPanel(session, "analysispanel",
  #                     selected = "analysis_selected")
  # })



# Load the sample data for AmpBS-Seq ---------------------------------------------------------------

observeEvent(input$AMP_sampleData,{
  updateSelectInput(session,"selectGenomeBisulfit",selected = "hg19")
  updateSelectInput(session,"uiInputChoiceBisulfit",selected = "coord")
  updateTextInput(session,"genomicCoord",value="chr19:43203328-43203689")
})


# BisBlast ----------------------------------------------------------------

  blast_coord <- eventReactive(input$submit_blast,{
    input$blast_in
  })

  output$blast_table <- DT::renderDataTable(
    primer_bs_blast(blast_coord(),path_to_genome=paste0("/var/ressources/",input$selectGenome_blast,"_bowtie_ct_ga_indexed.fa")))
    #MAC#primer_bs_blast(blast_coord(),path_to_bowtie="/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/bowtie",path_to_genome=paste0("/Users/maximilianschoenung/Documents/AG_Lipka/Computational/software/bowtie-1.2.2-macos-x86_64/genomes/",input$selectGenome_blast,"_bowtie_ct_ga_indexed.fa")))

  # observeEvent(input$submit_blast, {
  #   primer_bs_blast(input$blast_in,path_to_genome=paste0("/var/ressources/",input$selectGenome_blast,"_bowtie_ct_ga_indexed.fa"))
  # })

# Interactive data from help ----------------------------------------------

  observeEvent(input$interactive_AMP,{
    updateTabsetPanel(session, "page",
                      selected = "AmpliconBisulfiteSequencing")
    #click("AMP_sampleData")
    #updateSelectInput(session,"selectGenomeBisulfit",selected = "hg19")
    #updateSelectInput(session,"uiInputChoiceBisulfit",selected = "coord")
    delay(3000,updateTextInput(session,"genomicCoord",value="chr19:43203328-43203689"))
    delay(3000,click("SubmitBTNRegionBisulfit"))
    updateTabsetPanel(session, "BSampTabset",
                      selected = "designedPrimer")

  })

  observeEvent(input$interactive_MA,{
    updateTabsetPanel(session, "page",
                      selected = "MassArray")
    click("MA_sampleData")
    delay(3000,click("SubmitBTNRegion"))
    updateTabsetPanel(session, "MassArrayTabset",
                      selected = "inputRegion_tabset")

  })




# Analysis Pipeline -------------------------------------------------------

  # file_upload <- eventReactive(input$submit_analysis,{
  #   req(input$input_analysis_files)
  #   file_analysis <- input$input_analysis_files
  #   file_analysis <- file_analysis$datapath
  #   print(file_analysis)
  # })
  #
  # meta_upload <- eventReactive(input$submit_analysis,{
  #   req(input$input_analysis_sample)
  #   file_sample <- input$input_analysis_sample
  #   file_sample <- file_sample$datapath
  #   file <- read.delim(file_sample,header=T)
  #   rownames(file) <- file[,1]
  #   file
  # })

  methSet <- eventReactive(input$submit_analysis,{
    file_analysis <- input$input_analysis_files
    print(file_analysis$name)
    file_analysis2 <- file_analysis$datapath
    print(file_analysis2)
    file_sample <- input$input_analysis_sample
    file_sample <- file_sample$datapath
    load_amp(file_analysis2,file_sample,file_analysis$name)
  })

  meta <- eventReactive(input$submit_analysis,{
    file_sample <- input$input_analysis_sample
    file_sample <- file_sample$datapath
    read.delim(file_sample,row.names = 1)
  })

  bed <- eventReactive(input$submit_analysis,{
    file_analysis <- input$input_analysis_regions
    file_analysis <- file_analysis$datapath
    cg_in <- read.delim(file_analysis,header = F)
    GenomicRanges::GRanges(seqnames = cg_in$V1, ranges = IRanges::IRanges(cg_in$V2,width=cg_in$V3-cg_in$V2+1) ,cg_is=cg_in$V4)
    })
  
  format_bed <- reactive({
    req(bed())
    df <- as.data.frame(bed())
    paste0(df$seqnames,":",df$start,"-",df$end)
  })

  output$dt_overview <- renderTable({
    data.frame("Paramter"=c("Number of Samples","Meta-Data","Sites (Total)","Sites (Coverage-Filtered)","Regions (bed-Input)"),
    "Value"=c(length(colnames(methSet())),
    stringr::str_flatten(colnames(methSet()@colData),", "),
    nrow(methSet()@elementMetadata),
    nrow(BiSeq::filterByCov(methSet(),as.numeric(100),global=FALSE)@elementMetadata),
    nrow(bed()@elementMetadata)))
  },align="c")

  cg_table <- reactive({
    cg_tab <- cbind("Covered CpG"=BiSeq::covStatistics(methSet())$Covered_CpG_sites,
                    "Median Coverage"=BiSeq::covStatistics(methSet())$Median_coverage)
    rownames(cg_tab) <- meta()$UPN
    cg_tab
  })


  output$cov_sites <- renderTable({
    cg_table()
  },rownames = T)



  # output$qc_totalReads <- renderPlot({
  #   total.reads <- BiSeq::totalReads(methSet())
  #   total.reads <- data.frame("Sample"=colnames(total.reads),"Sites"=colSums(total.reads))
  #   ggplot(total.reads,aes(Sample,log10(Sites)))+
  #     geom_bar(stat="identity")+
  #     theme_minimal()+
  #     ylab("log10(Total Reads per Sample)")+
  #     theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # })


  ### Generate qc Plots

  amp_qc_raw <- reactive({
    qc_plots(methSet())
  })

  output$qc_totalReads <- renderPlot({
    plot(amp_qc_raw()[[1]])
  })

  output$qc_covSites <- renderPlot({
    plot(amp_qc_raw()[[2]])
  })

  output$qc_cor <- renderPlot({
    plot(amp_qc_raw()[[3]])
  })

  output$qc_covCG <- renderPlot({
    plot(amp_qc_raw()[[4]])
  })


  ### Coverage FIltering and QC Plots

  methSet_filt <- reactive({
    Mset <- BiSeq::filterByCov(methSet(),as.numeric(input$cutoff),global=FALSE)
    IRanges::subsetByOverlaps(Mset,bed())
  })

  methSet_filtR <- reactive({
    IRanges::subsetByOverlaps(methSet_filt(),bed())
  })

  amp_qc_filt <- reactive({
    qc_plots(methSet_filt())
  })

  output$cf_covSites <- renderPlot({
    plot(amp_qc_filt()[[2]])
  })

  output$cf_cor <- renderPlot({
    plot(amp_qc_filt()[[3]])
  })

  output$cf_covCG <- renderPlot({
    plot(amp_qc_filt()[[4]])
  })

  ## pca plot
  ### renderUI seems to make it slow...
  output$select_pca = renderUI({
    req(methSet_filtR())
    selectInput('pca_color', 'Coloring of PCA', colnames(methSet_filtR()@colData))
  })

  output$pca_ampBS <- renderPlot({
    clustering(methSet_filtR(),input$pca_color)
  })

  output$select_heatmap <-  renderUI({
    req(methSet_filtR())
    selectInput('select_heatmap', 'Heatmap Annotation', colnames(methSet_filtR()@colData),selected = "UPN")
  })

  output$heatmap_rows <-  renderUI({
  checkboxInput("clusterRow","Cluster Rows")})

  output$heatmap_cols <-  renderUI({
    checkboxInput("clusterCol","Cluster Columns")
    })

  output$heatmap_colNames <-  renderUI({
    checkboxInput("colnames","Show column names")
  })

  output$heatmap_rowNames <-  renderUI({
    checkboxInput("rownames","Show row names")
  })

  output$heatmap_coloring <-  renderUI({
    selectInput("coloring","Choose Color Palette",
                c("blue-yellow-red"="byr",
                  "green-white-red"="gwr",
                  "blue-white-red"="bwr"))
  })

  output$heatmap_ampBS <- renderPlot({
    if(input$coloring=="byr"){
      col_selected <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    }else if(input$coloring=="gwr"){
      col_selected <- colorRampPalette(c("forestgreen","white","red"))(100)
    }else if(input$coloring=="bwr"){
      col_selected <- colorRampPalette(c("navy","white","red"))(100)
    }
    req(methSet_filtR())
    amp_heatmap(methSet_filtR(),as.character(input$select_heatmap),cluster_cols=input$clusterCol,cluster_rows=input$clusterRow,
                show_colnames=input$colnames,show_rownames=input$rownames,border_color=NA,
                color=col_selected)
  })
  
  output$select_region <-  renderUI({
    req(methSet_filtR())
    selectInput('select_regionHeatmap', 'Select Amplicon Region', format_bed())
  })
  
  output$co_region_ui <-  renderUI({
    checkboxInput("co_region","Apply Coverage Cut-Off")
  })
  
  output$heatmap_region <- renderPlot({
    req(methSet())
    print(input$select_regionHeatmap)
    region_heatmap(methSet(),input$select_regionHeatmap,co=input$co_region,co_num=as.numeric(input$cutoff))
  })


###### Download of AmpBS-Seq Data


    output$downloadBut_ampBS <- renderUI({
      downloadButton("amp_report", "Download Report")
    })

    output$amp_report <- downloadHandler(
      filename = "report_ampBS_analysis.html",
      content = function(file_amp){

        tempReport <- file.path(tempdir(), "report_ampBS_analysis.Rmd")
        file.copy("report_ampBS_analysis.Rmd", tempReport, overwrite = TRUE)

        #Set up parameters to pass to Rmd document
        params <- list(set = methSet(),
                       set_filt = methSet_filt(),
                       set_filtR = methSet_filtR(),
                       cov_co = input$cutoff,
                       bed_regions = bed(),
                       meta_df = meta())

        rmarkdown::render(tempReport, output_file = file_amp,
                          params = params,
                          envir = new.env(parent = globalenv()))
      }
    )


  # enable download of count and coverage matrix and report
  # make sample data available
  # write help section

# Load sample data for the pipeline ---------------------------------------

    output$Pipeline_sampleData <- downloadHandler(
      filename = "sampleData_pipeline.zip",
      content = function(file_pipeline){
        file.copy("www/sampleData_pipeline.zip", file_pipeline)
      }
    )




# Warnings AmpBS Panel ----------------------------------------------------

  ## check if the UPN column and the rownames exist in the meta sheet
  # check if the meta names match with the sample files
    
    ## first generate a reactive value to handle error messages when an action button is pressed
    output$clicked_ampBS <- renderText("no")
    
    outputOptions(output, "clicked_ampBS", suspendWhenHidden=FALSE)
    
    observeEvent(input$submit_analysis,{
    file_analysis <- input$input_analysis_files
    output$clicked_ampBS <- renderText("yes")
    if(is.null(input$input_analysis_files)){
      showNotification(paste("Please Choose Input Files"), duration = 0,type="error")
      output$clicked_ampBS <- renderText("no")}
    else if(is.null(input$input_analysis_sample)){
      showNotification(paste("Please Choose a Sample Sheet"), duration = 0,type="error")
      output$clicked_ampBS <- renderText("no")
      }
    else if(is.null(input$input_analysis_regions)){
      showNotification(paste("No regions have been selected. Please insert your query regions in bed-file format."), duration = 0,type="warning")
      output$clicked_ampBS <- renderText("no")
    }else if(length(file_analysis$name)!=nrow(meta())){
      showNotification(paste("Number of samples in sample sheet does not match uploaded files."), duration = 0,type="error")
      output$clicked_ampBS <- renderText("no")
    }else if(!c("UPN")%in%colnames(meta())){
      showNotification(paste("The sample sheet must contain a unique patient number (UPN) column."), duration = 0,type="error")
      output$clicked_ampBS <- renderText("no")
    }else if(sum(as.numeric(file_analysis$name%in%rownames(meta())))!=length(file_analysis$name)){
      showNotification(paste("First column in the sample sheet do not exactly match the filenames of the Bismark coverage files."), duration = 0,type="error")
      output$clicked_ampBS <- renderText("no")
    }else{
      updateTabsetPanel(session, "analysispanel",
                        selected = "analysis_selected")
    }
    # ## check the quality of the uploaded files
    
    })
    

# Warning MassArray -------------------------------------------------------
    output$clicked_ma <- renderText("no")
    outputOptions(output, "clicked_ma", suspendWhenHidden=FALSE)
    observeEvent(input$SubmitBTNRegion,{
    region <- inReg()
    chr.ind <- stringr::str_split_fixed(region,":",2)[,1]
    num.ind <- stringr::str_split_fixed(region,":",2)[,2]
    num.start <- as.numeric(stringr::str_split_fixed(num.ind,"-",2)[,1])
    num.end <- as.numeric(stringr::str_split_fixed(num.ind,"-",2)[,2])
    
    output$clicked_ma <- renderText("yes")
    if(!chr.ind%in%paste0("chr",c(1:22,"X","Y"))){
      showNotification(paste("The chromosome index does not exist. Please insert in this format 'chr19:43204328-43203689'."), duration = 0,type="error")
      output$clicked_ma <- renderText("no")
    }else if(num.start>num.end){
      showNotification(paste("Chromosome start cannot be greater than end. Please insert in this format 'chr19:43204328-43203689'."), duration = 0,type="error")
      output$clicked_ma <- renderText("no")
    }else{
        updateTabsetPanel(session, "MassArrayTabset",
                          selected = "inputRegion_tabset")
    }
      })
    
    

# Warning AmpBS-Seq -------------------------------------------------------
    output$clicked_bs <- renderText("no")
    outputOptions(output, "clicked_bs", suspendWhenHidden=FALSE)
    observeEvent(input$SubmitBTNRegionBisulfit,{
      if(input$uiInputChoiceBisulfit == "fasta"){
        ## check for valid sequence and genome coordinate style header
      }else if(input$uiInputChoiceBisulfit == "coord"){
        ## check for proper coordinate
      }else if(input$uiInputChoiceBisulfit == "cpgID"){  
        ## check if this is actually on the array IF nrow dt == 0 (sonst rechnet das zu lange)
      }else{
        ## check primer design parameter IF nrow dt == 0 
        ## also give a warning if this happens that the region might either be CG rich or the parameters do not fit
      }
    })
    

}


