
###################################################################
#                                                                 #
#                                                                 #
#                                                                 #
#                 MassArray Functions                            #
#                                                                 #
#                                                                 #
#                                                                 #
###################################################################

dbQuery <- function(chr,begin,end,genome,S=T,R=F){
  start_time_q <- Sys.time()
  my_db = DBI::dbConnect(RPostgreSQL::PostgreSQL(), user="postgres", password="postgres",
                         host="c010-shiny2", port=5432, dbname="testdb")
  on.exit(DBI::dbDisconnect(my_db))
  if(S==T){
   query <-  as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM snp WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))"))))
  }else if(R==T){
    query <- as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM cgi WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))"))))
  }
  end_time_q <- Sys.time()
  message(paste0("Time for the Postgres DB Query:",c(end_time_q - start_time_q)))
  if(nrow(query)!=0){
    track <- Gviz::AnnotationTrack(query)
  }else{
    track <- Gviz::AnnotationTrack(NULL)  
    }
  return(track)
}




doColoring <- function(begin,seq,SNPposstart,SNPposwidth,rbegin,rwidth,linebreak,startPrimer, endPrimer){
  seq <- unlist(strsplit(seq,""))
  str <-""

  for(i in 1:length(seq)){
    if (!is.null(startPrimer) && i == startPrimer){
      str <- paste0(str,"<mark>")
    }
    printedSNP <-FALSE
    for(z in 1:length(SNPposstart)){
      #SNP starting positions are returned by UCSC one position before the actual start of SNP
      #width is always given with 2 bases too much

      if(length(SNPposstart)>0 && (begin+i-1) >= SNPposstart[z] && (begin+i-1) < SNPposstart[z]+SNPposwidth[z] ){
        str <-paste0(str,"<font color='violet'>",seq[i],"</font>")
        printedSNP =TRUE
        break
      }
    }
    if(!printedSNP){
      one <- FALSE
      two<- FALSE
      if(i<length(seq)){
        one <- (seq[i] =="C" && seq[i+1] == "G")
      }
      if(i > 1){
        two <- seq[i] =="G" && seq[i-1] =="C"
      }
      if(one||two){
        str <-paste0(str,"<font color='red'>",seq[i],"</font>")
      }else{
        printed =FALSE
        for(x in 1:length(rbegin)){
          #repeat starting position is returned by UCSC browser one base before the actual begin of the repeat
          #width of repeat is returned by UCSC browser by one base too much
          if(length(rbegin)>0 && ((begin+i-1) >= rbegin[x]) && ((begin+i-1) < rbegin[x]+rwidth[x])){#if base is part of a repeat color it green
            str <- paste0(str,"<font color='green'>",seq[i],"</font>")
            printed = TRUE
            break
          }
        }
        if(!printed){#else just add letter
          str <- paste0(str,seq[i])
        }
      }
    }
    if(!is.null(endPrimer) && i == endPrimer){
      str <- paste0(str,"</mark>")
    }
    if((i%%linebreak) == 0){
      str <- paste0(str,"<br/>")
    }
  }
  return(str)
}


doColoring_rev <- function(begin,seq,SNPposstart,SNPposwidth,rbegin,rwidth,linebreak,startPrimer, endPrimer){
  seq <- unlist(strsplit(seq,""))
  str <-""

  for(i in 1:length(seq)){
    if (!is.null(startPrimer) && i == startPrimer){
      str <- paste0(str,"<mark>")
    }
    printedSNP <-FALSE
    for(z in 1:length(SNPposstart)){
      #SNP starting positions are returned by UCSC one position before the actual start of SNP
      #width is always given with 2 bases too much

      if(length(SNPposstart)>0 && (begin+i-1) >= SNPposstart[z] && (begin+i-1) < SNPposstart[z]+SNPposwidth[z] ){
        str <-paste0(str,"<font color='violet'>",seq[i],"</font>")
        printedSNP =TRUE
        break
      }
    }
    if(!printedSNP){
      one <- FALSE
      two<- FALSE
      if(i<length(seq)){
        one <- (seq[i] =="G" && seq[i+1] == "C")
      }
      if(i > 1){
        two <- seq[i] =="C" && seq[i-1] =="G"
      }
      if(one||two){
        str <-paste0(str,"<font color='red'>",seq[i],"</font>")
      }else{
        printed =FALSE
        for(x in 1:length(rbegin)){
          #repeat starting position is returned by UCSC browser one base before the actual begin of the repeat
          #width of repeat is returned by UCSC browser by one base too much
          if(length(rbegin)>0 && ((begin+i-1) >= rbegin[x]) && ((begin+i-1) < rbegin[x]+rwidth[x])){#if base is part of a repeat color it green
            str <- paste0(str,"<font color='green'>",seq[i],"</font>")
            printed = TRUE
            break
          }
        }
        if(!printed){#else just add letter
          str <- paste0(str,seq[i])
        }
      }
    }
    if(!is.null(endPrimer) && i == endPrimer){
      str <- paste0(str,"</mark>")
    }
    if((i%%linebreak) == 0){
      str <- paste0(str,"<br/>")
    }
  }
  return(str)
}



calcMeltTemp <- function(seq){

  amountA <-stringr::str_count(seq,"A")
  amountT <-stringr::str_count(seq,"T")
  amountG <-stringr::str_count(seq,"G")
  amountC <-stringr::str_count(seq,"C")
  amountN <-stringr::str_count(seq,"N")
  amountY <-stringr::str_count(seq,"Y")
  assertthat::assert_that((amountN+amountA+amountC+amountG+amountT+amountY) == nchar(seq))
  return(((amountA+amountT)*2+(amountC+amountG+amountY)*4)-5)
}

# discoverSequence2 <-function(chromosome, begin, end, genome){
#   #get sequence with the help of chromosome and range
# 
#   if(genome == "hg38"){
#     g <-BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#   }else{
#     if(genome == "mm10"){
#       g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#     }else{
#       (genome == "hg19")
#       g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#     }
#   }
#   w <- end-begin+1
#   gr1 <- GenomicRanges::GRanges(seqnames = c(paste("chr",chromosome,sep = "")),ranges =IRanges::IRanges(start=begin,width=w))
#   names(gr1) <- head("A",1)
#   return(rtracklayer::getSeq(g,gr1))
# }

discoverSequence <-function(chromosome, begin, end, genome){
#discoverSequence("chr1",1565681,1565881,"hg19")
  genome.bed <<- tempfile("in_bedtools")
  genome_out.bed <<- tempfile("out_bedtools")
  bed <- matrix(c(chromosome,begin,end),byrow = F,ncol=3)
  write.table(bed,file=genome.bed,sep="\t",quote=F,row.names = FALSE,col.names = FALSE)
  system(paste0("bedtools getfasta -fi /var/ressources/fasta/",genome,".fa -bed ",genome.bed," -fo ",genome_out.bed))
  #system(paste0("bedtools getfasta -fi /Users/maximilianschoenung/Documents/AG_Lipka/Computational/annotation/hg19/",genome,".fa -bed ",genome.bed," -fo ",genome_out.bed))
  amp <- read.delim(genome_out.bed,strip.white = T,stringsAsFactors = F,header = F)[2,1]
  
  
  return(toupper(amp))
}


bisulfitConvert <- function(seq){
  #seq <- toupper(seq)
  seq <- gsub("CG","YG",seq)
  seq <- gsub("C","T",seq)
  seq <- gsub("YG","CG",seq)
  return(seq)
}


bisulfitConvert_rev <- function(seq){
  #seq <- toupper(seq)
  seq <- gsub("GC","GY",seq)
  seq <- gsub("C","T",seq)
  seq <- gsub("GY","GC",seq)
  return(seq)
}




checkLetters <- function(sequence){
  #turn sequence from a DNAStringObject into a vector of characters
  sequence1 <- unlist(strsplit(sequence,""))
  validLetters <- c("A", "C", "G", "T", "N","a","c","g","t","n")
  return(all(sequence1 %in% validLetters))
}



findCG <- function(sequence){
  return(stringr::str_locate_all(sequence,"CG"))
}




doAmpliconPrediction <- function(seq){
  return(MassArray::ampliconPrediction(seq))
}


makeReverseComplement <- function(seq){
  return(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
}

makeComplement <- function(seq){
  x <- Biostrings::DNAString(seq)
  x <- Biostrings::complement(x)
  return(as.character(x))
}


# Mass-Array Fragments ----------------------------------------------------

mass_output <- function(input_ranges,genome){
  string_sub <- stringr::str_split_fixed(input_ranges,":",2)
  string_sub1 <- string_sub[,2]
  string_sub1 <- stringr::str_split_fixed(string_sub1,"-",2)
  string_sub1 <- apply(string_sub1,2,as.numeric)
  if(string_sub1[1]>string_sub1[2]){
    select_ranges <- paste0(string_sub[,1],":",string_sub1[2],"-",string_sub1[1])
  }else{
    select_ranges <- input_ranges
  }
  chr <- as.numeric(stringr::str_remove_all(GenomicRanges::GRanges(select_ranges)@seqnames,"chr"))
  start <- GenomicRanges::GRanges(select_ranges)@ranges@start
  end <- GenomicRanges::GRanges(select_ranges)@ranges@start+GenomicRanges::GRanges(select_ranges)@ranges@width-1
  # extract the sequence
  sequence <- as.character(discoverSequence(chr,start,end,genome))
  if(string_sub1[1]>string_sub1[2]){
    sequence <- as.character(makeReverseComplement(Biostrings::DNAString(sequence)))
  }
  # count the CpGs in sequence
  cgs <- MassArray::countCGs(sequence)
  # Analyze the putative fragments
  frags <- MassArray::inSilicoFragmentation(sequence)
  # Amplicon prediction
  preds <- MassArray::ampliconPrediction(sequence,plot = F)
  # Fragment output
  y <- frags
  for(i in 1:length(y)){
    if(i==1){
      x <- y[[i]]
      df <- data.frame("ID"=x@ID,
                       "Sequence"=x@sequence,
                       "Position"=x@position,
                       "Length"=x@length,
                       "CpGs"=x@CpGs,
                       "MW"=x@MW,
                       "Collision"=x@collisions,
                       "Cleavage"=x@type,
                       "Strand"=x@direction)
    }else{
      x <- y[[i]]
      df <- rbind(df,data.frame("ID"=x@ID,
                                "Sequence"=x@sequence,
                                "Position"=x@position,
                                "Length"=x@length,
                                "CpGs"=x@CpGs,
                                "MW"=x@MW,
                                "Collision"=x@collisions,
                                "Cleavage"=x@type,
                                "Strand"=x@direction))}}
  print(list("Sequence"=sequence,
             "CpGs"=cgs,
             "AmpliconPreds"= preds,
             "Fragments"=df[which(df$MW>1499),]))

}


###################################################################
#                                                                 #
#                                                                 #
#                                                                 #
#                 BS-AMP-seq Functions                            #
#                                                                 #
#                                                                 #
#                                                                 #
###################################################################


bs_convert <- function(x){
  amplicons <- lapply(x,function(x)toupper(x))
  amplicons <- lapply(amplicons,function(x)stringr::str_replace_all(x,"C","t"))
  return(amplicons)
}

primer_design <- function(fastas,
                          namesx=names(fastas),
                          CGtarget=paste0(nchar(fastas[[1]][1])/2),
                          primer_size=20,primMaxSize=22,primMinSize=15,
                          primTm=55,primMinTm=50,primMaxTm=62,
                          amp_min=150,amp_max=300,
                          exclude_cpg=T,primerNumber=10,
                          advSet=T){


  #assess the CGs for later exclusion
  cg_list <- lapply(fastas,function(x)gregexpr("tG",x))


  #create input file
  input.df <- matrix()
  input.df <- rbind(input.df,
                    paste0("PRIMER_OPT_SIZE=",primer_size),
                    paste0("PRIMER_MIN_SIZE=",primMinSize),
                    paste0("PRIMER_MAX_SIZE=",primMaxSize),
                    paste0("PRIMER_NUM_RETURN=",primerNumber),
                    paste0("PRIMER_OPT_TM=",primTm),
                    paste0("PRIMER_MIN_TM=",primMinTm),
                    paste0("PRIMER_MAX_TM=",primMaxTm),
                    paste0("PRIMER_PRODUCT_SIZE_RANGE=",amp_min,"-",amp_max),
                    "PRIMER_FILE_FLAG=0")
  for (i in 1:length(fastas)){
    input.df <- rbind(input.df,
                      paste0("PRIMER_SEQUENCE_ID=",namesx[i]),
                      paste0("COMMENT=",attributes(fastas)$Amplicons[[i]]),
                      paste0("SEQUENCE=",fastas[[i]][1]),
                      if(advSet==T){
                        if(stringr::str_detect(namesx[i],"Minus")){
                          paste0("TARGET=",nchar(fastas[[i]][1])-as.numeric(CGtarget),",2")
                        }else{
                          paste0("TARGET=",CGtarget,",2")}
                        })
    if (exclude_cpg == TRUE){
      for(j in 1:length(cg_list[[i]][[1]])){
        input.df <- rbind(input.df,
                          paste0("EXCLUDED_REGION=",cg_list[[i]][[1]][j],",2"))
      }
    }
    input.df <- rbind(input.df,
                      "=")
  }
  in.temp <<- tempfile("in")
  write.table(input.df[2:nrow(input.df),],file=in.temp,sep="\t",quote=F,row.names = FALSE,col.names = FALSE)

  #create an output directory
  output <<- tempfile("out")

  #run Primer3 bash code
  if(Sys.info()["sysname"]=="Linux"){
    primer3 <- paste0(getwd(),"/software/primer3_linux/primer3_core")
  }else{
  primer3 <- paste0(getwd(),"/software/primer3_mac/primer3_core")}
  system(
    paste(primer3,"< ",in.temp," > ",output))

  #read output file
  primers <- read.table(output)
  return(primers)

}



cpg_coordinates <- function(x,epic=F){
  if(epic==T){
  return(data.frame(ampliconDesignData::hm850k[x]))
  }else{
  return(data.frame(ampliconDesignData::hm450k[x]))}
}


generate_amplicons_coord <- function(x,g,extend=0,bs_convert=TRUE){
 #test: generate_amplicons_coord(as.data.frame("chr1:1565681-1565281"))
  
  # #select the genome of choice
  # if(g=="hg19"){
  #   genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  # }else if(g=="hg38"){
  #   genome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  # }else{
  #   genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  # }
  
  ############################################
  # import the genome data of choice using bedtools


  # generate the amplicons and save in a list
  # amp<-list()
  # amp_att <- list()
  # for(i in 1:nrow(x)){
  # 
  #   chr <- as.character(x$seqnames)[i]
  #   start <- as.numeric(x$start[i])
  #   end <- as.numeric(x$end[i])
  #   amp[[i]] <- genome[[paste0(chr)]][(start-extend):(end+extend)]
  #   names(amp)[i] <- paste0(chr,":",start,"-",end)
  #   amp_att[[i]] <- paste0(chr,":",start-extend,"-",end+extend)
  # 
  # }

  ############################################
  # import the genome data of choice using bedtools
  in.bed <<- tempfile("in_bedtools")
  out.bed <<- tempfile("out_bedtools")
  splitted <- stringr::str_split_fixed(x[,1],":",2)
  splitted.pos <- stringr::str_split_fixed(splitted[,2],"-",2)
  bed <- data.frame(splitted[,1],(as.numeric(splitted.pos[,1])-extend),(as.numeric(splitted.pos[,2])+extend))
  write.table(bed,file=in.bed,sep="\t",quote=F,row.names = FALSE,col.names = FALSE)
  system(paste0("bedtools getfasta -fi /var/ressources/fasta/",g,".fa -bed ",in.bed," -fo ",out.bed))
  amp.tab <- read.delim(out.bed,strip.white = T,stringsAsFactors = F,header = F)
  amp <- as.list(amp.tab[(seq(1,nrow(amp.tab),2)+1),1])
  names(amp) <- stringr::str_remove_all(amp.tab[(seq(1,nrow(amp.tab),2)),1],">")
  
  
  #generate the second strand
  amp <- c(amp,secStrand(amp))
  amp <- lapply(amp,toupper)


  #BS convert
  if(bs_convert==TRUE){
    amp<- bs_convert(amp)
  }

  #add amplicon genomic region as attributes
  attr(amp,"Amplicons") <- rep(stringr::str_remove_all(amp.tab[(seq(1,nrow(amp.tab),2)),1],">"),2)

  #return the amplicons
  return(amp)
}



cutTable <- function(d,handles_f,handles_r){

  output.primer <- read.table(output)
  ret <- matrix(nrow=10,ncol =17)

  #divide file at each new amplicon
  output_split <- stringr::str_split(toString(output.primer$V1),c("PRIMER_SEQUENCE_ID[_]?[0-9]?[0-9]?[0-9]?[0-9]?"))
  text <- output_split[[1]]
  #for each amplicon produce a matrix that contains the primer info and add the matrices together
  for(i in 2:length(text)){
    tmp <- ret
    ret <- produceOutput(strsplit(text[i],"[,=]")[[1]],i-1)
    if(!is.na(tmp[1][1])){
      ret <- rbind(tmp,ret)
    }
  }
  if(!is.null(d)){
    ret[,4] <- paste0(tolower(handles_f),ret[,4])
    ret[,5] <- paste0(tolower(handles_r),ret[,5])
  }

  return(ret)
}


produceOutput <- function(output_split,numberOfAmplicon){

  #save primer id from output
  id <- toString(output_split[2])
  ucsc <- toString(output_split[4])
  ucsc_ranges <- GenomicRanges::GRanges(ucsc)

  #cut string at every new primer location
  
  cut <- grep("^ PRIMER_LEFT*.*.SEQUENCE$",output_split)
  cut <- append(cut,length(output_split))
  if(length(cut) > 0){

    x <- 1
    mat <- matrix(dimnames = list(c(),c("Amplicon","Primer ID","Strand","Primer sequences left","Primer sequences right","Primer left begin","Primer left length","Primer right end","Primer right length","GC percent left","GC percent right","Melting temp left","Melting temp right","Amplicon size","Region","UCSC Primer Left","UCSC Primer Right","Amplicon","#CpGs")),nrow = length(cut)-1,ncol = 19)
    count <- 1
    for(i in cut){
      if( x != 1){
        splitted <- output_split[x:i-1]
        #save which amplicon is currently worked on (1st,2nd,...)
        mat[,1] <- numberOfAmplicon
        #save the ID of the amplicon also
        mat[,2] <- id
        #strand
        mat[,3] <- if(stringr::str_detect(id,"Minus")==TRUE){"-"}else{"+"}
        #save the primer sequences
        mat[count,4] <- toString(splitted[grep("PRIMER_LEFT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?SEQUENCE",splitted,ignore.case = T)+1])
        mat[count,5] <- toString(splitted[grep("PRIMER_RIGHT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?SEQUENCE",splitted,ignore.case = T)+1])
        #save the start and end positions of the primers
        mat[count,6] <- toString(splitted[grep("PRIMER_LEFT[_]?[0-9]?[0-9]?[0-9]?[0-9]?",splitted,ignore.case = T)+1][2])
        mat[count,7] <- toString(splitted[grep("PRIMER_LEFT[_]?[0-9]?[0-9]?[0-9]?[0-9]?",splitted,ignore.case = T)+2][2])

        mat[count,8] <- toString(splitted[grep("PRIMER_RIGHT[_]?[0-9]?[0-9]?[0-9]?[0-9]?",splitted,ignore.case = T)+1][2])
        mat[count,9] <- toString(splitted[grep("PRIMER_RIGHT[_]?[0-9]?[0-9]?[0-9]?[0-9]?",splitted,ignore.case = T)+2][2])
        if(length(grep("Minus",id))!=0){
          rev.start <- as.numeric(mat[count,6])
          rev.length <- as.numeric(mat[count,7])
          fwd.start <- as.numeric(mat[count,8])
          fwd.length <- as.numeric(mat[count,9])
          mat[count,6] <- as.numeric(ucsc_ranges@ranges@width)-rev.start-rev.length
          mat[count,8] <- as.numeric(ucsc_ranges@ranges@width)-fwd.start+fwd.length-2
        }
        #save the GC percent of the primers
        mat[count,10] <- toString(splitted[grep("PRIMER_LEFT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?GC",splitted,ignore.case = T)+1])
        mat[count,11] <- toString(splitted[grep("PRIMER_RIGHT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?GC",splitted,ignore.case = T)+1])
        #save the predicted melting temperature
        mat[count,12] <- toString(splitted[grep("LEFT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?TM",splitted,ignore.case = T)+1])
        mat[count,13] <- toString(splitted[grep("RIGHT_[0-9]?[0-9]?[0-9]?[0-9]?[_]?TM",splitted,ignore.case = T)+1])
        #save the size of the amplicons
        mat[count,14] <- toString(splitted[grep("PRIMER_PRODUCT_SIZE[_]?[0-9]?[0-9]?[0-9]?[0-9]?",splitted,ignore.case = T)+1])
        #format to UCSC style
        mat[count,15] <- ucsc
        mat[count,16] <- paste0(as.character(ucsc_ranges@seqnames@values),":",ucsc_ranges@ranges@start+as.numeric(mat[count,6]),"-",ucsc_ranges@ranges@start+as.numeric(mat[count,6])-1+as.numeric(mat[count,7]))
        mat[count,17] <- paste0(as.character(ucsc_ranges@seqnames@values),":",ucsc_ranges@ranges@start+as.numeric(mat[count,8])+1-as.numeric(mat[count,9]),"-",ucsc_ranges@ranges@start+as.numeric(mat[count,8]))
        # if(length(grep("Minus",id))!=0){
        #   mat[count,15] <- paste0(as.character(ucsc_ranges@seqnames@values),":",ucsc_ranges@ranges@start,"-",ucsc_ranges@ranges@start+ucsc_ranges@ranges@width)
        #   end <- ucsc_ranges@ranges@start+ucsc_ranges@ranges@width
        #   mat[count,16] <- paste0(as.character(ucsc_ranges@seqnames@values),":",end-as.numeric(mat[count,6])-as.numeric(mat[count,7])-1,"-",end-as.numeric(mat[count,6])-1)
        #   mat[count,17] <- paste0(as.character(ucsc_ranges@seqnames@values),":",end-2-as.numeric(mat[count,8]),"-",end-as.numeric(mat[count,8])+as.numeric(mat[count,9])-2)
        # }
        mat[count,18] <- paste0(stringr::str_split_fixed(mat[count,16],"-",2)[1],"-",stringr::str_split_fixed(mat[count,17],"-",2)[2])
        #count number of CpG in Amplicon -> this function might be a bit slow
        ranges_cg <- mat[count,18]
        chr_loc <- stringr::str_split_fixed(ranges_cg,":",2)
        chr_loc <- stringr::str_split_fixed(chr_loc,"-",2)
        if(length(grep("Minus",id))!=0){
          mat[count,19] <- stringr::str_count(as.character(discoverSequence(chr_loc[1,1],
                                                                            as.numeric(chr_loc[2,2]),
                                                                            as.numeric(chr_loc[2,1]),"hg19"))
                                              ,"CG")}else{
                                                mat[count,19] <- stringr::str_count(as.character(discoverSequence(chr_loc[1,1],
                                                                                                                  as.numeric(chr_loc[2,1]),
                                                                                                                  as.numeric(chr_loc[2,2]),"hg19")),"CG")}
        count <- count+1

      }
      x <- i
    }
  }

  return(mat)
}

# Create a second strand --------------------------------------------------

secStrand <- function(x){
  save.names <- names(x)
  dna.rev <- lapply(x, MassArray::revComplement)
  names(dna.rev) <- paste0(save.names,"_Minus")
  return(dna.rev)
}



gviz_primer <- function(genome,chr,begin,end,primer_fwd.start,primer_fwd.length,
                        primer_rev.end,primer_rev.length,strand,primer=T){
  
  start_time <- Sys.time()
  
  #AmpliconDesign::gviz_primer("hg19","chr1",1565681,1566182,39,20,326,20,"+")
  #gviz_primer("hg19","chr19",43203328,43203689,63,20,324,23,"+")
  print("Please wait while R is calculating")
  # print(paste0(c(genome,chr,begin,end,primer_fwd.start,primer_fwd.length,
  #                primer_rev.end,primer_rev.length,strand)))

  
  ## DATABASE QUERY ##
  start_time_q <- Sys.time()
  my_db = DBI::dbConnect(RPostgreSQL::PostgreSQL(), user="postgres", password="postgres",
                         host="c010-shiny2", port=5432, dbname="testdb")
  on.exit(DBI::dbDisconnect(my_db))

  query_list <- list(
    "SNP"=as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM snp WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))")))),
    "CGI"=as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM cgi WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))")))),
    "genes"=as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM genes WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))")))),
    "repeats"=as.data.frame(dplyr::tbl(my_db,dplyr::sql(paste0("SELECT * FROM repeats WHERE genome='",genome,"' AND chromosome='",chr,"' AND ((start>=",begin," AND start<=",end,") OR (\"end\">=",begin," AND \"end\"<=",end,"))"))))
  )
  
  end_time_q <- Sys.time()
  message(paste0("Time for the Postgres DB Query:",c(end_time_q - start_time_q)))
  
  
  #call genome of interest
  # if (genome=="hg19"){
  #   genome.track <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  # 
  # } else if(genome=="hg38"){
  #   genome.track <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  # 
  # } else if(genome=="mm10"){
  #   genome.track <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  # }
  
  #create a sequence track
  #sTrack <- Gviz::SequenceTrack(Biostrings::DNAStringSet(discoverSequence(chr,begin,end,genome)))
  #create a genome axis track
  gtrack <- Gviz::GenomeAxisTrack()

  #find CpGs within Track
  seq.cg <- Biostrings::DNAString(discoverSequence(chr,begin,end,genome))
  #names(seq.cg) <- chr1
  #sTrack <- Gviz::SequenceTrack(seq.cg)
  
  if(strand=="-"){
    cg.df <- Biostrings::matchPattern("CG",Biostrings::reverseComplement(seq.cg))}else{
      cg.df <- Biostrings::matchPattern("CG",seq.cg)
    }
  cg.gr <- GenomicRanges::GRanges(chr,strand = strand,ranges = IRanges::IRanges(cg.df@ranges@start+begin,width=2))
  cg.plot <- Gviz::AnnotationTrack(cg.gr,name="CpG",shape="box",fill="grey")

  #create a plot for the primers
  if(primer==T){
  primer_start <- c(primer_fwd.start+2,primer_rev.end+3-primer_rev.length)
  primer_length <- c(primer_fwd.length,primer_rev.length)
  primer.gr <- GenomicRanges::GRanges(chr,strand = strand,ranges = IRanges::IRanges(primer_start+begin,width=primer_length))
  primer.plot <-Gviz:: AnnotationTrack(primer.gr,name="Primer",shape="box",fill="brown4")
  }
  
  #SNP track
  if(nrow(query_list$SNP)!=0){
    snp.plot <- Gviz::AnnotationTrack(query_list$SNP,name="SNP",shape="box",fill="palegreen4")
  }else{
    snp.plot <-  Gviz::AnnotationTrack(NULL,name="SNP",showID=T)  
  }
  
  #Repeat track
  if(nrow(query_list$repeats)!=0){
   rep.plot <- Gviz::AnnotationTrack(query_list$repeats,name="Repeats",shape="box",fill="sienna1")
  }else{
   rep.plot <-  Gviz::AnnotationTrack(NULL,name="Repeats",showID=T)  
  }
  
  #SNP track
  if(nrow(query_list$CGI)!=0){
    cgi.plot <- Gviz::AnnotationTrack(query_list$CGI,name="CGI",shape="box",fill="navyblue")
  }else{
    cgi.plot <-  Gviz::AnnotationTrack(NULL,name="CGI",showID=T)  
  }

  #add transcripts
  if(nrow(query_list$genes)!=0){
    genetrack <-  Gviz::GeneRegionTrack(query_list$genes,name="Transcripts",showID=T)
  }else{
    genetrack <-  Gviz::AnnotationTrack(NULL,name="Transcripts",showID=T)  
    }

  Gviz::displayPars(snp.plot) <- list(size=5)
  if(primer==T){
  Gviz::displayPars(primer.plot) <- list(size=5)
  }
  
  Gviz::displayPars(cg.plot) <- list(size=5)
  Gviz::displayPars(cgi.plot) <- list(size=5)
  Gviz::displayPars(rep.plot) <- list(size=5)
  #Gviz::displayPars(genetrack) <- list(fontsize=10)
  #gene <- coMET::knownGenes_UCSC(genome, chr, begin, end, showId=TRUE)

  end_time <- Sys.time()
  message(paste0("Time for the plot calculation:",c(end_time - start_time)))
    
  #return(Gviz::plotTracks(list(sTrack,gtrack,cg.plot,primer.plot,genetrack,cgi.plot,rep.plot,snp.plot),from=begin,to=end,cex=0.8,fontcolor.exon=20))
  if(primer==T){
  return(Gviz::plotTracks(list(gtrack,cg.plot,primer.plot,genetrack,cgi.plot,rep.plot,snp.plot),from=begin,to=end,cex=0.8,fontcolor.exon=20))
  }else{
    return(Gviz::plotTracks(list(gtrack,cg.plot,genetrack,cgi.plot,rep.plot,snp.plot),from=begin,to=end,cex=0.8,fontcolor.exon=20))  
  }
  end_time2 <- Sys.time()
  message(paste0("Time for the plot printing:",c(end_time2 - start_time)))
}




###################################################################
#                                                                 #
#                                                                 #
#                                                                 #
#                 PrimerBlast Functions                           #
#                                                                 #
#                                                                 #
#                                                                 #
###################################################################

primer_bs_blast <- function(primer_vec,
                            path_to_bowtie="/srv/shiny-server/software/bowtie-1.2.2/bowtie",
                            path_to_genome="/var/ressources/hg19_bowtie_ct_ga_indexed.fa",
                            k=10){

# bowtie /C010-projects/Maxi_Schönung/Computational/genomes/genomes_bs_blast/hg19_bowtie_ct_ga_indexed.fa -c CACAaCCTCCCCAAaTaCTa -k 20 --large-index

  align <<- tempfile()
  system(paste(path_to_bowtie,path_to_genome,"-c",paste(primer_vec,collapse = ","),"-k",k,"--large-index >",align))
  primer_align <- read.delim(align,header = F)[,-c(1,6,7)]
  colnames(primer_align) <- c("Strand","Chromosome","Start","Sequence","Mismatches")
  return(primer_align)
}


# AmpBS-Seq Analysis ---------------------------------------------------------------


#' load_amp
#'
#' @param files
#' @param meta
#' @import BiSeq
#'
#' @return
#' @export
#'
#' @examples
load_amp <- function(files,meta,files_orig){
  meta <- read.delim(meta,row.names = 1)
  #order.df <- do.call(rbind,stringr::str_split(files,"/"))
  #order.df[,ncol(order.df)]
  BiSeq::readBismark(files[match(rownames(meta),files_orig)],colData=meta)
}


#' qc_plots
#'
#' @param set
#'
#' @return
#' @export
#' @import BiSeq
#'
#' @examples
qc_plots <- function(set){

  id <- as.character(set@colData$UPN)
  id <- factor(id, levels = id)

  #### total reads
  total.reads <- BiSeq::totalReads(set)
  total.reads <- data.frame("Sample"=id,"Sites"=colSums(total.reads))
  p1 <- ggplot(total.reads,aes(Sample,log10(Sites)))+
    geom_bar(stat="identity")+
    theme_minimal()+
    ylab("log10(Total Reads per Sample)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #### Detected CpGs
  cg.cov <- BiSeq::covStatistics(set)
  cg.cov <- data.frame("Sample"=id,"CpG"=cg.cov$Covered_CpG_sites)
  p2 <- ggplot(cg.cov,aes(Sample,CpG))+
    geom_bar(stat="identity")+
    theme_minimal()+
    ylab("Covered CpG sites")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


  #### Numer of reads per sample versus detected CpG sites
  reads_cg <- cbind(total.reads,"CpG"=cg.cov$CpG)
  p3 <- ggplot(reads_cg,aes(log10(Sites),CpG))+
    geom_point()+
    theme_minimal()+
    geom_smooth(method="lm")+
    ylab("Covered CpG Sites")+
    xlab("log10(Reads)")

  #### Coverage per CpG
  cov_cg <- BiSeq::totalReads(set)
  colnames(cov_cg) <- as.character(id)
  cov_cg.melt <- reshape2::melt(cov_cg)
  cov_cg.melt <- cov_cg.melt[cov_cg.melt$value>0,]
  p4 <- ggplot(cov_cg.melt,aes(Var2,log10(value)))+
    geom_boxplot()+
    theme_minimal()+
    xlab("Sample")+
    ylab("log10(Coverage per CpG)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  return(list(p1,p2,p3,p4))
}


#' clustering
#'
#' @param set
#' @param by
#'
#' @return
#' @export
#'
#' @examples
clustering <- function(set,by){
  # get beta matrix
  meth.clean.rel <- BiSeq::rawToRel(set)
  meth.cluster <- na.omit(SummarizedExperiment::assay(meth.clean.rel))
  #calculate pca
  pca <- prcomp(t(meth.cluster))
  variance <- (pca$sdev)^2 /sum(pca$sdev^2)
  p1 <- ggplot(as.data.frame(pca$x),aes(PC1,PC2))+
    geom_point()+
    theme_bw()+
    xlab(paste0("PC1 (",round(variance[1]*100,2),"% of Variance)"))+
    ylab(paste0("PC2 (",round(variance[2]*100,2),"% of Variance)"))+
    theme(legend.title = element_blank())
  if(by=="zero"){
    print(p1)
  }else{
    print(p1+geom_point(aes(color=set@colData[,by])))
  }
}


#' amp_heatmap
#'
#' @param set
#' @param col_heatmap
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
amp_heatmap <- function(set,col_heatmap,...){
  meth.clean.rel <- BiSeq::rawToRel(set)
  ids <- as.data.frame(SummarizedExperiment::rowRanges(set))
  meth.cluster <- SummarizedExperiment::assay(meth.clean.rel)
  rownames(meth.cluster) <- paste0(ids$seqnames,":",ids$start)
  meth.cluster <- na.omit(meth.cluster)
  pheatmap::pheatmap(meth.cluster,
                     clustering_distance_cols = "manhattan",
                     clustering_distance_rows = "manhattan",
                     clustering_method = "ward.D2",
                     annotation_col = as.data.frame(set@colData[col_heatmap]),...)
}


region_heatmap <- function(set,bed_ranges,co=T,co_num){
  set1 <- IRanges::subsetByOverlaps(set,GenomicRanges::GRanges(bed_ranges))
  if(co==T){
    set1 <- BiSeq::filterByCov(set1,as.numeric(co_num),global=FALSE)
  }
  meth.set1 <- BiSeq::rawToRel(set1)
  meth.set1 <- SummarizedExperiment::assay(meth.set1)
  ids <- as.data.frame(SummarizedExperiment::rowRanges(set1))
  rownames(meth.set1) <- paste0(ids$seqnames,":",ids$start)
  pheatmap::pheatmap(meth.set1,cluster_rows = F,cluster_cols = F)
}
