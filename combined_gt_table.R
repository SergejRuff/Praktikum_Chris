all_data <- list()


for (i in 1:length(sequences)){

  ##############
  # preprocessing
  ###############

  # extract specific fasta file from list
  sequence_ <- sequences[[i]]

  # extract the nucleotide sequence
  seq <- sequence_[1]

  seq<- toupper(seq)

  # get name of taxa
  NameofTaxon <- getName(sequence_)

  # find the ORFs in + direction with length of 300 or more. Using ORFik-package
  orfs<- findORFs(seq, minimumLength = 98,startCodon = "ATG")      # !!!!!!!! might need to change length later. 248 for 750
  # !!!! Startcodon changes results. Default uses alternative codons as well.


  # Convert the Data to a data frame
  orf_df <- as.data.frame(orfs@unlistData)

  # extract last position of 5´-UTR
  firstposition <- min(orf_df$start)-1

  # extract first position of 3´-UTR
  lastposition <- max(orf_df$end)+1


  # extract 5´- and 3´-UTR
  five_UTR <- substring(seq, first = 1, last=firstposition)
  three_UTR <- substring(seq, first =lastposition )

  #RNAstructure requires nucleotides to be upper case
  five_UTR <- toupper(five_UTR)

  three_UTR <- toupper(three_UTR)


  #cat(paste("for virus",NameofTaxon))
  #print(orfs) # print info about available ORFs.

  # Create a data frame with the information for this run
  data <- data.frame(
    "Virus" = NameofTaxon,
    "UTRs" = c("5´-UTR Position", "3´-UTR Position"),
    "position" = c(paste("1-",firstposition), paste(lastposition,"-",nchar(seq))),
    "length"= c(firstposition, nchar(seq) - lastposition + 1)
  )

  # Add the data frame to the list
  all_data[[i]] <- data


}

# Combine all data frames into one
combined_data <- do.call(rbind, all_data)

# Create a gt table using the combined data
combined_table <- gt(combined_data) %>%
  tab_header(title = "UTR Positions and Lengths") %>%
  opt_align_table_header(align = "left")

# Print the combined table
print(combined_table)
