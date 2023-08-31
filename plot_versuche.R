setwd(cT_path)
ct_files <- list.files(path=cT_path,pattern="*.ct",full.names = FALSE)

# go through each ct file in folder and plot them
for (i in 1:length(ct_files)){

  #extract virusname from ct_files for plot_title.
  vnames <- new_variable <- gsub("\\.ct$", "", ct_files[i])
  print(ct_files)

  foldplot <- aptPlotCT(ct_files[i],nt=TRUE,labTF=TRUE, main=vnames) # you have to setwd to the folder containing the ct-files.
  print(foldplot)

  p#lot_path2 <- file.path(folder_path_, "Alle_Foldplots") # check for existing folder and creat if not there.

  if (!dir.exists(plot_path2)) {
    dir.create(plot_path2, recursive = TRUE)
  }

  # Export the plot as a PNG file
  #plot_filename <- paste0(plot_path2,"/","foldplot_",vnames,".png")
  #png(filename=plot_filename, width = 480, height = 480)
  #dev.off()  # Close the PNG device
}

aptPlotCT("3_UTR of Eriocheir_sinensis_ronivirus.ct")

aptPlotCT("3_UTR of Eriocheir_sinensis_ronivirus.ct",nt=TRUE,labTF=TRUE, main="3_UTR of Eriocheir_sinensis_ronivirus",add=TRUE,dp=0.5)

test <- loadCt("3_UTR of Eriocheir_sinensis_ronivirus.ct")
RNAstrPlot(test)


setwd("A:/Praktikum_Chris/R/code")
