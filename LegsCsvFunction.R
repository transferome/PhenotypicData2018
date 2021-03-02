
# Execute two lines below, readline command and then desired filepath beneath, this converts the forward slash to a double forward slash
tempPath <- readline()
Z:\SelectionSimulans2018\Generation_6b
# Run this command which will convert double forward slashes to single backslashes, L should be the desired filepath, but with forward slashes converted to backslashes
Path <- gsub("\\\\","/",tempPath)
Path




Legs  <- function(Path, Type) {
  
  
######################################################
##### Find all folders that contain leg files ########
######################################################

## Lists first leg file plus its extension
folders  <- list.files(path = Path, full.names = TRUE, recursive = TRUE, pattern = ".*0001\\.legs$")
## Lists the Leg folder and it's path that contain the leg files listed in folders
tpaths   <-  dirname(folders)
## Lists all folders that contain Leg folders that contain leg files
tpaths1  <- dirname(tpaths)

## Gives the filepath and file extension of any Legs.csv file
folders1    <- list.files(path = Path, full.names = TRUE, recursive = TRUE, pattern = ".*Legs\\.csv$")
## This removes the extension (the filename) so all that is left is the folders contained Legs.csv files
csvpaths    <- dirname(folders1)

## This finds which folders that contain a Legs folder with .legs files, do not contain Legs.csv's
tpaths2 <- setdiff(tpaths1,csvpaths)

## The paths we care about are ones specified above, but need /Legs added to them for the next step of the function
paths   <- paste(tpaths2, "/Legs", sep = "")

######################################################
######################################################


######################################################
############# Leg Path Loop ##########################
######################################################

######
###### This will run another loop for every file path found to contain leg files

for (h in 1:length(paths)) {
    x  <- paths[h]
    
    #### Want to go back a folder for eventually outputting the .csv file so
    z  <- gsub("(.*)/.*","\\1",paths[h])

##### Lists the entire file path of every leg file in path[h]
files  <- list.files(path = x, full.names = TRUE, pattern = "\\.legs$")
##### Lists only the filename of every leg file in path[h]
files2 <- list.files(path = x, full.names = FALSE, pattern = "\\.legs$")

##### A dummy vector to hold extracted information
hold   <- rep(NA, 73)
##### A dummy matrix that will store the hold vector for every legs file in path[h]
dm <- matrix(ncol = length(hold), nrow = length(files))


########################################################
########### Leg Data Loop ##############################
########################################################

for (i in 1:length(files)) {
  
  #### This loads in a leg data file
  data <- dget(files[i], keep.source = FALSE )
  
  #### Name of Leg file
  legname <- files2[i]
  hold[1] <- legname
  
  #### Image name which will be unique to inividuals
  indiv <- data$ascEntry$imageName
  hold[2] <- indiv
  
  #### This is to find the imaging station
  #### Remove everything up to including the  "_" in imagename
  tempM    <- gsub(".*_","",indiv)
  hold[70] <- substr(tempM, 1 , 1)
  
  #### This is for the imager
  hold[71] <- data$ascEntry$processorName
  
  #### Extracts Date
  hold[72] <-data$ascEntry$date
  
  #### Extracts Time
  hold[73] <-data$ascEntry$time
  
  #### This gives the genotype
  geno  <- data$ascEntry$genotype
  hold[3] <- geno
  
  #### This gives the sex
  bisex <- data$ascEntry$sex
  hold[4] <- bisex
  
  #################################################################
  ######## If statements to handle filenames ######################
  ######## This handles blocks,rnai, species, and selection   #####
  #################################################################
  
  #### This gets the file name
  a <- data$ascEntry$filename
  ###### If Type argument is not specified, it equals 1
  ###### Else Type equals whatever is specified, which should be a 2 for species
  ###### E.g. to a run a block Legs(Path), to run species Legs(Path, 2)
  
  if (missing(Type)) {
    type <- 1
  } else {
    type <- Type 
  }
  
  ###### If type = 1, the block number will be extracted from the filename argument
  ###### If type = 2, for the species routine than it will extract the species from filename
  if (type == 1) {
  blocknum <- as.numeric(strsplit(a, "\\D+")[[1]][-1])
  hold[5] <- blocknum
  } 
  if (type == 2) {
  species <- sub(".asc.*","", a)
  hold[5] <- species
  }
  if (type == 3) {
  rnai    <- sub(".asc.*","", a)
  hold[5] <- rnai
  }
  if (type == 4) {
  a <- basename(a)
  hold[5] <- sub(".asc.*","",a)
  }
  ###############################################
  ###### End of if statements for filename ######
  ###############################################
  
  #### This gets the scale
  Scale <- data$ascEntry$scale
  hold[6] <- Scale
  
  ##############################################
  ####### Beginning of if, else statement ######
  #######        for error status         ######
  ##############################################
  
  
  if (data$results$errorStatus == "none") {
    
    #############################################
    #############################################
    ###### Find x,y values for each segment #####
    #############################################
    #############################################
    
    ##################
    ###### Front #####
    ##################
    
    ### Base points for front leg
    b1x <- data$ascEntry$legMarks[[1]][1]
    b1y <- data$ascEntry$legMarks[[1]][2]
    hold[7] <- b1x
    hold[8] <- b1y
    
    ### x and y coordinate vectors for front leg
    flx <- data$results$legAxes[[1]][[1]]
    fly <- data$results$legAxes[[1]][[2]]
    
    ### Knee coordinates for front leg
    k1x <- data$results$kneeCoords[[1]][1]
    k1y <- data$results$kneeCoords[[1]][2]
    hold[9]  <- k1x
    hold[10] <- k1y
    
    ### Ankle coordinates for front leg
    a1x <- data$results$ankleCoords[[1]][1]
    a1y <- data$results$ankleCoords[[1]][2]
    hold[11] <- a1x
    hold[12] <- a1y
    
    ### End points for front leg
    e1x <- data$results$legEndCoords[[1]][1]
    e1y <- data$results$legEndCoords[[1]][2]
    hold[13] <- e1x
    hold[14] <- e1y
    
    ### Find what element of f1x matches with the Knee x coordinate
    
    f1  <- which(flx == k1x)
    n1  <- which(fly == k1y)
    w1  <- which(f1 %in% n1)
    f1  <- f1[w1]
    ### Femur goes from 1:f1
    
    ### Find what element of flx matches with the Ankle x coordinate
    f2  <- which(flx == a1x)
    n2  <- which(fly == a1y)
    w2  <- which(f2 %in% n2)
    f2  <- f2[w2]
    ### Tibia goes from f1:f2
    
    ### The Tarsi go from f2, to the last element of the vector
    f3  <- length(flx)
    ### Tarsi goes from f2:f3
    
    ##################
    ##### Middle #####
    ##################
    
    ### Base points for middle leg
    b2x <- data$ascEntry$legMarks[[2]][1]
    b2y <- data$ascEntry$legMarks[[2]][2]
    hold[15] <- b2x
    hold[16] <- b2y
    
    ### x and y coordinates for middle leg
    mlx <- data$results$legAxes[[2]][[1]]
    mly <- data$results$legAxes[[2]][[2]]
    
    ### Knee coordinates for middle leg
    k2x <-data$results$kneeCoords[[2]][1]
    k2y <-data$results$kneeCoords[[2]][2]
    hold[17] <- k2x
    hold[18] <- k2y
    
    ### Ankle coordinates for middle leg
    a2x <- data$results$ankleCoords[[2]][1]
    a2y <- data$results$ankleCoords[[2]][2]
    hold[19] <- a2x
    hold[20] <- a2y
    
    ### End points for middle leg
    e2x <- data$results$legEndCoords[[2]][1]
    e2y <- data$results$legEndCoords[[2]][2]
    hold[21] <- e2x
    hold[22] <- e2y
    
    ### Find what element of m1x matches with the Knee x coordinate
    m1  <- which(mlx == k2x)
    o1  <- which(mly == k2y)
    v1  <- which(m1 %in% o1)
    m1  <- m1[v1]
    ### Femur goes from 1:m1
    
    ### Find what element of flx matches with the Ankle x coordinate
    m2  <- which(mlx == a2x)
    o2  <- which(mly == a2y)
    v2  <- which(m2 %in% o2)
    m2  <- m2[v2]
    ### Tibia goes from m1:m2
    
    ### The Tarsi go from f2, to the last element of the vector
    m3  <- length(mlx)
    ### Tarsi goes from m2:m3
    
    
    ##################
    ##### Back #######
    ##################
    
    ### Base points for middle leg
    b3x <- data$ascEntry$legMarks[[3]][1]
    b3y <- data$ascEntry$legMarks[[3]][2]
    hold[23] <- b3x
    hold[24] <- b3y
    
    ### x and y coordinate for back leg
    blx <- data$results$legAxes[[3]][[1]]
    bly <- data$results$legAxes[[3]][[2]]
    
    # Knee coordinates for back leg
    k3x <- data$results$kneeCoords[[3]][1]
    k3y <- data$results$kneeCoords[[3]][2]
    hold[25] <- k3x
    hold[26] <- k3y
    
    # Ankle coordinates for back leg
    a3x <- data$results$ankleCoords[[3]][1]
    a3y <- data$results$ankleCoords[[3]][2]
    hold[27] <- a3x
    hold[28] <- a3y
    
    ### End points for middle leg
    e3x <- data$results$legEndCoords[[3]][1]
    e3y <- data$results$legEndCoords[[3]][2]
    hold[29] <- e3x
    hold[30] <- e3y
    
    ### Find what element of m1x matches with the Knee x coordinate
    b1  <- which(blx == k3x)
    z1  <- which(bly == k3y)
    zz1 <- which(b1 %in% z1)
    b1  <- b1[zz1]
    ### Femur goes from 1:b1
    
    ### Find what element of flx matches with the Ankle x coordinate
    b2  <- which(blx == a3x)
    z2  <- which(bly == a3y)
    zz2 <- which(b2 %in% z2)
    b2  <- b2[zz2]
    ### Tibia goes from b1:b2
    
    ### The Tarsi go from f2, to the last element of the vector
    b3  <- length(blx)
    ### Tarsi goes from b2:b3
    
    ########################################
    ########################################
    ### End of segment elements ############
    ########################################
    ########################################
    
    
    ########################################
    ##### Segment Lengths Point by Point ###
    ########################################
    ############### Commented out.  This subsets the x,y coordinate vectors by leg segment, needed for calculating point to point distances.
    
    ######## Front Leg Femur x,y
    #d1x <- flx[1:f1]
    #d1y <- fly[1:f1]
    ##### Front leg Tibia x,y
    #d2x <- flx[f1:f2]
    #d2y <- fly[f1:f2]
    ##### Front leg Tarsi x,y
    #d3x <- flx[f2:f3]
    #d3y <- fly[f2:f3]
    
    ######## Middle Leg Femur x,y
    #d4x <- mlx[1:m1]
    #d4y <- mly[1:m1]
    ##### Middle leg Tibia x,y
    #d5x <- mlx[m1:m2]
    #d5y <- mly[m1:m2]
    ##### Middle leg Tarsi x,y
    #d6x <- mlx[m2:m3]
    #d6y <- mly[m2:m3]
    
    ######## Back Leg Femur x,y
    #d7x <- blx[1:b1]
    #d7y <- bly[1:b1]
    ##### Back leg Tibia x,y
    #d8x <- blx[b1:b2]
    #d8y <- bly[b1:b2]
    ##### Back leg Tarsi x,y
    #d9x <- blx[b2:b3]
    #d9y <- bly[b2:b3]
    
    #xlist <- list( d1x, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x)
    #ylist <- list( d1y, d2y, d3y, d4y, d5y, d6y, d7y, d8y, d9y)
    
    ########################################
    ######################################## End of x,y coordinates listed by segment
    ########################################
    
    ########################################
    ########################################
    ######### Mean, Max, Min Widths ########
    ########################################
    ########################################
    
    ##### front leg width
    flw <- data$results$legWidth[[1]]
    
    ##### Front leg Femur width
    ffw  <- flw[1:f1]
    ##### Front leg Tibia width
    ftw  <- flw[f1:f2]
    ##### Front leg Tarsi width
    ftsw <- flw[f2:f3]
    
    hold[31]  <- max(ffw)
    hold[32]  <- min(ffw)
    hold[33] <- mean(ffw)
    
    hold[34]  <- max(ftw)
    hold[35]  <- min(ftw)
    hold[36] <- mean(ftw)
    
    hold[37]  <- max(ftsw)
    hold[38]  <- min(ftsw)
    hold[39] <- mean(ftsw)
    
    
    
    ##### middle leg width
    mlw <- data$results$legWidth[[2]]
    
    ##### middle leg Femur width
    mfw  <- mlw[1:m1]
    ##### middle leg Tibia width
    mtw  <- mlw[m1:m2]
    ##### middle leg Tarsi width
    mtsw <- mlw[m2:m3]
    
    hold[40]  <- max(mfw)
    hold[41]  <- min(mfw)
    hold[42] <- mean(mfw)
    
    hold[43]  <- max(mtw)
    hold[44]  <- min(mtw)
    hold[45] <- mean(mtw)
    
    hold[46]  <- max(mtsw)
    hold[47]  <- min(mtsw)
    hold[48] <- mean(mtsw)
    
    ##### back leg width
    blw <- data$results$legWidth[[3]]
    
    ##### back leg Femur width
    bfw  <- blw[1:b1]
    ##### back leg Tibia width
    btw  <- blw[b1:b2]
    ##### back leg Tarsi width
    btsw <- blw[b2:b3]
    
    hold[49]  <- max(bfw)
    hold[50]  <- min(bfw)
    hold[51] <- mean(bfw)
    
    hold[52]  <- max(btw)
    hold[53]  <- min(btw)
    hold[54] <- mean(btw)
    
    hold[55]  <- max(btsw)
    hold[56]  <- min(btsw)
    hold[57] <- mean(btsw)
    
    ########################################
    ########################################
    ########### End of Widths ##############
    ########################################
    ########################################
    
    #############################################
    ### Calculating point to point distances ####
    #############################################
    
    ################# Commented out
    #disth <- rep(NA, 9)
    
    #for (j in 1:9) {
    #dx   <- xlist[[j]]
    #dy   <- ylist[[j]]
    #dvec <- rep(NA, (length(dx) - 1))
    
    #for (k in 2:length(dx)) {
    #dist       <- sqrt( (dx[k] - dx[k-1] )^2 + ( dy[k] - dy[k-1] )^2 )
    #dvec[k-1]  <- dist
    #}
    #disth[j] <- sum(dvec)
    #}
    #disth
    #hold[58] <- disth[1]
    #hold[59] <- disth[2]
    #hold[60] <- disth[3]
    #hold[61] <- disth[4]
    #hold[62] <- disth[5]
    #hold[63] <- disth[6]
    #hold[64] <- disth[7]
    #hold[65] <- disth[8]
    #hold[66] <- disth[9]
    
    ########################################
    ######### Bills Lengths ################
    ########################################
    #Total lengths
    bt1 <- data$results$legLengths[1]
    bt2 <- data$results$legLengths[2]
    bt3 <- data$results$legLengths[3]
    
    
    # Bill's segment lengths
    # Length of front femur
    bd1 <- data$results$kneePositions[1]
    # Length of front tibia
    bd2 <- data$results$anklePositions[1] - data$results$kneePositions[1]
    # Length of front tarsal segments
    bd3 <- bt1 - data$results$anklePositions[1]
    # Legnth of middle femur
    bd4 <- data$results$kneePositions[2]
    # Length of middle tibia
    bd5 <- data$results$anklePositions[2] - data$results$kneePositions[2]
    # Length of middle tarsal segments
    bd6 <- bt2 - data$results$anklePositions[2]
    # Length of back femur
    bd7 <- data$results$kneePositions[3]
    # Length of back tibia
    bd8 <- data$results$anklePositions[3] - data$results$kneePositions[3]
    # Length of back tarsal segments
    bd9 <- bt3 - data$results$anklePositions[3]
    
    
    ###### When the point by point distances are not commented out this goes from 67 - 78
    hold[58] <- bd1
    hold[59] <- bd2
    hold[60] <- bd3
    hold[61] <- bd4
    hold[62] <- bd5
    hold[63] <- bd6
    hold[64] <- bd7
    hold[65] <- bd8
    hold[66] <- bd9
    hold[67] <- bt1
    hold[68] <- bt2
    hold[69] <- bt3
    
    ##################################
    ###### End of Bill's Lengths #####
    ##################################
    
    
  } else {
    ## When the error status is anything other than "none", it is a bad spline
    ## No data can be collected and so elements 7:69 of the hold vector will be
    ## Filled with NA's
    hold[7:69]  <- NA
  }
  ##### Fill in the i'th row of the dummy matrix with the current hold vector
  ##### Every row of dm will be filled with information from each leg file from the current folder
  dm[i,]  <- hold
  
}
 #############################################
 ######## End of Leg Data Loop ###############
 #############################################

#### If else statement for how to name the columns of the csv file
#### Difference is between block and species naming
if (type == 1) {
  newcol <- c("legfile","indi", "geno", "sex","block", "scale", "b1x","b1y","k1x","k1y","a1x","a1y","e1x","e1y","b2x","b2y","k2x","k2y","a2x","a2y","e2x","e2y","b3x","b3y","k3x","k3y","a3x","a3y","e3x","e3y","maxF1","minF1","meanF1","maxT1","minT1","meanT1","maxTs1","minTs1","meanTs1","maxF2","minF2","meanF2","maxT2","minT2","meanT2","maxTs2","minTs2","meanTs2","maxF3","minF3","meanF3","maxT3","minT3","meanT3","maxTs3","minTs3","meanTs3","F1len","T1len","Ts1len","F2len","T2len","Ts2len","F3len","T3len","Ts3len","len1","len2","len3","station","imager","date", "time")
} 
if (type == 2) {
  newcol <- c("legfile","indi", "geno", "sex","species", "scale", "b1x","b1y","k1x","k1y","a1x","a1y","e1x","e1y","b2x","b2y","k2x","k2y","a2x","a2y","e2x","e2y","b3x","b3y","k3x","k3y","a3x","a3y","e3x","e3y","maxF1","minF1","meanF1","maxT1","minT1","meanT1","maxTs1","minTs1","meanTs1","maxF2","minF2","meanF2","maxT2","minT2","meanT2","maxTs2","minTs2","meanTs2","maxF3","minF3","meanF3","maxT3","minT3","meanT3","maxTs3","minTs3","meanTs3","F1len","T1len","Ts1len","F2len","T2len","Ts2len","F3len","T3len","Ts3len","len1","len2","len3","station","imager","date", "time")
}
if (type == 3) {
  newcol <- c("legfile","indi", "geno", "sex","rnai", "scale", "b1x","b1y","k1x","k1y","a1x","a1y","e1x","e1y","b2x","b2y","k2x","k2y","a2x","a2y","e2x","e2y","b3x","b3y","k3x","k3y","a3x","a3y","e3x","e3y","maxF1","minF1","meanF1","maxT1","minT1","meanT1","maxTs1","minTs1","meanTs1","maxF2","minF2","meanF2","maxT2","minT2","meanT2","maxTs2","minTs2","meanTs2","maxF3","minF3","meanF3","maxT3","minT3","meanT3","maxTs3","minTs3","meanTs3","F1len","T1len","Ts1len","F2len","T2len","Ts2len","F3len","T3len","Ts3len","len1","len2","len3","station","imager","date", "time")
}
if (type == 4) {
  newcol <- c("legfile","indi", "geno", "sex","Generation", "scale", "b1x","b1y","k1x","k1y","a1x","a1y","e1x","e1y","b2x","b2y","k2x","k2y","a2x","a2y","e2x","e2y","b3x","b3y","k3x","k3y","a3x","a3y","e3x","e3y","maxF1","minF1","meanF1","maxT1","minT1","meanT1","maxTs1","minTs1","meanTs1","maxF2","minF2","meanF2","maxT2","minT2","meanT2","maxTs2","minTs2","meanTs2","maxF3","minF3","meanF3","maxT3","minT3","meanT3","maxTs3","minTs3","meanTs3","F1len","T1len","Ts1len","F2len","T2len","Ts2len","F3len","T3len","Ts3len","len1","len2","len3","station","imager","date", "time")
}
#### End of if else for naming columns

#### Renames the columns of the dummy matrix according to the if else statements above
colnames(dm) <- newcol

#### Prints the current block or species name.  This is for possible troubleshooting.
#### Blocks/species that were successful should have there name printed before error
print(dm[1,5])

#### If else statements for the .csv file name
#### if names it Block#Legs.csv
#### else names it SpeciesLegs.csv, RNaiNameLegs.csv, or Generation#Legs.csv
if (type == 1) {
  f            <- paste("/Block",dm[1,5],"Legs",".csv", sep = "")
} else {
  f            <- paste("/",dm[1,5],"Legs",".csv", sep = "")
}
#### End of if else for naming the .csv output files

#### The output will be the combination of z, which is the folder outside of the Legs folder
#### And f, which is the name of the .csv file specified above
output       <- paste(z,f,sep = "" )

#### Writes the csv file
write.csv(dm, output, row.names = FALSE )
}
##################################################
####### End of Leg Path Loop #####################
##################################################
}

###### For Partial Diallel 
Legs ( Path )
###### For Species
Legs( Path, 2 )
###### For RNAi
Legs( Path, 3 )
###### For Selection
Legs( Path, 4 )

