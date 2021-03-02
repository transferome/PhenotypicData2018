########################## Selection for Mel and Simulans 2018 ##################################


select <- function(x,y,z) {

species <- 'Sim'
if (species == 'Sim') {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/'
} else {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/'
}


### This function will keep only the last n characters of a string
substrRight <- function(w,n) {
  substr(w, nchar(w)-n+1, nchar(w))
}

generation <- '6b'
replicate  <- substrRight(generation,1)

if (species == 'Sim' & replicate == 'a') {
  male_means   <- read.table( '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/SimMaleMeansA.csv', sep = ',', header = TRUE)
  female_means <- read.table( '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/SimFemaleMeansA.csv', sep = ',', header = TRUE)
}
if (species == 'Sim' & replicate == 'b') {
  male_means   <- read.table( '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/SimMaleMeansB.csv', sep = ',', header = TRUE)
  female_means <- read.table( '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/SimFemaleMeansB.csv', sep = ',', header = TRUE)
}
if (species == 'Mel' & replicate == 'a') {
  male_means   <- read.table( '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/MelMaleMeansA.csv', sep = ',', header = TRUE)
  female_means <- read.table( '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/MelFemaleMeansA.csv', sep = ',', header = TRUE)
}
if (species == 'Mel' & replicate == 'b') {
  male_means   <- read.table( '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/MelMaleMeansB.csv', sep = ',', header = TRUE)
  female_means <- read.table( '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/MelFemaleMeansB.csv', sep = ',', header = TRUE)
}



#legclassified file path and name
legpath <- paste0(par_dir,'Generation_',generation,'/Gen',generation,'legclassified.csv')
#cpr file path and name
cprpath <- paste0(par_dir,'Generation_',generation,'/Output_Landmarks+Veins+Outline.dat')
#wng file path and name
wngpath <- paste0(par_dir,'Generation_',generation,'/Wings/','WGen',generation,'.asc')


#### read in the 3 data files, leg data, cpr data, and wing data (wing data is used to align the original tif name into the cpr file, so that the cpr file and leg file can be merged by this unique ID)
legdat     <- read.table(legpath, sep=",", header = TRUE)
cprdat     <- read.table(cprpath, sep="\t", header = TRUE)
wngdat     <- read.table(wngpath, sep="\t", header = FALSE)

#### add .tif to the CPFile so that cpr data can be merged with the wing data file
cprdat$CPFile <- paste0(cprdat$CPFile,'.tif')

#### merge cprdat with wngdat by variables that they have in common
# V1 = the wing tif name, V6 = imager(same as perp), V7 = Date, V9 = Genotype, V10 = Sex
cpr_wng <- merge(cprdat, wngdat, by.x = c("CPFile", "Perp","Date","Tags","Sex"), by.y = c("V1","V6","V7","V9","V10"), all = TRUE)

#### These columns are all 1's and can be removed
column_check1 <- which(colnames(cpr_wng) == "V12") #118
column_check2 <- which(colnames(cpr_wng) == "V13") #119

if (column_check1 == 118 & column_check2 == 119) {
  print('First Data Merge Appears Okay')
} else {
  print ('WARNING! Data Merge Appears out of Order')
}

#### Remove columns 118, and 119
cpr_wng <- cpr_wng[ , - c(118,119)]


###merge by V14(original tif name, the unique tif ID), then V11(scale), Tags, Perp, Date, Sex
dat <- merge(cpr_wng, legdat, by.x = c('V14','V11','Tags','Perp','Sex'), by.y = c('indi','scale','geno','imager','sex'), all = TRUE)
colnames(dat)[colnames(dat) == 'V11'] <- 'orgscale'
colnames(dat)[colnames(dat) == 'V14'] <- 'individual'


### the individual ID can be used to get the vial number the fly is in
### make a temporary string that will remove the .tif from the individual ID
tempstr <- gsub("\\..*", "", dat$individual)


### Take the last three characters from the tempstr, and create the flynum vector
flynum <- substrRight(tempstr,3)

### Convert character vector into an integer
flynum <- as.integer(flynum)

### Add flynum to the data
dat$flynum <- flynum
### in the dat data file, V2,V3,V4,V5 are the original wing coordinates, and V8 is the time the image was taken


######Data merge check
column_check3 <- which(colnames(dat) == "Sex") #5
column_check4 <- which(colnames(dat) == "x32") #77
column_check5 <- which(colnames(dat) == "flynum") #193
if (column_check3 == 5 & column_check4 == 77 & column_check5 == 193) {
  print('Final Data Merge Appears Okay')
} else {
  print ('WARNING! Final Data Merge Appears out of Order')
}

temp_out <- paste0(par_dir,'Generation_',generation,'/',species,generation,'NAunremoved.csv')
######Write a csv file that has all of the NA's and Removed Legs still included
write.csv(dat, temp_out, row.names = FALSE)

### remove all NAs from the data set
dat <- na.omit(dat)

### remove all individuals that are categorized as a 1 in the removed column
dat <- dat[dat$removed != 1 , ]


male_dat   <- subset(dat, Sex == 'M')
female_dat <- subset(dat, Sex == 'F')

dat.list <- list(male_dat,female_dat)
ratios.list <- list()
wing_deviation.list <- list()

for (i in 1:length(dat.list)) {
  #### for first iteration the wing_means is from the male generation 0 file, for second iteration the wing_means is from the female generation 0 file
    if (i == 1) {
    wing_means <- male_means
    } else {
    wing_means <- female_means
    }
  #### name the male/female dat file as data
  data <- dat.list[[i]]
  
  #### double checks that the leg data is where it is supposed to be
  column_check6 <- which(colnames(data) == 'F1len') #172
  column_check7 <- which(colnames(data) == 'len3') #183
    if (column_check6 == 172 & column_check7 == 183) {
      print('Leg Data is in the Proper Columns')
    } else {
      print ('Warning! Leg Data Does Not Appear to be in the Correct Columns')
    }

  ### multiply the leg sgements lengths that are in pixels, by the scale to get the mm's
  temp_lengths <- data[, c(172:183)] * data$orgscale
  ###rename these length variables
  names(temp_lengths) <- c('F1','T1','Ts1','F2','T2','Ts2','F3','T3','Ts3','L1','L2','L3')
  ### calculate the ratio of each leg segment for each leg
  front_ratios <- as.data.frame(cbind(temp_lengths$F1/temp_lengths$L1,temp_lengths$T1/temp_lengths$L1, temp_lengths$Ts1/temp_lengths$L1 ))
  middle_ratios <- as.data.frame(cbind(temp_lengths$F2/temp_lengths$L2,temp_lengths$T2/temp_lengths$L2, temp_lengths$Ts2/temp_lengths$L2 ))
  back_ratios <- as.data.frame(cbind(temp_lengths$F3/temp_lengths$L3,temp_lengths$T3/temp_lengths$L3, temp_lengths$Ts3/temp_lengths$L3 ))
  ### column bind these ratio dataframes
  ratio_dat <- cbind(front_ratios, middle_ratios, back_ratios)
  ### name these as such
  names(ratio_dat) <- c('F1','T1','Ts1','F2','T2','Ts2','F3','T3','Ts3')

  ratios.list[[i]] <- ratio_dat


  ### check on taking out the wing coordinates
  column_check8 <- which(colnames(data) == 'x1')
  column_check9 <- which(colnames(data) == 'y49')
   if (column_check8 == 15 & column_check9 == 112) {
  print('Wing Data is in the Proper Columns')
   } else {
  print ('Wing Data Does Not Appear to be in the Correct Columns')
    }

  ### take out the wing coordinates from the dataframe
  wng_cor <- data[, c(15:112)]

  ### dummy matrix of zero that the for loop will fill in with the deviation of each individuals coordinate from the mean of that coordinate
  wng_dev_mat <- matrix(0, nrow = nrow(wng_cor), ncol = ncol(wng_cor)) 
  ### loop subtracts wng_cor column by that columns mean value
  for (j in 1:ncol(wng_cor)) {
   temp_vec <- wng_cor[, j] - wing_means[,j]
  wng_dev_mat[,j] <- temp_vec
   }

  ### wing square deviation matrix
  wng_sqr_dev_temp <- wng_dev_mat %*% t(wng_dev_mat)
  ### the diagonal of this matrix, is the sum of the squared deviations of each coordinate for each individual
  wng_deviations <- diag(wng_sqr_dev_temp)
  wing_deviation.list[[i]] <- wng_deviations

}

male_dat <- cbind(male_dat, ratios.list[[1]])
female_dat <- cbind(female_dat, ratios.list[[2]])
male_dat$wngdev   <- wing_deviation.list[[1]]
female_dat$wngdev <- wing_deviation.list[[2]]

### This is the selection vector for legs
leg_sel_vec <- c(0.2357, 0.2357, -0.47140, 0.2357, 0.2357, -0.47140, 0, 0, 0)
leg_sel_vec1 <- (leg_sel_vec + 0.018) #balances negative and positive scores
len <- sqrt(t(leg_sel_vec1) %*% leg_sel_vec1) # this is the original length
leg_sel_vec1 <- as.vector(leg_sel_vec1)/as.vector(len) #standardizes to unit length 1


dat.list <- list(male_dat, female_dat)
sel_index.list <- list()
for (i in 1:length(dat.list)) {
  data <- dat.list[[i]]

  ### check on taking out the wing coordinates
  column_check10 <- which(colnames(data) == 'wngdev') #203
  column_check11 <- which(colnames(data) == 'F1') #194
    if (column_check10 == 203 & column_check11 == 194) {
      print('Data Still In proper columns')
    } else {
      print ('WARNING! Columns Out of Order')
    }
  
  sel_index <- as.matrix(data[ , c(194:202)]) %*% leg_sel_vec1
  
  if (i ==1) {
    hist(sel_index, main = 'Males')
  } else {
    hist(sel_index, main = 'Females')
  }
  
  print('Check Histograms for Weirdness')
  sel_index.list[[i]] <- sel_index
  
}

male_dat$dir_index   <- sel_index.list[[1]]
female_dat$dir_index <- sel_index.list[[2]]


m_stabindex <- ifelse(male_dat$Tags %in% c('up1','up2'), (male_dat$dir_index - (z*male_dat$wngdev)), (male_dat$dir_index + (z*male_dat$wngdev)))
f_stabindex <- ifelse(female_dat$Tags %in% c('up1','up2'), (female_dat$dir_index - (z*female_dat$wngdev)), (female_dat$dir_index + (z*female_dat$wngdev)))

male_dat$stab_index   <- m_stabindex
female_dat$stab_index <- f_stabindex

temp_dat <- rbind(male_dat, female_dat)

gene      <- rep(generation, nrow(temp_dat))
spec     <- rep(species, nrow(temp_dat))
mult     <- rep(z, nrow(temp_dat))
##adds generation to the tempdat
temp_dat$gen <- gene
temp_dat$species <- spec
temp_dat$multiplier <- mult

temp_output <- paste0(par_dir,'Generation_',generation,'/PreSelect',species,generation,'.csv')
###Pre selection file
write.csv(temp_dat, temp_output, row.names=FALSE )

####################subset males by treatment
U1M <- subset(male_dat, male_dat$Tags == "up1")
U2M <- subset(male_dat, male_dat$Tags == "up2")
D1M <- subset(male_dat, male_dat$Tags == "down1")
D2M <- subset(male_dat, male_dat$Tags == "down2")
CM  <- subset(male_dat, male_dat$Tags == "control")
####################subset females by treatment
U1F <- subset(female_dat, female_dat$Tags == "up1")
U2F <- subset(female_dat, female_dat$Tags == "up2")
D1F <- subset(female_dat, female_dat$Tags == "down1")
D2F <- subset(female_dat, female_dat$Tags == "down2")
CF  <- subset(female_dat, female_dat$Tags == "control")

####Path for the select flies output folder
create_csv_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies')
####Creates a folder based off the output path created above
dir.create(create_csv_output)

####Output paths and file names for the different treatmens x sex
U1M_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/U1M.csv')
U2M_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/U2M.csv')
D1M_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/D1M.csv')
D2M_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/D2M.csv')
CM_output  <- paste0(par_dir,'Generation_',generation,'/SelectFlies/CM.csv')
U1F_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/U1F.csv')
U2F_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/U2F.csv')
D1F_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/D1F.csv')
D2F_output <- paste0(par_dir,'Generation_',generation,'/SelectFlies/D2F.csv')
CF_output  <- paste0(par_dir,'Generation_',generation,'/SelectFlies/CF.csv')

CM$Select <- rep(0, nrow(CM))
CMselect <- CM[order(CM[ , 204], decreasing=TRUE), ]
write.csv(CMselect[ , c(193,204,3,5,206)], file = CM_output, row.names = FALSE)
U1M$Select <- rep(0, nrow(U1M))
U1Mselect <- U1M[order(U1M[ , 204], decreasing=TRUE), ]
write.csv(U1Mselect[ , c(193,204,3,5,206)], file = U1M_output, row.names = FALSE)
U2M$Select <- rep(0, nrow(U2M))
U2Mselect <- U2M[order(U2M[ , 205], decreasing=TRUE), ]
write.csv(U2Mselect[ , c(193,205,3,5,206)], file = U2M_output, row.names = FALSE)
D1M$Select <- rep(0, nrow(D1M))
D1Mselect <- D1M[order(D1M[ , 204], decreasing=FALSE), ]
write.csv(D1Mselect[ , c(193,204,3,5,206)], file = D1M_output, row.names = FALSE)
D2M$Select <- rep(0, nrow(D2M))
D2Mselect <- D2M[order(D2M[ , 205], decreasing=FALSE), ]
write.csv(D2Mselect[ , c(193,205,3,5,206)], file = D2M_output, row.names = FALSE)

U2M_stab <- U2Mselect[1:15,]
s_u2m_num <- U2M_stab$flynum
U2M_dir  <- U2M[order(U2M[ , 204], decreasing=TRUE), ][1:15,]
d_u2m_num <- U2M_dir$flynum
different_u2m <- unique(s_u2m_num[! s_u2m_num %in% d_u2m_num])
u2m_diff <- length(different_u2m)

D2M_stab <- D2Mselect[1:15,]
s_d2m_num <- D2M_stab$flynum
D2M_dir  <- D2M[order(D2M[ , 204], decreasing=FALSE), ][1:15,]
d_d2m_num <- D2M_dir$flynum
different_d2m <- unique(s_d2m_num[! s_d2m_num %in% d_d2m_num])
d2m_diff <- length(different_d2m)


CF$Select <- rep(0, nrow(CF))
CFselect <- CF[order(CF[ , 204], decreasing=TRUE), ]
write.csv(CFselect[ , c(193,204,3,5,206)], file = CF_output, row.names = FALSE)
U1F$Select <- rep(0, nrow(U1F))
U1Fselect <- U1F[order(U1F[ , 204], decreasing=TRUE), ]
write.csv(U1Fselect[ , c(193,204,3,5,206)], file = U1F_output, row.names = FALSE)
U2F$Select <- rep(0, nrow(U2F))
U2Fselect <- U2F[order(U2F[ , 205], decreasing=TRUE), ]
write.csv(U2Fselect[ , c(193,205,3,5,206)], file = U2F_output, row.names = FALSE)
D1F$Select <- rep(0, nrow(D1F))
D1Fselect <- D1F[order(D1F[ , 204], decreasing=FALSE), ]
write.csv(D1Fselect[ , c(193,204,3,5,206)], file = D1F_output, row.names = FALSE)
D2F$Select <- rep(0, nrow(D2F))
D2Fselect <- D2F[order(D2F[ , 205], decreasing=FALSE), ]
write.csv(D2Fselect[ , c(193,205,3,5,206)], file = D2F_output, row.names = FALSE)

U2F_stab <- U2Fselect[1:15,]
s_u2f_num <- U2F_stab$flynum
U2F_dir  <- U2F[order(U2F[ , 204], decreasing=TRUE), ][1:15,]
d_u2f_num <- U2F_dir$flynum
different_u2f <- unique(s_u2f_num[! s_u2f_num %in% d_u2f_num])
u2f_diff <- length(different_u2f)

D2F_stab <- D2Fselect[1:15,]
s_d2f_num <- D2F_stab$flynum
D2F_dir  <- D2F[order(D2F[ , 204], decreasing=FALSE), ][1:15,]
d_d2f_num <- D2F_dir$flynum
different_d2f <- unique(s_d2f_num[! s_d2f_num %in% d_d2f_num])
d2f_diff <- length(different_d2f)

stabilizing_report <- cbind('gen' = generation,  'down F' = d2f_diff, 'down M' = d2m_diff, 'up F' = u2f_diff, 'up M' = u2m_diff)
print ('The Number of Flies Not Selected Due to the Stabilizing Selection')
print (stabilizing_report)

stab_rep_output <- paste0(par_dir,'Stabilizing.csv')
file_check <- file_test('-f',stab_rep_output)

if (file_check == FALSE) {
  write.table(stabilizing_report, file = stab_rep_output, sep = ',', row.names = FALSE)
} else {
write.table(stabilizing_report, file = stab_rep_output, sep = ',', row.names = FALSE, col.names = FALSE, append = TRUE)
}


}


select('Mel','16a', 80)


save(select, file = '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/SelectFunction2018StabUpdatePost.rda')
