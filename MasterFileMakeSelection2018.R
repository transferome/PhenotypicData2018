####### Creating a final selection data file.  Has all individual data, plut selection (0,1,2), and includes NAs and Removes
MasterMake <- function (x , y) {
  
  
species <- x
if (species == 'Sim') {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/'
} else {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/'
}

substrRight <- function(x,n) {
  substr(x, nchar(x)-n+1, nchar(x))
}

generation <- y
replicate  <- substrRight(generation,1)

####### Data that includes NA and removed individuals
NA_csv_file_name <- paste0(species,generation,'NAunremoved.csv')
NA_csv_file_path <- paste0(par_dir,'Generation_',generation,'/',NA_csv_file_name)

####### Data that includes selection index etc... but does not include NA and removed individuals
pre_select_file_name <- paste0('PreSelect',species,generation,'.csv')
pre_select_file_path <- paste0(par_dir,'Generation_',generation, '/', pre_select_file_name)

###### Read in these two files
NAdat <- read.table(NA_csv_file_path, sep=',', header=TRUE)

###### Read in the data that was created prior to the selection files
Preselect <- read.table(pre_select_file_path, sep=',', header=TRUE)

###### Get the select files
spath <- paste0(par_dir,'Generation_',generation, '/SelectFlies/')

###### These will be the female files' paths
fcon_path <- paste0(spath,'CF.csv')
fd1_path  <- paste0(spath,'D1F.csv')
fd2_path  <- paste0(spath,'D2F.csv')
fu1_path  <- paste0(spath,'U1F.csv')
fu2_path  <- paste0(spath,'U2F.csv')

###### These will be the male files' paths
mcon_path <- paste0(spath,'CM.csv')
md1_path  <- paste0(spath,'D1M.csv')
md2_path  <- paste0(spath,'D2M.csv')
mu1_path  <- paste0(spath,'U1M.csv')
mu2_path  <- paste0(spath,'U2M.csv')

###### Read in the female files, and substring the sex column to make it F instead of FALSE
fcon <- read.table(fcon_path,sep=',',header=TRUE)
fcon$Sex <- substring(fcon$Sex,1,1)
fd1 <- read.table(fd1_path,sep=',',header=TRUE)
fd1$Sex <- substring(fd1$Sex,1,1)
fd2 <- read.table(fd2_path,sep=',',header=TRUE)
fd2$Sex <- substring(fd2$Sex,1,1)
fu1 <- read.table(fu1_path,sep=',',header=TRUE)
fu1$Sex <- substring(fu1$Sex,1,1)
fu2 <- read.table(fu2_path,sep=',',header=TRUE)
fu2$Sex <- substring(fu2$Sex,1,1)

###### Read in the male files
mcon <- read.table(mcon_path,sep=',',header=TRUE)
md1 <- read.table(md1_path,sep=',',header=TRUE)
md2 <- read.table(md2_path,sep=',',header=TRUE)
mu1 <- read.table(mu1_path,sep=',',header=TRUE)
mu2 <- read.table(mu2_path,sep=',',header=TRUE)


###### Remove second column from all of these dataframes
fcon <- fcon[, -2]
fd1 <- fd1[, -2]
fd2 <- fd2[, -2]
fu1 <- fu1[, -2]
fu2 <- fu2[, -2]
mcon <- mcon[, -2]
md1 <- md1[, -2]
md2 <- md2[, -2]
mu1 <- mu1[, -2]
mu2 <- mu2[, -2]

select.list <- list(fcon,mcon,fd1,fd1,fu1,fu2,md1,md2,mu1,mu2)

####### This code will give warnings if 15 flies were not selected for any of the treatments, or if any zeroes exist within the control select files

for (i in 1:length(select.list)) {
  if (i <= 2) {
    data <- select.list[[i]]
    sele <- data$Select
    v1   <- 1 == sele
    totz <- sum(v1)
       if (totz >= 1) {
         temp_sex <- data$Sex[1]
         warn     <- paste(temp_sex, 'Control Should Not Contain Zeroes', '!!!!!!!', sep = ' ')
         stop(warn)
       }
  
} else {
   data <- select.list[[i]]
   sele <- data$Select
   v1   <- 1 == sele
   s_num <- sum(v1)
       if ( s_num != 15) {
         treat <- levels(data$Tags)
         temp_sex <- data$Sex[1]
         warn  <- paste('You did not select 15 flies from', treat, temp_sex, '!!!!!!!', sep = ' ')
         stop(warn)
       }
}
}


##### rbind these all together to make the selection data
s_dat <- rbind(fcon,fd1,fd2,fu1,fu2,mcon,md1,md2,mu1,mu2)

##### merge this to the preselect dataframe
m1_dat <- merge(Preselect, s_dat, by.x = c('flynum','Tags','Sex'), by.y = c('flynum','Tags','Sex'), all = TRUE)

##### Get the cpr NAs, the leg NAs, and the leg removed flies
legna <- NAdat[is.na(NAdat$removed), ]
temp_na_omit <- na.omit(NAdat)
remove <- temp_na_omit[temp_na_omit$removed == 1, ]
cprna <- NAdat[is.na(NAdat$x1), ]

if (nrow(cprna) == '0') {
  bad_dat <- rbind(legna, remove)
} else {
  bad_dat <- rbind(cprna, legna, remove)
}

#### load plyr package for necessary rbind.fill function
library(plyr)

#### final data
final_dat <- rbind.fill(m1_dat, bad_dat)


final_dat_output <- paste0(par_dir, 'Generation_', generation,'/')
final_dat_name   <- paste0(species,generation,'Master.csv')
output_dat       <- paste0(final_dat_output,final_dat_name)

write.csv(final_dat, file = output_dat, row.names = FALSE)

}

MasterMake('Sim','6b')

save(MasterMake, file = 'Z:/SelectionMelanogaster2018/MakeAMasterFile.rda')
