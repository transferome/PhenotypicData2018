library(plyr)

species <- 'Sim'
if (species == 'Sim') {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionSimulans2018/'
} else {
  par_dir <- '/mnt/topaz/phenotypic_selection/SelectionMelanogaster2018/'
}

gen_h <- 'Generation_'

#### Number of generations into the experiment
gen_num <- 15
#### Vector of post selection generations
postselection_generations <- c("15")

##########################################################
# Combines Data From All Master Files, Removes flies with Wing NAs, Leg NAs, and Bad Outlier Legs
##########################################################
make_data <- function(number_of_generations) {
## Repeats a, b for the number of generations + 1, this is because of Generation 0
  let_rep  <- rep(letters[1:2], number_of_generations + 1)
## Repeats 0 to number_of_generations each twice
  num_rep  <- rep(0:number_of_generations, each = 2)
## Combining let_rep + num_rep resulsts in 0a, 0b, 1a, 1b, etc...
  gen_name <- paste0(num_rep, let_rep)
## Adds Generation_ to each element of gen_name
  gen_fol  <- paste0(gen_h, gen_name,'/')
## Adds the parent directory path to the Generation_# folder
  sel_fol  <- paste0(par_dir,gen_fol)
## The Master file name will be the species + the gen_name + Master.csv
  mast_name <- paste0(species, gen_name, 'Master.csv')
## This will result in a vector of the filepath and filename for each Master.csv
  mast_path <- paste0(sel_fol,mast_name)
  print('Check that Master File Paths Look Okay')
  print(mast_path)
##########################################################

##########################################################
########## Read in Data for Each Master File
########## Remove NA Data and Badleg Data
##########################################################
## This makes an empty list that will hold the Master.csv data
  data_list <- list()
## This loop will read in each Master.csv file in the mast_path vector and save it to the data list
  for (i in 1:length(mast_path)) {
    data_list[[i]] <- read.table(mast_path[i], sep = ',', header=TRUE)
  }
## This will rbind all the dataframes in data_list, but fill in missing info with NAs
## This is because Gen0 includes some variables that the subsequent generations don't
  dat <- rbind.fill(data_list)
## Check that Data is in the usual order
  column_check1 <- which(colnames(dat) == "date") #209
  column_check2 <- which(colnames(dat) == "time") #210
  if (column_check1 == 209 & column_check2 == 210) {
    print('Rbind Looks Normal. Data Merged Correctly.')
  } else {
    print ('Data Looks to Be Out of Order')
  }
## Remove the extra date and time columns
  dat <- dat[ , -c(209,210)]
## Remove any individual that has an NA in the removed column.These are flies whose legs did not spline
  dat <- dat[!is.na(dat$removed), ]
## Remove any individual that has a value of 1 in the removed column.These are flies whose legs splined but were incorrect. BadLandmark/BadImage/BadSpline
  dat <- dat[dat$removed != 1, ]
## Remove any individual that has an NA in the x1 columns.These are individuals whose wings did not spline.
  dat <- dat[!is.na(dat$x1), ]

## this will create a vector of what replicate an individual is in
  rep_vec <- gsub('[[:digit:]]+','', dat$gen)
## this will create a vector of what generation number an individual is in
  gen_vec <- gsub('[[:alpha:]]+','',dat$gen)

  dat$rep        <- rep_vec
  dat$gen_number <- gen_vec
  return(dat)
}
dat <- make_data(gen_num)
###########################################################
###########################################################

#Currently out of place wing data creation
# 
wng_outliers <- boxplot(dat$wngdev)$out
wngoutrow    <- which(dat$wngdev %in% wng_outliers)
wngout <- dat[wngoutrow, ]

chk_wing <- wngout[order(wngout$wngdev, decreasing = TRUE), ]
write.csv(chk_wing, file = 'z:/SelectionMelanogaster2018/CheckWings.csv', row.names = FALSE)
which.max(wngout$wngdev)
# # 
boxplot(dat$dir_index)$out
# # 
dat[which.max(dat$wngdev),]
dat[which.max(dat$dir_index),]

###########################################################
gen0 <- dat[dat$gen_number == 0, ]
# standard deviation of population at Generation 0 for replicate a
sd0a <- with(gen0[gen0$rep == 'a', ], sd(dir_index))
# standard deviation of population at Generation 0 for replicate b
sd0b <- with(gen0[gen0$rep == 'b', ], sd(dir_index))



###########################################################
############ Generation 0 means, variances, N for selected and unselected flies
###########################################################
gen_zero_means <- function(gen0_data) {
## Directional_Index mean for Generation 0, mean of controls only, and mean of all other flies
  gen0_select_means  <- with(gen0_data[gen0_data$treatment != 'control', ], aggregate(dir_index, by =list(gen), mean))
  select_mean_0a <- gen0_select_means[1,2]
  select_mean_0b <- gen0_select_means[2,2]
  gen0_control_means <- with(gen0_data[gen0_data$treatment == 'control', ], aggregate(dir_index, by =list(gen), mean))
  control_mean_0a <- gen0_control_means[1,2]
  control_mean_0b <- gen0_control_means[2,2]

## Directional_Index of the selected flies
  g0smeans <- with(gen0_data[gen0_data$Select == 1, ], aggregate(dir_index, by = list(gen, treatment), mean))
  colnames(g0smeans) <- c('generation','treatment','Smean')
  g0M      <- c(control_mean_0a, control_mean_0b, rep(c(select_mean_0a, select_mean_0b), times = 4))
  g0smeans$Mean <- g0M
  g0s <- g0smeans$Smean - g0smeans$Mean
  g0smeans$S <- g0s
  return(g0smeans)
}
gen0means <- gen_zero_means(gen0)
gen_zero_varN <- function(gen0_data) {
## Directional_Index variance for Generation 0, variance of controls only, and variance of all other flies
  gen0_select_var  <- with(gen0_data[gen0_data$treatment != 'control', ], aggregate(dir_index, by =list(gen), var))
  select_var_0a <- gen0_select_var[1,2]
  select_var_0b <- gen0_select_var[2,2]
  gen0_control_var <- with(gen0_data[gen0_data$treatment == 'control', ], aggregate(dir_index, by =list(gen), var))
  control_var_0a <- gen0_control_var[1,2]
  control_var_0b <- gen0_control_var[2,2]

## Directional_Index of the selected flies
  g0svar <- with(gen0_data[gen0_data$Select == 1, ], aggregate(dir_index, by = list(gen, treatment), var))
  colnames(g0svar) <- c('generation','treatment','SVp')
  g0v      <- c(control_var_0a, control_var_0b, rep(c(select_var_0a, select_var_0b), times = 4))
  g0svar$Vp <- g0v

## Directional_Index sample size for Generation 0, mean of controls only, and sample size of all other flies
  gen0_select_size  <- with(gen0_data[gen0_data$treatment != 'control', ], aggregate(dir_index, by =list(gen), length))
  select_size_0a <- gen0_select_size[1,2]
  select_size_0b <- gen0_select_size[2,2]
  gen0_control_size <- with(gen0_data[gen0_data$treatment == 'control', ], aggregate(dir_index, by =list(gen), length))
  control_size_0a <- gen0_control_size[1,2]
  control_size_0b <- gen0_control_size[2,2]

  g0size <- c(control_size_0a, control_size_0b, rep(c(select_size_0a, select_size_0b), times = 4))
  g0svar$N <- g0size

  gen0_selected <- with(gen0_data[ gen0_data$Select == 1, ], aggregate(dir_index, by =list(gen,treatment), length))
  g0svar$Ns <- gen0_selected$x
  return(g0svar)
}
gen0varN <- gen_zero_varN(gen0)

##########################################################
# Generation 1+ means, variances, N for selected and unselected
##########################################################
selection_generations <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14")
gen_selection <- dat[dat$gen_number %in% selection_generations, ]

gen_sel_means <- function(gen_plus){
##########################################################################
# Mean of index selection generations
##########################################################################
  gen_s_means <- with(gen_plus[gen_plus$Select == 1, ], aggregate(dir_index, by = list(gen, Tags), mean))
  colnames(gen_s_means) <- c('generation','treatment','Smean')
  gen_pop_means <- with(gen_plus, aggregate(dir_index, by = list(gen,Tags),mean))
  gen_s_means$Mean <- gen_pop_means$x
  gen_s_s <- gen_s_means$Smean - gen_s_means$Mean
  gen_s_means$S <- gen_s_s
  return(gen_s_means)
}
genSmeans <- gen_sel_means(gen_selection)
# function below does variance and sample size N
gen_sel_varN <- function(gen_plus){
  gen_s_var <- with(gen_plus[gen_plus$Select == 1, ], aggregate(dir_index, by = list(gen, Tags), var))
  colnames(gen_s_var) <- c('generation','treatment','SVp')
  gen_pop_var <- with(gen_plus, aggregate(dir_index, by = list(gen,Tags),var))
  gen_s_var$Vp <- gen_pop_var$x
  # population size/sample size
  genp_size <- with(gen_plus, aggregate(dir_index, by = list(gen,Tags),length))
  gen_s_var$N <- genp_size$x
  genp_selected <- with(gen_plus[gen_plus$Select == 1, ], aggregate(dir_index, by = list(gen, Tags), length))
  gen_s_var$Ns <- genp_selected$x
  return(gen_s_var)
}
genSvarN <- gen_sel_varN(gen_selection)

##########################################################
# Generation 15+ post selection means, variances, N
##########################################################
gen_postsel <- dat[dat$gen_number %in% postselection_generations, ]
gen_post_means <- function(gen_post) {
# populations means, assuming Mean of selected is same as mean of population so S = 0
  gen_post_means <- with(gen_post, aggregate(dir_index, by = list(gen, Tags), mean))
  colnames(gen_post_means) <- c('generation','treatment','Smean')
  gen_post_means$Mean <- gen_post_means$Smean
  gen_post_s <- gen_post_means$Smean - gen_post_means$Mean
  gen_post_means$S <- gen_post_s
  return(gen_post_means)
}
genPmeans <- gen_post_means(gen_postsel)
gen_post_varN <- function(gen_post) {
# calculates variance of unselected and selected individuals, in this case variance selected = variance unselected
  gen_post_var <- with(gen_post, aggregate(dir_index, by = list(gen, Tags), var))
  colnames(gen_post_var) <- c('generation','treatment','SVp')
  gen_post_var$Vp <- gen_post_var$SVp
# sample size for post selection generations
  genpost_size <- with(gen_post, aggregate(dir_index, by = list(gen,Tags),length))
  gen_post_var$N <- genpost_size$x
  genpost_selected <- rep(15, length(gen_post_var$N))
  gen_post_var$Ns <- genpost_selected
  return(gen_post_var)
}
genPvarN <- gen_post_varN(gen_postsel)

##########################################################
# Combine Generation 0 and 1+ Means, Var, and N, Organize and Clean data
##########################################################
organize_dat <- function(g0smeans, gen_s_means, gen_post_means, g0svar, gen_s_var, gen_post_var) {
  sel_dat <- rbind(g0smeans,gen_s_means, gen_post_means)
  sel_dat$rep    <- gsub('[[:digit:]]+','', sel_dat$generation)
  sel_dat$gennum <- gsub('[[:alpha:]]+','', sel_dat$generation)
  # Order data by generation, number of generations, and then reset the rownames
  sel_dat <- sel_dat[order(sel_dat$generation), ]
  rownames(sel_dat) <- NULL
  # Add variance and population size combine variance data first
  vdat <- rbind(g0svar,gen_s_var, gen_post_var)
  vdat <- vdat[order(vdat$generation), ]
  rownames(vdat) <- NULL
  sel_dat$SVp <- vdat$SVp
  sel_dat$Vp  <- vdat$Vp
  sel_dat$N   <- vdat$N
  sel_dat$Ns  <- vdat$Ns
  return(sel_dat)
}

s_dat <- organize_dat(gen0means, genSmeans, genPmeans, gen0varN, genSvarN, genPvarN)
head(s_dat)
tail(s_dat)

##########################################################
# Needed for Simulans
##########################################################
#rbind(existingDF[1:r,],newrow,existingDF[-(1:r),])
if(species == 'Sim') {
  newrow <- list('11a','down1',0.009701307,0.009701307,0,'a','11', 0, 8.929542e-05, 65, 0 )
  names(newrow) <- c('generation','treatment','Smean','Mean','S','rep','gennum','SVp','Vp','N','Ns')
  s_dat <- rbind(s_dat[1:111,],newrow,s_dat[-(1:111),])
  rownames(s_dat) <- NULL
}


######################################################
######Calculating the Cumulative S, Response, and SE of Response
######################################################
create_response_data <- function(sdat, sd_0a, sd_0b){
##### element vectors for each treatment of each replicate
##### Each number in these vectors, will be the row number for a treatment and replicate
##### Example d1avec will list the row numbers for down1 replicate a
  # control a vec
  t1 <- seq(1 , nrow(sdat), 10)
  # d1avec
  t2 <- seq(2 , nrow(sdat), 10)
  # d2avec
  t3 <- seq(3 , nrow(sdat), 10)
  # u1avec
  t4 <- seq(4 , nrow(sdat), 10)
  # u2avec
  t5 <- seq(5 , nrow(sdat), 10)
  # control b vec
  t6 <- seq(6 , nrow(sdat), 10)
  # d1bvec
  t7 <- seq(7 , nrow(sdat), 10)
  # d2bvec
  t8 <- seq(8 , nrow(sdat), 10)
  # u1bvec
  t9 <- seq(9 , nrow(sdat), 10)
  # u2bvec
  t10 <- seq(10 , nrow(sdat), 10)

# These are dummy vectors that will hold the S value for each treatment/rep for every generation
  conaS <- rep(NA, length(t1))
  d1aS  <- rep(NA, length(t1))
  d2aS  <- rep(NA, length(t1))
  u1aS  <- rep(NA, length(t1))
  u2aS  <- rep(NA, length(t1))
  conbS <- rep(NA, length(t1))
  d1bS  <- rep(NA, length(t1))
  d2bS  <- rep(NA, length(t1))
  u1bS  <- rep(NA, length(t1))
  u2bS  <- rep(NA, length(t1))

# vector of S values for each treatment put into the dummy vectors above
  for (i in 1:length(t1)) {
    conaS[i] <- sdat[t1[i],5]
    d1aS[i]  <- sdat[t2[i],5]
    d2aS[i]  <- sdat[t3[i],5]
    u1aS[i]  <- sdat[t4[i],5]
    u2aS[i]  <- sdat[t5[i],5]
    conbS[i] <- sdat[t6[i],5]
    d1bS[i]  <- sdat[t7[i],5]
    d2bS[i]  <- sdat[t8[i],5]
    u1bS[i]  <- sdat[t9[i],5]
    u2bS[i]  <- sdat[t10[i],5]
  }

# cumsum function will find the cumulative S from the S vectors created above
# The final element in these vectors below will be the total cumulative S
# the cumulative selection differential for each treatment
  conaCS <- cumsum(conaS)
  d1aCS <- cumsum(d1aS)
  d2aCS <- cumsum(d2aS)
  u1aCS <- cumsum(u1aS)
  u2aCS <- cumsum(u2aS)
  conbCS <- cumsum(conbS)
  d1bCS <- cumsum(d1bS)
  d2bCS <- cumsum(d1bS)
  u1bCS <- cumsum(u1bS)
  u2bCS <- cumsum(u2bS)

# Take out the initial population mean value at generation 0
# Response at Generation 0, the starting value
  conaM <- sdat[t1[1], 4]
  d1aM  <- sdat[t2[1], 4]
  d2aM  <- sdat[t3[1], 4]
  u1aM  <- sdat[t4[1], 4]
  u2aM  <- sdat[t5[1], 4]
  conbM <- sdat[t6[1], 4]
  d1bM  <- sdat[t7[1], 4]
  d2bM  <- sdat[t8[1], 4]
  u1bM  <- sdat[t9[1], 4]
  u2bM  <- sdat[t10[1], 4]

# Dummy vectors will hold the value of the population mean at gen# - population mean at gen0
  conaR <- rep(NA, length(t1))
  d1aR <- rep(NA, length(t1))
  d2aR <- rep(NA, length(t1))
  u1aR <- rep(NA, length(t1))
  u2aR <- rep(NA, length(t1))
  conbR <- rep(NA, length(t1))
  d1bR <- rep(NA, length(t1))
  d2bR <- rep(NA, length(t1))
  u1bR <- rep(NA, length(t1))
  u2bR <- rep(NA, length(t1))
# Vectors of the response which is the mean of each generation substracted from the starting mean
  for (i in 1:length(t1)) {
    conaR[i] <- sdat[t1[i], 4] - conaM
    d1aR[i]  <- sdat[t2[i], 4] - d1aM
    d2aR[i]  <- sdat[t3[i], 4] - d2aM
    u1aR[i]  <- sdat[t4[i], 4] - u1aM
    u2aR[i]  <- sdat[t5[i], 4] - u2aM
    conbR[i] <- sdat[t6[i], 4] - conbM
    d1bR[i]  <- sdat[t7[i], 4] - d1bM
    d2bR[i]  <- sdat[t8[i], 4] - d2bM
    u1bR[i]  <- sdat[t9[i], 4] - u1bM
    u2bR[i]  <- sdat[t10[i], 4] - u2bM
  }

# Need to add a zero to the cumulative selection differential because the selection differential at generation 0 is 0
  conaCSt <- c(0, conaCS)
  d1aCSt <- c(0,d1aCS)
  d2aCSt <- c(0,d2aCS)
  u1aCSt <- c(0,u1aCS)
  u2aCSt <- c(0,u2aCS)
  conbCSt <- c(0, conbCS)
  d1bCSt <- c(0,d1bCS)
  d2bCSt <- c(0,d2bCS)
  u1bCSt <- c(0,u1bCS)
  u2bCSt <- c(0,u2bCS)

# What is the length of the Cumulative selection differential vectors
  len_temp <- length(d1aCSt)
# Not interested in the final value in the selection differential vector because that is for the next generation for which we have no data yet
  conaSD <- conaCSt[-len_temp]
  d1aSD <- d1aCSt[-len_temp] 
  d2aSD <- d2aCSt[-len_temp] 
  u1aSD <- u1aCSt[-len_temp] 
  u2aSD <- u2aCSt[-len_temp] 
  conbSD <- conbCSt[-len_temp]
  d1bSD <- d1bCSt[-len_temp] 
  d2bSD <- d2bCSt[-len_temp] 
  u1bSD <- u1bCSt[-len_temp] 
  u2bSD <- u2bCSt[-len_temp] 

# Now create a vector that will be the cumulative S for each treatment/rep through the generations
# dummy vector for cumulative S
  CumS <- rep(NA, nrow(sdat))
# Fills in the vector
  for (i in 1:length(t1)) {
    CumS[t1[i]] <- conaSD[i]
    CumS[t2[i]] <- d1aSD[i]
    CumS[t3[i]] <- d2aSD[i]
    CumS[t4[i]] <- u1aSD[i]
    CumS[t5[i]] <- u2aSD[i]
    CumS[t6[i]] <- conbSD[i]
    CumS[t7[i]] <- d1bSD[i]
    CumS[t8[i]] <- d2bSD[i]
    CumS[t9[i]] <- u1bSD[i]
    CumS[t10[i]] <- u2bSD[i]
  }

# Now create a vector that will be the response for each treatment/rep through the generations
  Resp <- rep(NA, nrow(sdat))

  for (i in 1:length(t1)) {
    Resp[t1[i]]  <- conaR[i]
    Resp[t2[i]]  <- d1aR[i]
    Resp[t3[i]]  <- d2aR[i]
    Resp[t4[i]]  <- u1aR[i]
    Resp[t5[i]]  <- u2aR[i]
    Resp[t6[i]]  <- conbR[i]
    Resp[t7[i]]  <- d1bR[i]
    Resp[t8[i]]  <- d2bR[i]
    Resp[t9[i]]  <- u1bR[i]
    Resp[t10[i]] <- u2bR[i]
  }

# Add these vectors to the dataframe 
  sdat$cumS <- CumS
  sdat$resp <- Resp

# Response divided by the standard deviation of the population at generation 0
  resp_sd_temp <- rep(NA, length(sdat$resp))

  for (i in 1:length(sdat$resp)) {
    val     <- sdat$resp[i]
    rep_val <- sdat$rep[i]
    if (rep_val == 'a') {
    resp_sd_temp[i] <- val/sd_0a
    } else {
    resp_sd_temp[i] <- val/sd_0b
    }
  }
# Add vectors to dataframe that are the reponse in units of SD and the selection intensity
  sdat$respsd <- resp_sd_temp
  sdat$intensity <- sdat$S/sqrt(sdat$Vp)

# Get standard error between the difference of the response generation and generation 0
  se.ca <- rep(0, length(t1))
  se.d1a <- rep(0, length(t2))
  se.d2a <- rep(0, length(t3))
  se.u1a <- rep(0, length(t4))
  se.u2a <- rep(0, length(t5))
  se.cb <- rep(0, length(t6))
  se.d1b <- rep(0, length(t7))
  se.d2b <- rep(0, length(t8))
  se.u1b <- rep(0, length(t9))
  se.u2b <- rep(0, length(t10))

# 
  for (i in 1:(length(t1) - 1 )) {
  se.ca[i+1] <- sqrt((sdat[t1[i+1],9]/sdat[t1[i+1],10]) + (sdat[t1[i],9]/sdat[t1[i],10]))
  se.d1a[i+1] <- sqrt((sdat[t2[i+1],9]/sdat[t2[i+1],10]) + (sdat[t2[i],9]/sdat[t2[i],10]))
  se.d2a[i+1] <- sqrt((sdat[t3[i+1],9]/sdat[t3[i+1],10]) + (sdat[t3[i],9]/sdat[t3[i],10]))
  se.u1a[i+1] <- sqrt((sdat[t4[i+1],9]/sdat[t4[i+1],10]) + (sdat[t4[i],9]/sdat[t4[i],10]))
  se.u2a[i+1] <- sqrt((sdat[t5[i+1],9]/sdat[t5[i+1],10]) + (sdat[t5[i],9]/sdat[t5[i],10]))
  se.cb[i+1] <- sqrt((sdat[t6[i+1],9]/sdat[t6[i+1],10]) + (sdat[t6[i],9]/sdat[t6[i],10]))
  se.d1b[i+1] <- sqrt((sdat[t7[i+1],9]/sdat[t7[i+1],10]) + (sdat[t7[i],9]/sdat[t7[i],10]))
  se.d2b[i+1] <- sqrt((sdat[t8[i+1],9]/sdat[t8[i+1],10]) + (sdat[t8[i],9]/sdat[t8[i],10]))
  se.u1b[i+1] <- sqrt((sdat[t9[i+1],9]/sdat[t9[i+1],10]) + (sdat[t9[i],9]/sdat[t9[i],10]))
  se.u2b[i+1] <- sqrt((sdat[t10[i+1],9]/sdat[t10[i+1],10]) + (sdat[t10[i],9]/sdat[t10[i],10]))
  }

# Now create a vector that will be the response for each treatment/rep through the generations
  R.se <- rep(NA, nrow(sdat))

  for (i in 1:length(t1)) {
    R.se[t1[i]]  <- se.ca[i]
    R.se[t2[i]]  <- se.d1a[i]
    R.se[t3[i]]  <- se.d2a[i]
    R.se[t4[i]]  <- se.u1a[i]
    R.se[t5[i]]  <- se.u2a[i]
    R.se[t6[i]]  <- se.cb[i]
    R.se[t7[i]]  <- se.d1b[i]
    R.se[t8[i]]  <- se.d2b[i]
    R.se[t9[i]]  <- se.u1b[i]
    R.se[t10[i]] <- se.u2b[i]
  }

  zero.errordat <- sdat[1:10, 9:10]
  zero.error <- sqrt(zero.errordat[,1])/sqrt(zero.errordat[,2])
  R.se[1:10] <- zero.error
  sdat$R.error <- R.se

#### response divided by the standard deviation of the population at generation 0
  resp_error_temp <- rep(NA, length(sdat$R.error))

  for (i in 1:length(sdat$R.error)) {
    val     <- sdat$R.error[i]
    rep_val <- sdat$rep[i]
    if (rep_val == 'a') {
      resp_error_temp[i] <- val/sd_0a
    } else {
      resp_error_temp[i] <- val/sd_0b
    }
  }

  sdat$R.error.sd <- resp_error_temp *1.96
  return(sdat)
}
s_df <- create_response_data(s_dat, sd0a, sd0b)
head(s_df)
tail(s_df)
#############################################################
################## Calculate Heritabilities #################
#############################################################
## the final generation information will be the max value of the generation numbers vector (gennum)
#final_gen <- max(sdat$gennum)
## final gen is now just 15
final_gen <- 15
## the response and cumulative S data here is all that is needed to calculate the heritability
hdat <- s_df[s_df$gennum == final_gen, ]
## remove the controls from hdat
hdat <- hdat[hdat$treatment != 'control', ]
## response/cumulative S is = h^2
herit <- hdat$resp/hdat$cumS
herit <- as.data.frame(herit)
## Add row names
herit$treatment <- c('down1 A', 'down2 A', 'up1 A', 'up2 A', 'down1 B', 'down2 B', 'up1 B', 'up2 B')
herit

###########################################################
# Standard Deviation at Generation zero, defining Gen zero data

#### TODO Future Wing Stuff
# standard deviation at Generation 0 for replicate a and b, for wngdev
sd_wng_0a <- with(gen0[gen0$rep == 'a', ], sd(wngdev))
sd_wng_0b <- with(gen0[gen0$rep == 'b', ], sd(wngdev))
wngdat <- with(dat, aggregate(wngdev, by =list(gen,Tags), mean))
wngdat_var <- with(dat, aggregate(wngdev, by =list(gen,Tags), var))
#remove generation 0 for now
# Mel A has 0.0003424 at Generation 0 with var of 3.572107x10-8
# Mel B has 0.0003340 at Generation 0 with var of 3.444124x10-8
wngdat <- wngdat[-c(1:2), ]
wngdat_var <- wngdat_var[-c(1:2), ]
names(wngdat) <- c('generation', 'treatment', 'wngdev')
names(wngdat_var) <- c('generation', 'treatment', 'var')
wngdat$rep <- gsub('[[:digit:]]+','', wngdat$generation)
wngdat$gen <- gsub('[[:alpha:]]+','', wngdat$generation)
wngdat_var$rep <- gsub('[[:digit:]]+','', wngdat_var$generation)
wngdat_var$gen <- gsub('[[:alpha:]]+','', wngdat_var$generation)
wngdat <- wngdat[order(wngdat$generation), ]
rownames(wngdat) <- NULL
wngdat_var <- wngdat_var[order(wngdat_var$generation), ]
rownames(wngdat_var) <- NULL
if(species == 'Sim') {
  newrow <- list('11a','down1',0.0009774728, 'a','11' )
  names(newrow) <- c('generation','treatment','wngdev','rep', 'gen')
  wngdat <- rbind(wngdat[1:101,],newrow,wngdat[-(1:101),])
  rownames(wngdat) <- NULL
}
if(species == 'Sim') {
  newrow <- list('11a','down1',0.00000002553220, 'a','11' )
  names(newrow) <- c('generation','treatment','var','rep', 'gen')
  wngdat_var <- rbind(wngdat_var[1:101,],newrow,wngdat_var[-(1:101),])
  rownames(wngdat_var) <- NULL
}

# Response divided by the standard deviation of the population at generation 0
wresp_sd_temp <- rep(NA, length(wngdat$wngdev))
for (i in 1:length(wngdat$wngdev)) {
  val     <- wngdat$wngdev[i]
  rep_val <- wngdat$rep[i]
  if (rep_val == 'a') {
    wresp_sd_temp[i] <- val/sd_wng_0a
  } else {
    wresp_sd_temp[i] <- val/sd_wng_0b
  }
}
# Add vectors to dataframe that are the reponse in units of SD and the selection intensity
wngdat$wngdevsd <- wresp_sd_temp
wngdat$N <- s_df$N[11:length(s_df$N)]
wngdat$var <- wngdat_var$var

##### Subset data into A and B replicate, and then by the different treatments
datA <- s_df[s_df$rep == 'a', ]
controlA <- datA[datA$treatment == 'control',]
d1A  <- datA[datA$treatment == 'down1', ]
d2A  <- datA[datA$treatment == 'down2', ]
u1A  <- datA[datA$treatment == 'up1', ]
u2A  <- datA[datA$treatment == 'up2', ]
datB <- s_df[s_df$rep == 'b', ]
controlB <- datB[datB$treatment == 'control',]
d1B  <- datB[datB$treatment == 'down1', ]
d2B  <- datB[datB$treatment == 'down2', ]
u1B  <- datB[datB$treatment == 'up1', ]
u2B  <- datB[datB$treatment == 'up2', ]


if (species == 'Sim') {
  row_remove <- which(grepl('11a',d1A$generation))
  d1A <- d1A[-row_remove, ]
}

if (species == 'Sim') {
  title <- 'D. simulans'
} else {
  title <- 'D. melanogaster'
}

a_resp_range1 <- range(datA$respsd + datA$R.error.sd)
a_resp_range2 <- range(datA$respsd - datA$R.error.sd)
a_resp_range <- range(a_resp_range1, a_resp_range2)
b_resp_range1 <- range(datB$respsd + datB$R.error.sd)
b_resp_range2 <- range(datB$respsd - datB$R.error.sd)
b_resp_range  <- range(b_resp_range1, b_resp_range2)
rng <- c(a_resp_range, b_resp_range)
rng_gen <- range(as.numeric(s_df$gennum))


#pdf('z:/SelectionSimulans2018/MelResponseGeneration.pdf')

### This makes a plot of the response by the generation
plot(1, xlim = c(min(rng_gen), max(rng_gen)), ylim = c(min(rng), max(rng)),
     type = 'n', xlab = 'Generation', 
  ylab = 'Response in Phenotypic SD', las = 1, main = substitute(italic(x), list(x=title) ))

points(u1A$gennum, u1A$respsd, pch = 15, col ='purple3')
arrows(as.numeric(u1A$gennum), u1A$respsd - u1A$R.error.sd, as.numeric(u1A$gennum), 
       u1A$respsd + u1A$R.error.sd, length = 0.03, angle = 90, col = 'purple3',code = 3 )

points(u2A$gennum, u2A$respsd, pch = 15, col = 'purple3')
arrows(as.numeric(u2A$gennum), u2A$respsd - u2A$R.error.sd, as.numeric(u2A$gennum), 
       u2A$respsd + u2A$R.error.sd, length = 0.03, angle = 90, col = 'purple3', lty = 6,code = 3 )

points(u1B$gennum, u1B$respsd, pch = 16, col = 'gray30')
arrows(as.numeric(u1B$gennum), u1B$respsd - u1B$R.error.sd, as.numeric(u1B$gennum), 
       u1B$respsd + u1B$R.error.sd, length = 0.03, angle = 90, col = 'gray30',code = 3 )

points(u2B$gennum, u2B$respsd, pch = 16, col = 'gray30')
arrows(as.numeric(u2B$gennum), u2B$respsd - u2B$R.error.sd, as.numeric(u2B$gennum), 
       u2B$respsd + u2B$R.error.sd, length = 0.03, angle = 90, col = 'gray30', lty = 6, code = 3 )



lines(u1A$gennum, u1A$respsd, lty = 1, lwd = 2, col ='purple3')
lines(u2A$gennum, u2A$respsd, lty = 6, lwd = 2, col = 'purple3')

lines(u1B$gennum, u1B$respsd, lty = 1, lwd = 2, col = 'gray30')
lines(u2B$gennum, u2B$respsd, lty = 6, lwd = 2, col = 'gray18')


points(d1A$gennum, d1A$respsd, pch = 15, col ='purple3')
arrows(as.numeric(d1A$gennum), d1A$respsd - d1A$R.error.sd, as.numeric(d1A$gennum), 
       d1A$respsd + d1A$R.error.sd, length = 0.03, angle = 90, col = 'purple3',code = 3 )

points(d2A$gennum, d2A$respsd, pch = 15, col = 'purple3')
arrows(as.numeric(d2A$gennum), d2A$respsd - d2A$R.error.sd, as.numeric(d2A$gennum), 
       d2A$respsd + d2A$R.error.sd, length = 0.03, angle = 90, col = 'purple3',lty = 6,code = 3 )

points(d1B$gennum, d1B$respsd, pch = 16, col = 'gray30')
arrows(as.numeric(d1B$gennum), d1B$respsd - d1B$R.error.sd, as.numeric(d1B$gennum), 
       d1B$respsd + d1B$R.error.sd, length = 0.03, angle = 90, col = 'gray30',code = 3 )

points(d2B$gennum, d2B$respsd, pch = 16, col = 'gray30')
arrows(as.numeric(d2B$gennum), d2B$respsd - d2B$R.error.sd, as.numeric(d2B$gennum), 
       d2B$respsd + d2B$R.error.sd, length = 0.03, angle = 90, col = 'gray30', lty = 6, code = 3 )

lines(d1A$gennum, d1A$respsd, lty = 1, lwd = 2,  col ='purple3')
lines(d2A$gennum, d2A$respsd, lty = 6, lwd = 2, col = 'purple3')
lines(d1B$gennum, d1B$respsd, lty = 1, lwd = 2, col = 'gray30')
lines(d2B$gennum, d2B$respsd, lty = 6, lwd = 2, col = 'gray30')

### control or not (ctrl+shift+c)
points(controlA$gennum, controlA$respsd, pch = 18, col ='purple3')
arrows(as.numeric(controlA$gennum), controlA$respsd - controlA$R.error.sd, as.numeric(controlA$gennum),
       controlA$respsd + controlA$R.error.sd, length = 0.03, angle = 90, col = 'purple3',code = 3 )
lines(controlA$gennum, controlA$respsd, lty = 3, lwd = 2, col ='purple3')
points(controlB$gennum, controlB$respsd, pch = 18, col ='gray30')
arrows(as.numeric(controlB$gennum), controlB$respsd - controlB$R.error.sd, as.numeric(controlB$gennum),
       controlB$respsd + controlB$R.error.sd, length = 0.03, angle = 90, col = 'gray30',code = 3 )
lines(controlB$gennum, controlB$respsd, lty = 3, lwd = 2, col ='gray30')

### stabilizing increase indicators
# if (species == 'Sim') {
#   points(3, -2.5, pch = 8, col = 'black')
#   points(12, -5.2, pch = 8, col = 'black')
# } else {
#   points(3, -3, pch = 8, col = 'black')
#   points(12, -6, pch = 8, col = 'black')
# }


legend('topleft', legend = c('Directional A', 'Direct + Stab A', 'Directional B', 'Direct + Stab B', 'Control A', 'Control B'), 
       col = c('purple3','purple3','gray30','gray30', 'purple3', 'gray30'), lty = c(1,6,1,6,3,3), lwd = c(2,2,2,2,2,2), cex = 0.78)

#dev.off()

### find the range of the selection intensity
a_int_range <- range(datA$intensity)
b_int_range <- range(datB$intensity)
rngi <- c(a_int_range, b_int_range)

a_var_range <- range(datA$Vp)
b_var_range <- range(datB$Vp)
rngv <- c(a_var_range, b_var_range)

### This makes a plot of the response by the generation
plot(1, xlim = c(min(rng_gen), max(rng_gen)), ylim = c(min(rngi), max(rngi)), type = 'n', xlab = 'Generation', 
     ylab = 'Selection Intensity', las = 1)
points(u1A$gennum, u1A$intensity, pch = 15, col ='purple3')
points(u2A$gennum, u2A$intensity, pch = 15, col = 'purple3')
points(u1B$gennum, u1B$intensity, pch = 16, col = 'gray30')
points(u2B$gennum, u2B$intensity, pch = 16, col = 'gray30')
lines(u1A$gennum, u1A$intensity, lty = 1, lwd = 2, col ='purple3')
lines(u2A$gennum, u2A$intensity, lty = 4, col = 'purple3')
lines(u1B$gennum, u1B$intensity, lty = 1, lwd = 2, col = 'gray30')
lines(u2B$gennum, u2B$intensity, lty = 4, col = 'gray18')
points(d1A$gennum, d1A$intensity, pch = 15, col ='purple3')
points(d2A$gennum, d2A$intensity, pch = 15, col = 'purple3')
points(d1B$gennum, d1B$intensity, pch = 16, col = 'gray30')
points(d2B$gennum, d2B$intensity, pch = 16, col = 'gray30')
lines(d1A$gennum, d1A$intensity, lty = 1, lwd = 2,  col ='purple3')
lines(d2A$gennum, d2A$intensity, lty = 4, col = 'purple3')
lines(d1B$gennum, d1B$intensity, lty = 1, lwd = 2, col = 'gray30')
lines(d2B$gennum, d2B$intensity, lty = 4, col = 'gray30')
legend('left', legend = c('Directional A', 'Direct + Stab A', 'Directional B', 'Direct + Stab B'), col = c('purple3','purple3','gray30','gray30'), lty = c(1,4,1,4), lwd = c(2,1,2,1), cex = 0.75)



s_dfsel <- s_df[s_df$gennum %in% c("0", selection_generations), ]

##### Subset data into A and B replicate, and then by the different treatments
datA <- s_dfsel[s_dfsel$rep == 'a', ]
d1A  <- datA[datA$treatment == 'down1', ]
d2A  <- datA[datA$treatment == 'down2', ]
u1A  <- datA[datA$treatment == 'up1', ]
u2A  <- datA[datA$treatment == 'up2', ]
datB <- s_dfsel[s_dfsel$rep == 'b', ]
d1B  <- datB[datB$treatment == 'down1', ]
d2B  <- datB[datB$treatment == 'down2', ]
u1B  <- datB[datB$treatment == 'up1', ]
u2B  <- datB[datB$treatment == 'up2', ]

#### Range of SD response
a_resp_range  <- range(datA$respsd)
b_resp_range  <- range(datB$respsd)
rng <- c(a_resp_range, b_resp_range)
#### Range of the cumulative S
a_cums_range <- range(datA$cumS)
b_cums_range <- range(datB$cumS)
rng_cums  <- c(a_cums_range, b_cums_range)



#### Graphs the response by the cumulative selection differential
plot(1, xlim = c(min(rng_cums), max(rng_cums)), ylim = c(min(rng), max(rng)), type = 'n', 
  xlab = 'Cumulative Selection Differential', ylab = 'Response in Phenotypic SD', las = 1, main = substitute(italic(x), list(x=title) ) )

points(u1A$cumS, u1A$respsd, pch = 15, col ='purple3')
points(u2A$cumS, u2A$respsd, pch = 15, col = 'purple3')
points(u1B$cumS, u1B$respsd, pch = 16, col = 'gray30')
points(u2B$cumS, u2B$respsd, pch = 16, col = 'gray30')
lines(u1A$cumS, u1A$respsd, lty = 1, lwd = 2, col ='purple3')
lines(u2A$cumS, u2A$respsd, lty = 4, lwd = 2, col = 'purple3')
lines(u1B$cumS, u1B$respsd, lty = 1, lwd = 2, col = 'gray30')
lines(u2B$cumS, u2B$respsd, lty = 4, lwd = 2, col = 'gray30')
points(d1A$cumS, d1A$respsd, pch = 15, col ='purple3')
points(d2A$cumS, d2A$respsd, pch = 15, col = 'purple3')
points(d1B$cumS, d1B$respsd, pch = 16, col = 'gray30')
points(d2B$cumS, d2B$respsd, pch = 16, col = 'gray30')
lines(d1A$cumS, d1A$respsd, lty = 1, lwd = 2, col ='purple3')
lines(d2A$cumS, d2A$respsd, lty = 4, lwd = 2, col = 'purple3')
lines(d1B$cumS, d1B$respsd, lty = 1, lwd = 2, col = 'gray30')
lines(d2B$cumS, d2B$respsd, lty = 4, lwd = 2, col = 'gray30')
legend('topleft', legend = c('Directional A', 'Direct + Stab A', 'Directional B', 'Direct + Stab B'), 
       col = c('purple3','purple3','gray30','gray30'), lty = c(1,4,1,4), lwd = c(2,2,2,2), cex = 0.9)

## calculating heritability from cumS and the response through linear regression
x1A <- c(d1A$cumS[-1],u1A$cumS)
y1A <- c(d1A$resp[-1], u1A$resp)
x2A <- c(d2A$cumS[-1],u2A$cumS)
y2A <- c(d2A$resp[-1], u2A$resp)
x1B <- c(d1B$cumS[-1],u1B$cumS)
y1B <- c(d1B$resp[-1], u1B$resp)
x2B <- c(d2B$cumS[-1],u2B$cumS)
y2B <- c(d2B$resp[-1], u2B$resp)

m1A <- lm(y1A ~ x1A)
m2A <- lm(y2A ~ x2A)
m1B <- lm(y1B ~ x1B)
m2B <- lm(y2B ~ x2B)

herit.info <- function(x) {
h2m <- x$coefficients[2]
st.e <- summary(x)$coeff[[4]]
conf.high <- h2m + 1.96*st.e
conf.low  <- h2m - 1.96*st.e
output1 <- paste0('[',round(conf.low,3),'-',round(conf.high,3),']' )
output <- paste(round(h2m, 3), output1, sep = ' ')
return(output)
}

h1A <- herit.info(m1A)
h2A <- herit.info(m2A)
h1B <- herit.info(m1B)
h2B <- herit.info(m2B)

if (species == 'Sim') {
  text(0.068, -4, paste('Heritability [95% CI]', '\n' ,'Directional A  = ', h1A , '\n','Direct + Stab A = ', h2A, '\n', 'Directional B  = ', h1B, '\n', 'Direct + Stab B = ', h2B) , cex = 0.95)
}
if (species == 'Mel') {
  text(0.06, -4.3, paste('Heritability [95% CI]', '\n' ,'Directional A  = ', h1A , '\n','Direct + Stab A = ', h2A, '\n', 'Directional B  = ', h1B, '\n', 'Direct + Stab B = ', h2B) , cex = 0.95)
}



################################################################

wdatA <- wngdat[wngdat$rep == 'a', ]
wcontrolA <- wdatA[wdatA$treatment == 'control',]
wd1A  <- wdatA[wdatA$treatment == 'down1', ]
wd2A  <- wdatA[wdatA$treatment == 'down2', ]
wu1A  <- wdatA[wdatA$treatment == 'up1', ]
wu2A  <- wdatA[wdatA$treatment == 'up2', ]
wdatB <- wngdat[wngdat$rep == 'b', ]
wcontrolB <- wdatB[wdatB$treatment == 'control',]
wd1B  <- wdatB[wdatB$treatment == 'down1', ]
wd2B  <- wdatB[wdatB$treatment == 'down2', ]
wu1B  <- wdatB[wdatB$treatment == 'up1', ]
wu2B  <- wdatB[wdatB$treatment == 'up2', ]

##### Find the range of the response and the number of generations
wa_resp_range <- range(wdatA$wngdev)
wb_resp_range <- range(wdatB$wngdev)
wrng <- c(wa_resp_range, wb_resp_range)
wrng_gen <- range(as.numeric(wngdat$gen))

### This makes a plot of the response by the generation
plot(1, xlim = c(min(wrng_gen), max(wrng_gen)), ylim = c(min(wrng), max(wrng)), type = 'n', xlab = 'Generation', 
     ylab = 'Mean Wing Deviation', cex.lab = 0.75,cex.axis=0.75,las = 1, main = substitute(italic(x), list(x=title)) )
points(wu1A$gen, wu1A$wngdev, pch = 15, col ='purple3')
points(wu2A$gen, wu2A$wngdev, pch = 15, col = 'purple3')
points(wu1B$gen, wu1B$wngdev, pch = 16, col = 'gray30')
points(wu2B$gen, wu2B$wngdev, pch = 16, col = 'gray30')
points(wcontrolA$gen, wcontrolA$wngdev, pch = 18, col ='black')
points(wcontrolB$gen, wcontrolB$wngdev, pch = 18, col = 'red')
lines(wu1A$gen, wu1A$wngdev, lty = 1, lwd = 2, col ='purple3')
lines(wu2A$gen, wu2A$wngdev, lty = 6, lwd = 2, col = 'purple3')
lines(wu1B$gen, wu1B$wngdev, lty = 1, lwd = 2, col = 'gray30')
lines(wu2B$gen, wu2B$wngdev, lty = 6, lwd = 2, col = 'gray30')
lines(wcontrolA$gen, wcontrolA$wngdev, lty = 1, lwd = 2, col = 'black')
lines(wcontrolB$gen, wcontrolB$wngdev, lty = 1, lwd = 2, col = 'red')
points(wd1A$gen, wd1A$wngdev, pch = 15, col ='purple3')
points(wd2A$gen, wd2A$wngdev, pch = 15, col = 'purple3')
points(wd1B$gen, wd1B$wngdev, pch = 16, col = 'gray30')
points(wd2B$gen, wd2B$wngdev, pch = 16, col = 'gray30')
lines(wd1A$gen, wd1A$wngdev, lty = 1, lwd = 4, col ='purple3')
lines(wd2A$gen, wd2A$wngdev, lty = 6, lwd = 4, col = 'purple3')
lines(wd1B$gen, wd1B$wngdev, lty = 1, lwd = 4, col = 'gray30')
lines(wd2B$gen, wd2B$wngdev, lty = 6, lwd = 4, col = 'gray30')
legend('topleft', legend = c('Directional A', 'Direct + Stab A', 'Directional B', 'Direct + Stab B', 'Control A', 'Control B'), 
       col = c('purple3','purple3','gray30','gray30', 'black', 'red'), lty = c(1,6,1,6,1,1), lwd = c(2,2,2,2,2,2), cex = 0.78)


##### Find the range of the response and the number of generations
wa_resp_range <- range(wdatA$var)
wb_resp_range <- range(wdatB$var)
wrng <- c(wa_resp_range, wb_resp_range)

### This makes a plot of the response by the generation
plot(1, xlim = c(min(wrng_gen), max(wrng_gen)), ylim = c(min(wrng), max(wrng)), type = 'n', xlab = 'Generation', 
     ylab = 'Variance in Wing Deviations', cex.lab = 0.75,cex.axis=0.75,las = 1, main = substitute(italic(x), list(x=title)) )
points(wu1A$gen, wu1A$var, pch = 15, col ='purple3')
points(wu2A$gen, wu2A$var, pch = 15, col = 'purple3')
points(wu1B$gen, wu1B$var, pch = 16, col = 'gray30')
points(wu2B$gen, wu2B$var, pch = 16, col = 'gray30')
points(wcontrolA$gen, wcontrolA$var, pch = 18, col ='black')
points(wcontrolB$gen, wcontrolB$var, pch = 18, col = 'red')
lines(wu1A$gen, wu1A$var, lty = 1, lwd = 2, col ='purple3')
lines(wu2A$gen, wu2A$var, lty = 6, lwd = 2, col = 'purple3')
lines(wu1B$gen, wu1B$var, lty = 1, lwd = 2, col = 'gray30')
lines(wu2B$gen, wu2B$var, lty = 6, lwd = 2, col = 'gray30')
lines(wcontrolA$gen, wcontrolA$var, lty = 1, lwd = 2, col = 'black')
lines(wcontrolB$gen, wcontrolB$var, lty = 1, lwd = 2, col = 'red')
points(wd1A$gen, wd1A$var, pch = 15, col ='purple3')
points(wd2A$gen, wd2A$var, pch = 15, col = 'purple3')
points(wd1B$gen, wd1B$var, pch = 16, col = 'gray30')
points(wd2B$gen, wd2B$var, pch = 16, col = 'gray30')
lines(wd1A$gen, wd1A$var, lty = 1, lwd = 2, col ='purple3')
lines(wd2A$gen, wd2A$var, lty = 6, lwd = 2, col = 'purple3')
lines(wd1B$gen, wd1B$var, lty = 1, lwd = 2, col = 'gray30')
lines(wd2B$gen, wd2B$var, lty = 6, lwd = 2, col = 'gray30')
legend('topleft', legend = c('Directional A', 'Direct + Stab A', 'Directional B', 'Direct + Stab B', 'Control A', 'Control B'), 
       col = c('purple3','purple3','gray30','gray30', 'black', 'red'), lty = c(1,6,1,6,1,1), lwd = c(2,2,2,2,2,2), cex = 0.78)


