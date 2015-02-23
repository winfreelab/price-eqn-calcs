# PRICE EQUATION "SPACE" CALCULATIONS
# this R code goes with the Winfree et al. (2015) Ecology Letters paper.  Please refer to this paper for further details, or as a citation for the use of this code:
# Winfree R, Fox J, Williams N, Reilly JR, Cariveau D (2015) Abundance of common species, not species richness, drives delivery of a real-world ecosystem service. Ecology Letters (in press)
# this version of the code was written by James Reilly, Feb 2015

# to use this code with your own data, you will need to replace "EXAMPLE_CROP" with your own system name and create the following three input data files:
	#1) an individuals data file for each year with sites as the rows and species as the columns, named like this: "price_EXAMPLE_CROP_individuals_2010.csv"
	#2) an data file with single-visit pollen deposition observations by pollinator ID code: "price_EXAMPLE_CROP_pollen_IDcode.csv"
	#3) a correspondence table that relates species to pollinator ID codes: "price_EXAMPLE_CROP_species_IDcode.csv"
# you will also need to adjust the code so the years match those of your study on lines 39-60.

# Note on a common error:
# If there is an error like: "Error in baseline$pollen[j] = as.numeric(pollen_summary[subset(speciesID,  : replacement has length zero",
# there is probably a data mismatch between speciesID file and individuals file.  Run this to find which species did it:  splist[j]

###############

# required R packages must be installed:
library(reshape)
library(plyr)

# this code runs the "space" version of the price calculations:
this_code = "space"

this_study = "EXAMPLE_CROP"


# set whether to normalize the output here:
normalize = "TRUE"

pollen_stat = "mean"	# for sensitivity analysis of single visit pollen data--e.g. you can use q40 instead ~line 110
pollen_multiplier = 1 	# for sensitivity analysis of single visit pollen data


############################
# Data input
############################

	# specify here which individuals_year files you want to use:
	individuals_2010 <- read.csv(paste("R data/price_", this_study, "_individuals_2010.csv", sep=""), header=TRUE)
	individuals_2011 <- read.csv(paste("R data/price_", this_study, "_individuals_2011.csv", sep=""), header=TRUE)
	individuals_2012 <- read.csv(paste("R data/price_", this_study, "_individuals_2012.csv", sep=""), header=TRUE)

	individuals_datasets = list(individuals_2010,individuals_2011,individuals_2012)

	# make a list of years -- do this by hand:
	yearlist = c(2010, 2011, 2012)

	pollen_data <-  read.csv(paste("R data/price_", this_study, "_pollen_IDcode.csv", sep=""), header=TRUE)
	speciesID <-  read.csv(paste("R data/price_", this_study, "_species_IDcode.csv", sep=""), header=TRUE)

	# get a list of sites from one of the individuals datafiles (here getting it from the 2010 data):
	sitelist = as.character(individuals_2010$farm)

	# get a list of all species from the individuals files:
	splist = unique(c(
			names(individuals_2010)[2:ncol(individuals_2010)],
			names(individuals_2011)[2:ncol(individuals_2011)],
			names(individuals_2012)[2:ncol(individuals_2012)]
			))

############################
# Main program loops
############################

# test run loop calculates most functional sites and re-runs program using those as the baseline for that year
for (test_run in c(1,0)) {

	if (test_run==1) {
		suggested_sites = rep(sitelist[1], length(yearlist))
	} else {
		suggested_sites = suggestions$suggested_site
	}
	
	# set up an individuals dataframe with data for all species in all years at all farms -- there will be zeros when the species wasn't recorded.
	individuals_data = data.frame()

	for (y in 1:length(yearlist)) {

		individuals_temp = as.data.frame(matrix(0,length(sitelist),length(splist)+2))
		names(individuals_temp) = c("farm","year",splist)

		individuals_temp$farm = sitelist
		individuals_temp$year = rep(yearlist[y],length(sitelist))

		for (x in 1:length(splist)) {
			if (length(subset(names(individuals_datasets[[y]]), names(individuals_datasets[[y]])==as.character(splist[x])))>0) {
				individuals_temp[splist[x]] = individuals_datasets[[y]][splist[x]]
			} else {
				individuals_temp[splist[x]] = rep(0,length(sitelist))
			}
		}
		
		individuals_data = rbind(individuals_temp, individuals_data)
	}
	

	# get mean, etc for pollen grains by IDcode:
	mean_narm = function(x) {
		mean(x,na.rm=T)
	}
	q40 = function(x) {
		quantile(x,probs=c(0.4),na.rm=T)
	}
	q60 = function(x) {
		quantile(x,probs=c(0.6),na.rm=T)
	}

	pollen_summary = data.frame(
		mean = tapply(pollen_data$pollen_grains,pollen_data$IDcode,mean_narm),
		q40 = tapply(pollen_data$pollen_grains,pollen_data$IDcode,q40),
		q60 = tapply(pollen_data$pollen_grains,pollen_data$IDcode,q60)
	)	
	
	# set baseline site AND year here: (note: program calculates this itself)
	baseline_site = as.character(suggested_sites)

	# set up a master list to put the output dataframes in from each year
	output_list = list()

	# loop through the years:
	for (g in 1:length(yearlist)) {

		baseline_year = yearlist[g]

		# first calculate species presence and function for the baseline site.  This is done outside the loop of other sites, since they all need to be compared to this information
		baseline = data.frame(species = splist) 
		
		# create baseline data frame from individuals data, using names from the column headers in the individuals data.
		baseline_row = subset(individuals_data, farm==baseline_site[g] & year==baseline_year) 
		
		# look only at data for the baseline farm
		baseline$individuals = as.numeric(baseline_row[3:length(baseline_row)])
		
		# get the pollen per visit number for each species:
		for (j in 1:nrow(baseline)) {
			baseline$pollen[j] = 
			as.numeric(
				pollen_summary[
					as.character(subset(speciesID, genus_species==as.character(baseline$species[j]))$IDcode),
					pollen_stat
				]
			)*pollen_multiplier
		}
		
		baseline$func = baseline$individuals * baseline$pollen

		baseline_total_func = sum(baseline$func)
		baseline_total_species = nrow(subset(baseline, individuals>0))
		baseline_mean_func = baseline_total_func / baseline_total_species
		
		# set up a data frame to write the results into
		output = data.frame()

		# create data frame to store total function of each site
		sitefuncs = data.frame(farm=sitelist, func=rep(0,length(sitelist)), baseline_func=rep(0,length(sitelist)))

		# loop through the sites, calculating function, counting present species, and comparing to baseline

		for (h in 1:length(sitelist)) {
			
			# for each site, we create a temporary dataframe called "comparison", that we store the information in while we are comparing to baseline
			comparison = data.frame(species = splist)
			comparison_row = subset(individuals_data, farm==sitelist[h] & year==yearlist[g])

			comparison$individuals = as.numeric(comparison_row[3:length(comparison_row)])

			# get the mean pollen per visit for each species:
			for (k in 1:nrow(comparison)) {
				comparison$pollen[k] = 
				as.numeric(
					pollen_summary[
						as.character(subset(speciesID, genus_species==as.character(baseline$species[k]))$IDcode),
						pollen_stat
					]
				)*pollen_multiplier
			}

			comparison$func = comparison$individuals * comparison$pollen
			comparison_total_func = sum(comparison$func)
			comparison_total_species = nrow(subset(comparison, individuals>0))
			comparison_mean_func = comparison_total_func / comparison_total_species
			
			# compare each species row in the comparison site to that species row in the baseline site
			for (i in 1:nrow(comparison)) {
				# Is the species present at both sites?
				if (comparison$individuals[i] > 0 & baseline$individuals[i] > 0) {
					# The "common" column will be 1 if present at both, or 0 if not
					comparison$common[i] = 1
					comparison$zc[i] = baseline$func[i]
					comparison$zcprime[i] = comparison$func[i]
				} else {
					comparison$common[i] = 0
					comparison$zc[i] = NA
					comparison$zcprime[i] = NA
				}
			}

			# count up common species
			common_total_species = sum(comparison$common)
			# calculate mean for baseline species
			zc_mean = mean(comparison$zc, na.rm=T)
			zcprime_mean = mean(comparison$zcprime, na.rm=T)
				
			RICH_L = (common_total_species - baseline_total_species) * baseline_mean_func
			RICH_G = (comparison_total_species - common_total_species) * comparison_mean_func
			COMP_L = (zc_mean - baseline_mean_func) * common_total_species
			COMP_G = -(zcprime_mean - comparison_mean_func) * common_total_species
			ABUN = (zcprime_mean - zc_mean) * common_total_species
			T_T = comparison_total_func - baseline_total_func

			output_line = data.frame(comparison_site=sitelist[h], RICH_L=RICH_L, RICH_G=RICH_G, COMP_L=COMP_L, COMP_G=COMP_G, ABUN=ABUN,T_T=T_T)
			# add output from this site into the output dataframe
			output = rbind(output_line,output)
			
			sitefuncs$func[h] = comparison_total_func
			sitefuncs$baseline_func[h] = baseline_total_func
	
		} # end site loop

		# remove baseline site from table, since it is meaningless to compare it to itself
		output = subset(output, comparison_site!=baseline_site[g])


	# save the output into the master list
	output_list[[g]] = output

	mostfunctionalsite = as.character(sitefuncs$farm[which(sitefuncs$func==max(sitefuncs$func))])
	suggested_sites[g] = mostfunctionalsite

	} # end year loop

	names(output_list) = yearlist

	output_list
	suggestions = data.frame(year = yearlist, suggested_site = suggested_sites)

} # end test run loop

#suggestions

results = data.frame()
for (i in 1:length(yearlist)) {
	results = rbind(cbind(names(output_list)[i], output_list[[i]]), results)
}
names(results)[1] = "year"

# write the output to csv if desired:
write.table(results, file="price space output.csv", sep=",", col.names=T, row.names=F)
#results

# get rid of site-years in which any price term returns an NA:
results2 = results
results2$ALL = results2$RICH_L + results2$RICH_G + results2$COMP_L + results2$COMP_G + results2$ABUN
results2 = subset(results2, ALL!="NaN")

results2 = subset(results2, select=-c(T_T, ALL))

############################
# get *combined price term* results:
############################

results2$RICH_L_COMP_L = results2$RICH_L + results2$COMP_L
results2$RICH_G_COMP_G = results2$RICH_G + results2$COMP_G
results2$RICH_L_RICH_G_COMP_L_COMP_G = results2$RICH_L + results2$RICH_G + results2$COMP_L + results2$COMP_G


############################
# Normalize results
############################

# note we are normalizing across all 5 terms and summed terms so summed terms cannot exceed 1 and can be compared in magnitude to single terms

if (normalize == "TRUE") {
	results2$normalizer = apply(abs(results2[,3:10]), MARGIN=1, FUN=max)

	results2[,3:10] = results2[,3:10]/results2$normalizer
	results2 = subset(results2, select=-normalizer)
}

output2 = melt(results2, id=1:2, variable_name="price_term")
output2$price_term = factor(output2$price_term)

############################
# Summarize
############################

# print to screen summary stats for the normalized results, that is, summary stats behind the box plots:
price_summary = aggregate(data=output2, value ~ price_term + year, FUN=function(x) c(mean =mean(x), median=median(x), sd=sd(x)))
price_summary$value = round(price_summary$value,3)
price_summary


