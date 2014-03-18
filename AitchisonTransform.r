
# Ruth Wong
# Gloor Lab
# ruthgracewong@gmail.com

# run as follows:
# yoursummedtransformeddata <- aitchison.transform.reads(PARAMETERS)
# the aitchison.transform.reads function is defined as follows:
# aitchison.transform.reads <- function(filename="formatted_readcounts_subsyshier.txt",rounded=TRUE, subjects = -1, firstsubjectindex = 3, lastsubjectindex = -1, groupindex = -1,lengthindex=2,outputfolder="")

# two tables will be written to the outputfolder (which is relative to the current directory):
# AitchisonTransformedData.txt
# AitchisonTransformedDataForALDExInput.txt



#written by Andrew Fernandes
#perform Aitchison transform on one subject (ie. set of counts from one sample)
diffex.inputformat.aitchison.mean <- function( n, log=FALSE ) {
# Input is a vector of non-negative integer counts.
# Output is a probability vector of expected frequencies.
# If log-frequencies are requested, the uninformative subspace is removed.
	
# Example usage: If we observe 5 heads and 3 tails from coin-flipping,
#                then the expected proportion of heads versus tails
#                is 'aitchison.mean( c(5,3) )'.
	
    n <- round( as.vector( n, mode="numeric" ) )
    if ( any( n < 0 ) ) stop("counts cannot be negative")

#add pseudocount
    a <- n + 0.5
    sa <- sum(a)

#digamma is log space
#sample from digamma dist
    log.p <- digamma(a) - digamma(sa)
#subspace removal
    log.p <- log.p - mean(log.p)
		
    if ( log ) return(log.p)
	
    p <- exp( log.p - max(log.p) )
    p <- p / sum(p)
    return(p)
}




aldexformat <- function(subsys4,subjects) {
		subsys4aldexformat <- data.frame(matrix(nrow=nrow(subsys4),ncol=subjects))
		
		subsys4aldexformat <- subsys4[, 3:ncol(subsys4)]
		rownames(subsys4aldexformat) <- subsys4[, 1]
		colnames(subsys4aldexformat) <- colnames(subsys4)[3:ncol(subsys4)]
	return(subsys4aldexformat)
}



#sum Aitchison transformed data for each grouping
diffex.inputformat.extractAitchisontotals <- function(holder,firstsubjectindex,lastsubjectindex,groupindex){
	diffreads <- holder

		subjects <- colnames(diffreads)[firstsubjectindex:(lastsubjectindex)]
	#get columns of original data
	subject_cols <- c(firstsubjectindex:(lastsubjectindex))
	#get non-duplicated groupings
	grouping <- diffreads[which(!duplicated(diffreads[, (groupindex)])) , (groupindex)]
	#turn NA into Not Available
	grouping[which(is.na(grouping))] <- "Not Available"
	#create data structure for grouping totals
	subsys4 <- data.frame(matrix(ncol=0,nrow=length(grouping)))
	#add grouping categories in first column
	subsys4[1] <- grouping
	#assign original grouping column name (in case it contains information about grouping type)
	colnames(subsys4) <- colnames(diffreads)[groupindex]
	#assign indices of first instance of each grouping
	subsys4$indices <- which(!duplicated(diffreads[, (groupindex)]))
	
	
	
	
	
	#get sample names
	subjects <- colnames(diffreads)[firstsubjectindex:(lastsubjectindex)]

	#get indices of Aitchison-transformed data
	diffreads_transformed_subjects_col <- which(colnames(diffreads) %in% subjects)
	sort(diffreads_transformed_subjects_col)
	
	#create data structure for original data with Aitchison transformed sums appeneded in the last columns
	newsubsys4 <- data.frame(matrix(ncol = length(colnames(subsys4)) + length(subjects), nrow = length(rownames(subsys4))))
	colnames(newsubsys4) <- c(colnames(subsys4),paste(subjects,"sum",sep="_"))
	
	#get indices of where Aitchison-transformed data needs to be transferred
	subsys4_transformed_sum_col <- which(colnames(newsubsys4) %in% paste(subjects,"sum",sep="_"))
	sort(subsys4_transformed_sum_col)

	#insert old data into new data frame
	newsubsys4[, 1:length(colnames(subsys4))] <- subsys4

	#extract Aitchison transform sums
	for (i in 1:length(newsubsys4$indices)) {
		#if this is the last grouping, loop to end of data to find index of last Aitchison data of the last grouping
		if (is.na(newsubsys4$indices[i+1])) {
			lastIndex <- newsubsys4$indices[i]
			while (!is.na(diffreads[lastIndex,(groupindex)])) {
				lastIndex = lastIndex + 1
			}
			lastIndex <- lastIndex - 1
			newsubsys4 <- newsubsys4[1:i , ]
			i <- length(newsubsys4$indices)
		}
		#otherwise last index of Aitchison data for current grouping is one before the first index of data for the next grouping
		else {
			lastIndex = newsubsys4$indices[i+1]-1
		}
		#sum Aitchison data for current grouping, for each subject
		for (j in 1:length(subsys4_transformed_sum_col)){
			newsubsys4[i, subsys4_transformed_sum_col[j]] <- sum(diffreads[newsubsys4$indices[i]:lastIndex, diffreads_transformed_subjects_col[j]],na.rm = TRUE)
		}
	}
	return(newsubsys4)
}







#main function
# filename is the read count file, rounded is for whether the output should be integers or not,
# subjects is the number of subjects (the number of columns of read counts)
# firstsubjectindex, lastsubjectindex, lowest and highest indices of the columns with subject read count data - subjects are assumed to be consecutive
# groupcolumn is the column for the grouping data, length column is the column for the lengths of the reads
# outputfolder is where the transformed counts will be placed
aitchison.transform.reads <- function(filename="formatted_readcounts_subsyshier.txt",rounded=TRUE, subjects = -1, firstsubjectindex = 3, lastsubjectindex = -1, groupindex = -1,lengthindex=2,outputfolder=""){

	if (outputfolder !="") {
		outputfolder = paste("./",outputfolder,"/",sep="")
	}
	
	# file should have columns for refseq_id, length, subjects, and grouping, in that order
	originaldata <- read.table(filename, header = TRUE, sep= "\t", stringsAsFactors=F, quote = "", check.names = FALSE, comment.char = "")

	# d is data that will be messed with. Original data is left alone, just in case.
	d <- originaldata
	print("data successfully read in")
	
	#setting variables with information about the number and indices of the subjects

	if (lastsubjectindex == -1) {
		lastsubjectindex = ncol(d) - 1
	}
	if (subjects == -1) {
		subjects = lastsubjectindex-firstsubjectindex + 1
	}
	if (groupindex == -1) {
		groupindex = ncol(d)
	}

	print(paste("number of subjects is",subjects))

		# Aitchison transform, includes summing data to 1
		
		
		d[, firstsubjectindex:lastsubjectindex] <- apply(d[, firstsubjectindex:lastsubjectindex], 2, function(x) diffex.inputformat.aitchison.mean(x))
		print("data weighted by Aitchison transform")
		
		
		# length weight
		d[, c(firstsubjectindex:lastsubjectindex)] <- t(as.data.frame(apply(d[, c(firstsubjectindex:lastsubjectindex, lengthindex)], 1, function(x) x[1:subjects]/x[(subjects+1)])))
		print("Aitchison-transformed data length-weighted")

		#extract unique refseq_id/grouping combinations such that refseq_id which have the same grouping of interest (but are different for another grouping in the hierarchy) are not duplicated. Refseq_id which have multiple groupings of interest remain duplicated.
			diffreads <- d[which(!duplicated(paste(d[, (groupindex)],d[, 1],sep = "\t"))) ,]
			
	print("extracted rows with unique subsys4/refseq_id combinations")
	diffreads <- diffreads[order(diffreads[, groupindex]) , ]
	print("unique subsys/refID combos sorted by subsys4")
	
		print("extracting Aitchison transformed totals")
		subsys4 <- diffex.inputformat.extractAitchisontotals(diffreads,firstsubjectindex,lastsubjectindex,groupindex)
		print("done extracting Aitchison transformed totals.")
		
		
		columnsums <- apply(subsys4[3:ncol(subsys4)], 2, sum)
		
		subsys4[3:ncol(subsys4)] <- t(as.data.frame(apply(subsys4[3:ncol(subsys4)], 1, function(x) x/columnsums)))
		
		totalreads <- apply(originaldata[, firstsubjectindex:lastsubjectindex], 2, sum)
		
		subsys4[3:ncol(subsys4)] <- t(as.data.frame(apply(subsys4[3:ncol(subsys4)], 1, function(x) x*totalreads)))
		
		if (rounded) {
			subsys4[3:ncol(subsys4)] <- apply(subsys4[3:ncol(subsys4)],1:2,round)
		}
		
		write.table(subsys4[, c(1,c(3:ncol(subsys4)))], file = paste(outputfolder,"AitchisonTransformedData.txt",sep=""), row.names=FALSE, append = FALSE, quote = FALSE, sep = "\t")
		
		subsys4aldexformat <- aldexformat(subsys4,subjects)

		write.table(subsys4aldexformat, file = paste(outputfolder,"AitchisonTransformedDataForALDExInput.txt",sep=""), col.names=NA, append = FALSE, quote = FALSE, sep = "\t")
		
	return(subsys4)
}
