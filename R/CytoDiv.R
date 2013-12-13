cytoDiv <- function(df, para, Ncat = 16, log=TRUE, do.plot=FALSE){
	
	jet.colors <- colorRampPalette(c("white","red","red4","tomato4","black"))
	jet.colors <- colorRampPalette(c("blue4","royalblue4","deepskyblue3", "seagreen3", "yellow", "orangered2","darkred"))
	
	# Create log spaced grid of N categories
	
	if(log==TRUE){ 
		# limits of the different categories
		diff <- log(max(df[,para], na.rm=T), base=2)/(Ncat+2)
		cut.ext <- 2^((1:(Ncat+3) -1)*diff)
		breaks <- cut.ext[-c(1, Ncat+3)]
		
		# labels of the different categories
		cat.diff <- log(max(df[,para], na.rm=T), base=2)/(Ncat-1)
		categories <- 2^((1:(Ncat) -1)*cat.diff)
		}

	# Create lin spaced grid of N categories
	
	if(log==FALSE){ 
		# limits of the different categories
		cut.ext <- seq(1, max(df[,para], na.rm=T), length=Ncat+2)
		breaks <- cut.ext[-c(1, Ncat+2)]
		
		# labels of the different categories
		categories <- seq(1, max(df[,para], na.rm=T), length=Ncat)
		}
	
	m <- as.matrix(df[,para])

	if(dim(m)[2] == 1){
		mints <- data.frame(m1=cut(m[,para[1]],breaks=breaks, labels=FALSE))
		kk <- as.matrix(mints)
		
		# construct 4-d matrix of param space
		p <-array(0,c(Ncat))
				
		for(i in 1:dim(kk)[1]){
			k <-c(kk[i,])
			p[k[dim=1]]=p[k[dim=1]]+1
		#	print(i) #for debugging
			rm(k)
		}
		p <- p/sum(p) # calculate frequency
		
	if(do.plot==T & log==T) barplot(p, categories, log='x', main=paste("parameter =",para[1],"\nN =",Ncat), ylab="Probability", xlab=paste(para[1]), col='lightgrey', names.arg=round(categories),axes=T)
	if(do.plot==T & log==F) barplot(p, categories, main=paste("parameter =",para[1],"\nN =",Ncat), ylab="Probability", xlab=paste(para[1]), col='lightgrey', names.arg=round(categories),axes=T)
	
	}	

	if(dim(m)[2] == 2){
		mints <- data.frame(m1=cut(m[,para[1]],breaks=breaks, labels=FALSE),m2=cut(m[,para[2]],breaks=breaks, labels=FALSE))
		kk <- as.matrix(mints)
		
		# construct 4-d matrix of param space
		p <-array(0,c(Ncat,Ncat))
		
		for(i in 1:dim(kk)[1]){
			k <-c(kk[i,])
			p[k[dim=1],k[dim=2]]=p[k[dim=1],k[dim=2]]+1
		#	print(i) #for debugging
			rm(k)
		}
	
	p <- p/sum(p) # calculate frequency

	if(do.plot==T & log==T) persp(x=log(categories), y=log(categories), z=p,xlab=paste(para[1]), ylab=paste(para[2]),zlab="Probability", main=paste("N =",Ncat), ticktype='detailed', phi=45, theta=-35,col='lightgrey', shade=0.25)
		

	if(do.plot==T & log==F) persp(x=categories, y=categories, z=p,xlab=paste(para[1]), ylab=paste(para[2]),zlab="Probability", main=paste("N =",Ncat), ticktype='detailed', phi=45, theta=-35,col='lightgrey', shade=0.25)

	}
	
	if(dim(m)[2] == 3){
		mints <- data.frame(m1=cut(m[,para[1]],breaks=breaks, labels=FALSE),m2=cut(m[,para[2]],breaks=breaks, labels=FALSE),m3=cut(m[,para[3]],breaks=breaks, labels=FALSE))
		kk <- as.matrix(mints)
		
		# construct 4-d matrix of param space
		p <-array(0,c(Ncat,Ncat,Ncat))
			for(i in 1:dim(kk)[1]){
			k <-c(kk[i,])
			p[k[dim=1],k[dim=2],k[dim=3]]=p[k[dim=1],k[dim=2],k[dim=3]]+1
		#	print(i) #for debugging
			rm(k)
		}
		
		p <- p/sum(p) # calculate frequency
	
	if(do.plot==T ) print("No plot, too many dimensions ...")
		
	}
	
	if(dim(m)[2] == 4){
		mints <- data.frame(m1=cut(m[,para[1]],breaks=breaks, labels=FALSE),m2=cut(m[,para[2]],breaks=breaks, labels=FALSE),m3=cut(m[,para[3]],breaks=breaks, labels=FALSE), m4=cut(m[,para[4]],breaks=breaks, labels=FALSE))
		kk <- as.matrix(mints)
		
		# construct 4-d matrix of param space
		p <-array(0,c(Ncat,Ncat,Ncat,Ncat))

		for(i in 1:dim(kk)[1]){
			k <-c(kk[i,])
			p[k[dim=1],k[dim=2],k[dim=3],k[dim=4]]=p[k[dim=1],k[dim=2],k[dim=3],k[dim=4]]+1
		#	print(i) #for debugging
			rm(k)
		}
	
	if(do.plot==T ) print("No plot, too many dimensions ...")

	
	}
	
	
	
	# vectorize and remove null values from p matrix
	p <- cbind(p)

	# calculate diversity indices
	p_i = p[p>0]/sum(p[p>0]) # remove empty category
	N0 = sum(p_i^0) # Number of categories
	N1 = exp(-sum(p_i*log(p_i))) 
	H = log(N1) # Shannon-Wiener's Index - use natural log
	J = log(N1)/log(N0) # Evenness
	
	# output
	indices <- data.frame(cbind(N0, H, J))
	return(indices)
	
	# p <- pi[!(pi==0)] 
	# N0 <- sum(p^0) # Number of categories
	# H <- - sum(p*log(p)) # Shannon-Wiener Diversity index H'
	# N2 <- sum(p^2) # Simpson's index
	# D <- 1 - N2 # Simpson's Index of Diversity
	# J <- H/log(N0) # Evenness
	
	
}