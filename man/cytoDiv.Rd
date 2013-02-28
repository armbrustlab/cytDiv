\name{cytoDiv}
\alias{cytoDiv}
\title{Cytometric Diversity Indices}
\description{
Calculate the various cytometric diversity indices described in Li, W. 1997 Cytometric diversity in marine ultraphytoplankton. Limnology and Oceanography 42: 874 - 880. The proportional abundance of the ith category of organism is calculated using the KernSur function from GenKern package, which computes bivariate kernel density estimates.
}
\usage{
cytoDiv(df, x.var="fsc_small", y.var="chl_small", Ncat = 16, Nbit = 16, do.plot=FALSE, ...)
}
\arguments{
\item{df}{a flow cytometry file converted into a dataframe.}
\item{x.var}{name of the variable x. Default is "fsc_small" from SeaFlow instrument}
\item{y.var}{name of the variable y. Default is "chl_small" from SeaFlow instrument}
\item{Ncat}{Number of categories per variable. Default is 16}
\item{Nbit}{Number of the bit resolution the Analog-to-Digital Converter. Default is 16, i.e values range from 0 to 2^16}
\item{do.plot}{Option to plot the probability distribution of the bivariate data in 3D}
\item{...}{additional arguments to be passed to the KernSur function from the GenKern package}
}
\value{
A data frame that contains the Richness index N0 (number of categories), the Shannon-Wiener Diversity index H', the reciprocal of Simpson's index N2, the Simpson's Index of Diversity D and the Evenness index J'.
Indices calculation is based on Hill's method (Hill, M.O. 1973 Diversity and evenness: A unifying notation and its consequences. Ecology 54: 427 - 432). 
}
\note{
The highest number of distinct categories (Ncat) in cytometric classification depends on the resolution at which each variable is measured, which is controlled by the bit resolution (Nbit) of the Analog-to-Digital Converter. In the example below, measurements of single-cell forward light scatter and red fluorescence were collected using logarithmic amplification and recorded in relative units in a four-decade range spanned by 2^16 channels. For the calculation of the diversity indices, data resolution of each variable was reduced to 16 channels by successively binning the counts of 2^12 adjacent channels . In this case, the maximum number of possible categories in the light scatter-fluorescence domain was therefore 16x16 (=256). 
}
\examples{

## reading from a fcs dataframe

fcs.df <- system.file("extdata","fcs_dataframe.csv", 
				package="cytoDiv")
df <- read.csv(fcs.df)

## looking at the first rows of the data frame
head(df)

## plotting bivariate plot
rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))
plot(df[,2], df[,3], col= densCols(df[,2], df[,3], colramp = rainbow.cols))

# Remove the Internal Standard and noise before calculating the Cytometric Diversity Indices
cleaned.df <- subset(df, df[,4]  < df[,3] + 13000 & df[,3] > 5000)
plot(cleaned.df[,2], cleaned.df[,3], col= densCols(cleaned.df[,2], cleaned.df[,3], colramp = rainbow.cols))

#Calculating the Cytometric Diversity Indices
div <- cytoDiv(cleaned.df, x.var="FSC", y.var="F692.4", do.plot=TRUE)

print(div)

}