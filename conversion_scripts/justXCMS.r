library(xcms)
mzXMLfiles <- dir(full.names=TRUE,pattern="\\.mzXML")
for(inName in mzXMLfiles) {
	print(inName)
	outname = sub(".mzXML", ".xcms.csv", inName)
	print(outname)
	xset <- xcmsSet(inName, method='centWave', 
		ppm=2,
		peakwidth=c(10, 100), 
		snthresh=5,
		prefilter=c(3, 1000), 
		integrate=1,
		mzdiff=0.001,
		fitgauss="FALSE",
		verbose.columns=TRUE)
	write.csv(xset@peaks[, c("mz", "rt","maxo")], file=outname, quote=FALSE)
}


