# BaGFoot :
#
# Programmed by Songjoon Baek 2017
#
# Last Update : 4.18.2018

#library("hash");    # 2.2.6
#library("parallel");
#library("digest");   #0.6.9
#library("data.table"); # 1.10.0
#library("Cairo");      # 1.5-9
#library("aplpack");   # 1.3.0
#library('Rsamtools');  # 1.22.0
#library('BSgenome.Hsapiens.UCSC.hg19');  # 1.4.0
#library('BSgenome.Mmusculus.UCSC.mm9');  # 1.4.0

#library(doMC);
#library(foreach);
###################################bfoot.configuration #########################################

#MCCORES=getOption("mc.cores", 2L);	         # Number of CPU cores for parallel computation
MCCORES=8;
MCCORES_MOTIF = 8;     # How many motifs are calculated at the same time
#MCCORES_MOTIF = 1;     # How many motifs are calculated at the same time
FAST_AGGREGATION_MC_CORES = 8; # Number of CPU cores for less-memory dependent calculation

MAX_MOTIF_SITES = 10000000;    # Maximum number of FIMO sites for each motif
BFOOT_DEBUG = F;

#dyn.load("/opt/R/bagfootr/bagfoot/src/bagfoot.so");

#doMC::registerDoMC(cores=MCCORES);

###################################bfoot#########################################################

bagfootVersion <<- "0.9.7.02";
cat(sprintf('BaGFoot: v.%s\n', bagfootVersion));

printf <- function(...) { cat(sprintf(...)); };


extractBaseName<- function(filename) {
	bname = basename(filename);
	ext=regexpr("\\.",bname)[[1]];
	bnameNoExt = substring(bname, 1,ext-1);
	bnameNoExt;
}

countCutcountOnSites <- function(cutcount, sites) {
	if (sites=="" || is.na(sites) || is.null(sites)) {
		return (sum(sapply(cutcount, function(x) { sum(x$value);})));		
	}

	if (is.character(sites)) {
		hs <- readcsv(sites);
	} else {
		hs <- sites;
	}
  
	chroms = levels(hs$chr);
	sum_for_chr<- function(chr) {
#			print(chr);
			tdchr = cutcount[[chr]];
			hschr = hs[hs$chr==chr,];
			cc = vector("numeric", length=max(hs$ed+1));
			cc[tdchr$st] = tdchr$val;
			cc[tdchr$ed+1] = cc[tdchr$ed+1] - tdchr$val;
			cc = cumsum(cc);
			s=sum(sapply(seq_along(hschr$st),function(ii) sum(cc[hschr$st[ii]:hschr$ed[ii]])));
			rm(cc);
			s;
	}
	#browser();
	sum(unlist(parallel::mclapply(chroms, sum_for_chr,mc.cores=MCCORES)));
}

readCutCountOnSites<-function(cutcount_file="", site_file="") {

	cutcount = NULL;
	sites = NULL;

	reduced_bgr_filename = file.path(dirname(cutcount_file),sprintf("%s.bgr", paste(extractBaseName(cutcount_file),'on',extractBaseName(site_file), sep="_")));

	reduced_bgr_filename_zip = paste(reduced_bgr_filename,".gz", sep="");

	overlapping<-function(A,B)
	{
	  if (as.character(A$chr[1]) != as.character(B$chr[1])) {
		null;
	  } else {
		  R= vector("logical", length=nrow(A));

		  rangeA = 1:length(A$st);
		  rangeB = 1:length(B$st);

		  if ((length(rangeA) !=0) && (length(rangeB) !=0))
		  {
		  	  P=list(Ast=A$st[rangeA],Aed=A$ed[rangeA],Bst= B$st[rangeB],Bed= B$ed[rangeB]);
			  RR=  intersect_intervals_really_fast_parallel(P);
			  I1 = unique(RR);
		  	     R[rangeA[I1]]=TRUE;
		  }
		 # print(c(i,cend, cendcomp));
		  R;
	  }
	}

	filterCutcountByHotspot<-function(sites, chr) {
		cat(sprintf("reducing cutcount data for %s\n", chr));
		hotspot = sites[sites$chr==chr,c("chr","st","ed")];
		extendHotspot<-function(hotspot, extensionbp = 200)  {
			hotspot$st = hotspot$st - extensionbp;
			hotspot$ed = hotspot$ed + extensionbp;
			hotspot = hotspot[hotspot$st > 0,];
			hotspot;
		}
		if (nrow(hotspot)==0) {
			subcutcount=NULL;
		} else {
		#	browser();
			exh = extendHotspot(hotspot, extensionbp = 200);
			cutcountchr = cutcount[[chr]];

			ovlp = overlapping(cutcountchr, exh);
			subcutcount = cutcountchr[ovlp,];
		}
		subcutcount;
	};

	if (file.exists(reduced_bgr_filename_zip)) {
		cat(sprintf("reading the cutcount data file: %s\n", reduced_bgr_filename_zip));
		cutcount = readcutcount2(reduced_bgr_filename_zip);
	} else if (site_file == "" || is.na(site_file) || is.null(site_file)) {
		cat(sprintf('The site file is set to %s\n.', site_file));
		cutcount = readcutcount2(cutcount_file);
	} else {
		cat(sprintf("converting the cutcount data file: %s into %s\n", cutcount_file, reduced_bgr_filename));
		cutcount = readcutcount2(cutcount_file);
		sites = readhotspot(site_file);

	#	browser();
		reduced_cutcount= parallel::mclapply(names(cutcount), function(x) { filterCutcountByHotspot(sites, x); }, mc.cores=MCCORES);
		names(reduced_cutcount)<- names(cutcount);

		for (ii in seq_along(reduced_cutcount)) {
			chr = names(reduced_cutcount)[[ii]];
			dat = reduced_cutcount[[ii]];
			if (!is.null(dat)) {
				cat(sprintf("saving cutcount data (%s) for %s\n",reduced_bgr_filename, chr));
				dat$st = dat$st - 1;     # zero-based format
				if (ii==1) {
					write.table(dat, file=reduced_bgr_filename, sep=" ",   append=F, quote=F, row.names=F, col.names=F);
				} else {
					write.table(dat, file=reduced_bgr_filename, sep=" ",   append=T, quote=F, row.names=F, col.names=F);
				}
			}
		}
		#browser();
		system(sprintf('gzip %s', reduced_bgr_filename));  #gzip
		cutcount = reduced_cutcount;
	}
	cutcount;
}


readcutcount2 <- function(cutcountfile) {
	cutcount = readbgr(cutcountfile);
	rl = rle(as.vector(cutcount$chr));
	chrs = rl$values;
	ed=cumsum(rl$lengths);
	st= c(1,ed[-length(ed)]+1);
	cutcountPerChrom=lapply(1:length(chrs), function(ii) {cutcount[st[ii]:ed[ii],];});
	names(cutcountPerChrom) = chrs;
	cutcount={};
	cutcountPerChrom;
}


getOutputDirAndFilename<- function(outputdir_base, dataname,  sitefilename, halfwidth) {
	#dataname = datanamelist[[ii]];
	hotspotname = gsub(".csv","",basename(sitefilename));
	od = file.path(outputdir_base, sprintf("%s_on%s_%dbp", dataname,hotspotname, halfwidth));
	list(dir = od,
		 filename = file.path(od, sprintf('../%s_on%s_%dbp_out.txt', dataname, hotspotname,halfwidth)));
}


drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons <- function(
	sitefiledir='',
	cache='',
	outputdir='',
	datanamelist = list(),
	cutcountdatalist = list(),
	normalization="default",
	ncutcounts = NA,
	sitefilelist = list(),
	motiflistlist = list(),
	ftablefile = Fed_mm9_withMap,
	nuccodefilePattern = nuccode_hexamer_withMappability_35mer_mm9,
	halfwidthlist = c(50,100),
	comparisonIndices = list(c(1,2), c(3,2)),
	range= NA,
	yrange=c(-2.2,0.5),
	graphout=F) {

	threshold = 0;
	#browser();
	outputdir_base = outputdir;

 	if (!exists("ftablefile")) { # read once,
 	  freqtable<-readFrequencyTable(ftablefile);np = round(log(nrow(freqtable)-1)/log(4));
 	} else {
	  np = 6;   # by default
	}

  	chrs = paste("chr", c(1:22, "X","Y"),sep='');
	nuccodefiles= sapply(chrs,function(x) gsub("{chr}",x,nuccodefilePattern,fixed=T));
	doexist=file.exists(nuccodefiles);
	chrs= chrs[doexist];
	nuccodefiles=nuccodefiles[doexist];
	nuccodes= parallel::mclapply(nuccodefiles, function(nuccodefile) {
			nuccode=readNucleotideCodeForChromosomeForCuts(nuccodefile,np);
			nuccode; }, mc.cores=MCCORES,mc.preschedule=F);

	if (normalization=="default") {
		oldncutcounts= ncutcounts;
		#browser();
		ncutcounts = sapply(1:length(cutcountdatalist), function(kk) countCutcountOnSites(cutcountdatalist[[kk]],sitefilelist[[kk]]));
		
	}

	if (is.na(ncutcounts[1])) {
		stop("Error in drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons: ncutcounts should not be zero.")
	}

	for (halfwidth in halfwidthlist) {
		ll <- lapply(seq_along(datanamelist), function(ii) {
#				print(ii);
				dataname = datanamelist[[ii]];
				out=getOutputDirAndFilename(outputdir_base, dataname, sitefilelist[[ii]], halfwidth);
						run<-drawMotifAggPlotOnMotifSets(dataname = dataname,
								cutcount = cutcountdatalist[[ii]],
								ncutcount = ncutcounts[ii],
								hotspotfile = sitefilelist[[ii]],
								motiflist = read.table(motiflistlist[[ii]],stringsAsFactors=F) ,
								freqtablefile=ftablefile,
								nuccodes = nuccodes,
								sitefiledir=sitefiledir,
								cache=cache,
								outputdir=out$dir,
								halfwidth=halfwidth,
								logyrange=yrange,
								range=range,
								graphout=graphout);
						gc();
		}); #}, mc.cores=3);  # parallel::mclapply
	}  # for (half

	for (index in comparisonIndices) {
		for (halfwidth in halfwidthlist) {

			i1 = index[1]; i2 = index[2];

			dataname1 = datanamelist[[i1]];
			dataname2 = datanamelist[[i2]];
			motiflistfile1 =  motiflistlist[[i1]];
			motiflistfile2 =  motiflistlist[[i2]];
			out1 =  getOutputDirAndFilename(outputdir_base,dataname1 , sitefilelist[[i1]], halfwidth);
			out2 =  getOutputDirAndFilename(outputdir_base,dataname2 , sitefilelist[[i2]], halfwidth);

			title = sprintf('Footprinting Aggregation Plot (%s vs. %s)', dataname1, dataname2);

			datadir1 =  sprintf("%s_%dbp",dataname1, halfwidth);
			datadir2 =  sprintf("%s_%dbp",dataname2, halfwidth);
			MAoutputdir = file.path(outputdir, 'comparison');

			dlog(sprintf("\n\n\n\n\n\\nMAoutputdir=%s\n",MAoutputdir))
			#browser();

			MultipleAggregationPlot(
				dataname1=dataname1 ,
				dataname2=dataname2,
				motiflistfile1=motiflistfile1,
				motiflistfile2=motiflistfile2,
				outputfile1=out1$filename,
				outputfile2=out2$filename,
				title = title, threshold=threshold,
				datadir1=datadir1 ,
				datadir2=datadir2,
				outputdir =MAoutputdir,
#					scalefactor = c(scalefactor[i1], scalefactor[i2]),
				range=range,
				logyrange=yrange, graphout=graphout
			);
		}
	}

}

scatterDataOutGroups<-function(
	outputdir='',
	datanamelist = list(),
	motiflistlist = list(),
	sitefilelist = list(),
	halfwidthlist = c(50,100),threshold = 0,
	comparisonIndices = list(c(1,2), c(3,2)),
	range=NA) {

    out=c();
	for (index in comparisonIndices) {
		for (halfwidth in halfwidthlist) {

			i1 = index[1]; i2 = index[2];

			dataname1 = datanamelist[[i1]];
			dataname2 = datanamelist[[i2]];

			sitefilename = extractBaseName(basename(sitefilelist[[i1]]));


				if (length(motiflistlist)==1) {
				motiflistfile1 =  motiflistlist
			} else {
				motiflistfile1 = motiflistlist[[i1]]
			}

			out1=getOutputDirAndFilename(outputdir, datanamelist[[i1]], sitefilelist[[i1]], halfwidth);
			out2=getOutputDirAndFilename(outputdir, datanamelist[[i2]], sitefilelist[[i2]], halfwidth);
			outputfile1=out1$filename;
			outputfile2=out2$filename;

			title = sprintf('Footprinting Aggregation Plot (%s vs. %s)', dataname1, dataname2);

			datadir1 =  sprintf("%s_%dbp",dataname1, halfwidth);
			datadir2 =  sprintf("%s_%dbp",dataname2, halfwidth);
			MAoutputdir = file.path(outputdir, 'comparison');

			outputcsvfile = sprintf("%s_%s_On_%s_footprint_depth_table.csv", dataname1, dataname2, sitefilename);

			outfilename=ScatterDataOut(dataname1=dataname1,
						 dataname2=dataname2,
						 motiflistfile=motiflistfile1,
					outputfile1=outputfile1,
					outputfile2=outputfile2,
					title = title,threshold=threshold,datadir1=datadir1,datadir2=datadir2, outputcsv=
					outputcsvfile,
					range=range	);
			out=c(out, outfilename);
		}
	}
	out;
}

ScatterDataOut <- function(dataname1="", dataname2="", motiflistfile="",
					outputfile1="",outputfile2="",title = "",threshold=0,datadir1='',datadir2='', outputcsv=' ',range=NA) {

	motiflist=read.table(motiflistfile);

	calcObsMeans <- function (mobs1, mobs2,motifcenter, motiflogo, meanlog1, meanlog2, numsite) {
		motifoffset=2;
		halfwidth= (length(mobs1)-1)/2;
		motifwidth = nchar(motiflogo);
		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		hgt = 0;
		x1 = min(xr); x2 = max(xr);
		linea = x1 - motifoffset;
		lineb = x2 + motifoffset;
		motifregion = linea:lineb  + halfwidth + 1;
		flanking = setdiff(1:length(mobs1), motifregion);
		c1=log(sum(mobs2[motifregion])/sum(mobs1[motifregion]),2);
		c2=log(sum(mobs2[flanking])/sum(mobs1[flanking]),2);
		c3=log(sum(mobs2)/sum(mobs1),2);
		meancc1 = mean(mobs1);
	 	meancc2 = mean(mobs2);
		data.frame(lrmotif=c1, lrflanking=c2, lrtotal=c3, fdepth1= meanlog1, fdepth2= meanlog2,meancc1= meancc1, meancc2=meancc2, numsite);
	}

	calcprog<-function(ii) {
		motif = motiflist[ii,];
		name = as.character(motif$motif);
		#print(name);
		#rank = difforder[ii];

		sites1 = output1[ii,];
		sites2 = output2[ii,];

		numsite1 = sites1$numSites;
		numsite2 = sites2$numSites;

		motif = motiflist[ii,];
		name = as.character(motif$motif);
		motiflogo = as.character(motif$logo);
		motifcenter =motif$center;

		load(as.character(output1$graphdatafile[ii]));
		mobs1 = mobs;
		mexp1 = mexp;
		logratio1 =  logratio - baselogratio;
		meanlog1 = motiflogratio - baselogratio;

		load(as.character(output2$graphdatafile[ii]));
		mobs2 = mobs;
		mexp2 = mexp;
		logratio2 =  logratio - baselogratio;
		meanlog2 = motiflogratio - baselogratio;

		plot1yrange=0;

		halfwidth= (length(logratio1)-1)/2;
		filenamebase = make.names(name);


		calcObsMeans(mobs1, mobs2,motifcenter,motiflogo, meanlog1, meanlog2, numsite1);
	}

	if (!is.na(range[1])) {
		motiflist <- motiflist[range,];
	}

	output1 = read.table(outputfile1,sep='\t');  #200
	output2 = read.table(outputfile2,sep='\t');  #199

	selection= (output1$numSites>=threshold) & (output2$numSites>=threshold);

	output1 = output1[selection,];
	output2 = output2[selection,];
	motiflist = motiflist[selection,];
	cal<-do.call("rbind", lapply(1:nrow(motiflist), calcprog));
	tabl = data.frame(motiflist, cal);
	write.table(tabl, file=outputcsv)
	outputcsv;
}


runGfoot<-function(Gfoot, range=NA, graphout=F, halfwidths=c(50),normalization="default", yrange=c(-2.2,0.5)) {

	drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons(
		sitefiledir=Gfoot@motifDB@directory,
		cache= Gfoot@cachedir,
		outputdir=Gfoot@outputdir,
		datanamelist = list(Gfoot@control@name, Gfoot@treatment@name),
		cutcountdatalist = list(Gfoot@control@cutcount, Gfoot@treatment@cutcount),
		ncutcounts = c(Gfoot@control@count,Gfoot@treatment@count),
		normalization = normalization,
		sitefilelist = list(Gfoot@sitefile, Gfoot@sitefile),
		motiflistlist = list(Gfoot@motifDB@motiflistfile,Gfoot@motifDB@motiflistfile),
		ftablefile = Gfoot@gfootoption@biasfile,
		nuccodefilePattern = Gfoot@gfootoption@hexamer_pattern,
		halfwidthlist = halfwidths,
		comparisonIndices = list(c(1,2)),
		range=range,
		graphout=graphout,
		yrange=yrange
	);

	scatterDataOutGroups(
		outputdir=Gfoot@outputdir,
		datanamelist = list(Gfoot@control@name, Gfoot@treatment@name),
		motiflistlist =  list(Gfoot@motifDB@motiflistfile,Gfoot@motifDB@motiflistfile),
		sitefilelist = list(Gfoot@sitefile, Gfoot@sitefile),
		halfwidthlist =halfwidths,
		comparisonIndices = list(c(1,2)),
		range=range
	);
}

####
# OOP definition of

#' @export CutCount
CutCount <- setClass("CutCount",
					representation(file = "character",
							  cutcount = "list",
					   	      count = "numeric",
							  name = "character"),
					prototype = list(file="",
									 cutcount = list(),
									 count = 0,
									 name=""),
					validity = function(object) {
						if (file.exists(object@file) && object@count > 0) TRUE else "File not exist or count==0";
					});

#' @export
setGeneric(name="readCutCountSites",   #CutCount::readCutCountOnSites
			def = function(obj, site_filename) {
				standardGeneric("readCutCountSites");
			});
#' @export
setMethod(f="readCutCountSites",   #CutCount::readCutCountOnSites
		  signature="CutCount",
		  definition=function(obj,site_filename) {
			  if (length(obj@cutcount)==0) {
				  cat(sprintf("reading %s ...\n", obj@file));
				  obj@cutcount <- readCutCountOnSites(cutcount_file = obj@file,site_file = site_filename);
			  }
			  return(obj);
		});

#' @export MotifDB
MotifDB <- setClass("MotifDB",
					 representation(
							   motiflistfile = "character",
							   directory = "character"
							   ),
					  prototype = list(
						motiflistfile ='',
						directory='')
					);

#' @export GFootOption
GFootOption <- setClass("GFootOption",
					 representation(
							   biasfile = "character",
							   hexamer_pattern = "character"),
					  prototype = list(
						biasfile="",
						hexamer_pattern='')
					);

#' @export GFoot
GFoot <- setClass("GFoot",
					 representation(control = "CutCount",
							   treatment = "CutCount",
							   sitefile = "character",
							   motifDB = "MotifDB",
							   gfootoption = "GFootOption",
							   outputdir = "character",
							   cachedir = "character",
							   outputfiles="character"),
					  prototype = list(
						sitefile='',
						motifDB = MotifDB(),
						gfootoption= GFootOption(),
						outputdir = getwd(),
						cachedir = file.path(getwd(), "cache"),
						outputfiles='')
					);
#' @export
setGeneric(name="loadCutcount",    #GFoot::run
			def = function(obj) {
				standardGeneric("loadCutcount");
 			});
#' @export
setMethod(f="loadCutcount",		#GFoot::calcHeatmap
		  signature="GFoot",
		  definition=function(obj) {
				if (length(obj@control@cutcount)==0) {
					obj@control<-readCutCountSites(obj@control, obj@sitefile);   #read cutcount
				}
				if (length(obj@treatment@cutcount)==0) {
					obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
				}
				return(obj);
			});

#' @export
setGeneric(name="run",    #GFoot::run
			def = function(obj,range=NA,graphout=F,yrange=c(-2.2,0.5), mc.cores=MCCORES) {
				standardGeneric("run");
 			});

#' @export
setMethod(f="run",		#GFoot::run
		  signature="GFoot",
		  definition=function(obj, range=NA, graphout=F, yrange=c(-2.2,0.5), mc.cores= MCCORES) {
				#obj <- loadCutcount(obj);

				if (length(obj@control@cutcount)==0) {
				#	browser();
					obj@control<-readCutCountSites(obj@control, obj@sitefile);   #read cutcount
				}
				if (length(obj@treatment@cutcount)==0) {
					obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
				}
				#obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
				cat(sprintf("Num. of CPU cores to use: %d\n",mc.cores));
				MCCORES = mc.cores;
				obj@outputfiles=runGfoot(
					obj,              # GFoot obj
					range =  range,   # Range of CREB
					graphout = graphout,
					yrange= yrange
				);
				return(obj);
			});


######################################bfoot_calc##################################################



gfoot_ptm <<- proc.time();  # initial time

dlog <- function(output, ...) {
#~ #	cat(sprintf("[%g sec] %s", etime[["elapsed"]], output, ...));
}

#' @useDynLib bagfoot intersect_intervals
intersect_intervals_fast_index2C <- function(X1,X2,Y1,Y2)
{
	Isorted1= order(X1);
	Isorted2=  order(Y1);
	.Call("intersect_intervals", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}

#' @export
readFrequencyTable<-function(ftablefile) {
	cat(sprintf("reading %s ......\n", ftablefile));
	ftable = read.table(ftablefile);
	nr = nrow(ftable);
	ftable[nr,c("GPRatio","ObCutRatio")] = 0;
	ftable = cbind(ftable, Aj = ftable$ObCuts/ ftable$GenomicPositionCount);
	ftable$Aj[nr] = 0;
	# cat(sprintf("Num of cuts =  %d ......\n", sum(ftable$ObCuts)));
	ftable;
}


intersect_intervals_fast_parallel_index2 <- function(P) {

    R={};
    if (length(P$Ast) ==0  || length(P$Bst)==0) {
      return(R);
    }
    R=intersect_intervals_fast_index2C(P$Ast, P$Aed, P$Bst,P$Bed);
    R;
}

OverlappingIndexAParallel <- function(A,B,mc.cores=MCCORES)
{

  OverlappingIndexAParallelPerChrom<-function(i) {
   #   cat(paste(chrs[i],'\n'));
	  result = NULL;
      iA = 1:length(A$st);
      iB = 1:length(B$st);
      rangeA= iA[A$chr == chrs[i]];
      rangeB= iB[B$chr == chrs[i]];
      if ((length(rangeA) !=0) && (length(rangeB) !=0))
      {
      	  P=list(Ast=A$st[rangeA],Aed=A$ed[rangeA],Bst= B$st[rangeB],Bed= B$ed[rangeB]);
	      RR=  intersect_intervals_fast_parallel_index2(P);
	      if (!is.null(RR)) {
			result=rangeA[RR];
		  }
      }
	  result;
  }
  chrs=levels(factor(A$chr));

  if (mc.cores==1) {
	  kk=lapply(seq_along(chrs),OverlappingIndexAParallelPerChrom);
  } else {
	  kk=parallel::mclapply(seq_along(chrs),OverlappingIndexAParallelPerChrom , mc.cores=mc.cores);
  }
 # browser();
  unlist(kk);
}

CompareHotspotsFaster<-function(A,B, mc.cores=MCCORES) {
	idx = OverlappingIndexAParallel(A,B, mc.cores=mc.cores);
	A[idx,];
}

filesizenotzero <- function(filename) {
	return(ifelse(file.exists(filename),file.info(filename)$size > 0,FALSE ));
}

readbgr2 <- function(filename)
{
    if (filename=="")
    {
        R=NULL;
        R;
    } else {
        # browser();
        R= readLines(filename, n=7);
        i=1;
        while (i < 7)
        {
            if (substr(R[[i]],1,3)=="chr") break;
            i=i+1;
        }
        if (i >= 7)
        {
            R={};
            sprintf("%s is not a readable bed file.", A);
            return;
        } else
        {
            #R=read.csv(filename, sep="", header=FALSE,comment.char="", skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE);
            if (length(grep(".gz", filename))==1) {
				freadfilename = paste('zcat', filename);
			} else {
				freadfilename = filename;
			}

            R=data.frame(data.table::fread(freadfilename,  header=FALSE, skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE));
            names(R) <- c("chr","st","ed","value");
            R$chr <- factor(R$chr);
            R$st = R$st+1;
            R;
        }
    }
}


readcutcount3 <- function(cutcountfile) {
	cutcount = readbgr2(cutcountfile);
	cutcountPerChrom = split(cutcount, cutcount$chr);
	cutcount={};
	cutcountPerChrom;
}

readBGRAsList <- readcutcount3

drawMotifAggPlotOnMotifSets<- function(dataname = '',
	cutcount = '',
	ncutcount= 0,
	hotspotfile = '',
	motiflist = '',
	freqtablefile=NA,
	nuccodes = NA,
	sitefiledir='',
	cache='.',
	outputdir = '',
	halfwidth=50,
	strandwisePlot=F,
	logyrange=c(-2.2,0.5),
	range=NA,
	graphout=F
) {

	ftable = readFrequencyTable(freqtablefile);

	if (!file.exists(cache) && cache != '') {
		dir.create(cache, recursive=T);
	}

	if (cache=='') {
		cache = getwd();
	}

	if (!file.exists(outputdir) && cache != '') {
		dir.create(outputdir, recursive=T);
	}


	if (outputdir=='') {
		outputdir = getwd();
	}

	begintime = Sys.time();
#	out<- lapply(1:length(motiflist$motif), function(ii) {
	if (is.na(range[1])) {
		range=1:length(motiflist$motif);
	}

	jj=0;

	result = vector("list", length=length(motiflist$motif));

	if (ncutcount == 0) {
		stop("Error in drawMotifAggPlotOnMotifSets: ncutcount should not be zero.")
	}

	scalefactor = 100000000/ncutcount;
	cat(sprintf('num cutcount = %g,  scale factor = %g\n', ncutcount, scalefactor));


	runCalc<-function(ii, graphout=F) {
			#if (ii==255) {
			#	browser();
			#}
			outputfilename = make.names(sprintf('%s_%s_hw%dbp',motiflist$motif[ii], dataname,halfwidth));
			jpgfile1= file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts.jpg",outputfilename));


				cat(sprintf('No.%3d Motif:%s (%s %d bp)\n', ii, motiflist$motif[ii], motiflist$logo[ii], motiflist$motiflength[ii]));
				#cat(sprintf('Center=%d', ii, motiflist$motif[ii], motiflist$logo[ii], motiflist$motiflength[ii]));
				motif=list(
				name = as.character(motiflist$motif[ii]),
				sitefile=file.path(sitefiledir, motiflist$filename[ii]),
				motiflogo=as.character(motiflist$logo[ii]),
			#	logopic=NA,
				motifcenter=motiflist$center[ii]);

				motifname= motiflist$motif[ii];
				HS=list(
					 name='Hotspot',
					 hotspotfile=hotspotfile,
					 cutcountfile=cutcount,
					 cutcountsignature=digest::digest(cutcount));


				ChIP=list(   # Rep 1: GH1158 Rep 2: GH1159
				 name='ChIP',
				 hotspotfile=hotspotfile);


				AggregationPlotData= list(motifinfo=motif,
				  dhsinfo =HS,
				  chipinfo = ChIP,
				  ftablefile= freqtablefile,
				  ftable = ftable,
				  nuccodes =nuccodes);

				ss=DrawStamMotifAggPlotOnChIPBoundRegions(AggregationPlotData,
					output=outputfilename,
					halfwidth=halfwidth,
					title=sprintf('%s (%s)',dataname, motifname),
					yrange=logyrange,
					cache=cache,
					outputdir=outputdir,
					strandwisePlot=strandwisePlot,
					scalefactor = scalefactor,
					graphout=graphout);
				#	cat("warning: bfoot.R line 813:\n");
				ss;
				

		}

	if (MCCORES_MOTIF == 1) {
		result<-lapply(range, function(x) { runCalc(x, graphout=graphout); }); # parallel::mclapply
	} else {
		
		#result<- foreach::foreach(x = range, .combine=c) %dopar% runCalc(x, graphout=graphout);	
		result<-parallel::mclapply(range, function(x) { runCalc(x, graphout=graphout); },mc.cores=MCCORES_MOTIF );
		# ); # parallel::mclapply

	}
    #cat("warning: bfoot.R line 825:\n");
	rout<-do.call("rbind", result);
	hotspotname = gsub(".csv","",basename(hotspotfile));
	outputfilename = file.path(outputdir, sprintf('../%s_on%s_%dbp_out.txt', dataname,hotspotname,  halfwidth));
#	browser();
	write.table(rout,file=outputfilename,sep='\t');
	cat("warning: bfoot.R line 831:\n");
	#}
	rout;
}


quicksummaryListOrString<-function(obj) {
	if (is.character(obj)) {
		r=digest::digest(obj);
	} else {
		r=digest::digest(c(unlist(lapply(obj, nrow))));
	}
	r;
}


DrawStamMotifAggPlotOnChIPBoundRegions <- function(S, tag="", motifdir=c(1,-1),output='',  halfwidth=20,movingaverage=F, plot1yrange=0,  yrange=c(-2.5,1.5),title='', cache='.', outputdir='.',strandwisePlot=F, scalefactor=1, graphout=F) {



	S$cutcountname = S$dhsinfo$name;
	S$DHS = S$dhsinfo$hotspotfile;
	S$cutcount = S$dhsinfo$cutcountfile;
	S$cutcountsignature = S$dhsinfo$cutcountsignature
	S$CHIP = S$chipinfo$hotspotfile;
	S$chipname = S$chipinfo$name;

	S$motif=  S$motifinfo$sitefile;
	S$output = paste(S$output,tag, sep="");

	#dlog("DrawStamMotifAggPlotOnChIPBoundRegions():\n");
	#print(S$motif);

	#print(yrange);


	savefile = file.path(cache,sprintf('%s_%s_motifplotresult.dat',make.names(S$motifinfo$name), digest::digest(c(S$motif,S$CHIP,S$DHS,S$cutcountsignature,halfwidth))));

	readHotspotRegions <- function(filename) {
				dhs  =data.frame(data.table::fread(filename));
				names(dhs) <- c("ID","chr","st","ed", "MaxD","AveD","Zscore","pvalue");
				dhs = dhs[dhs$chr != "chrM",];   # eliminate chrM regions from the analysis
				dhs;
	}

	
	if (!file.exists( savefile) || BFOOT_DEBUG) {

		motif_chip_dhs=NULL;
			motif=data.frame(data.table::fread(S$motif))[,-1];
			motif = motif[,1:5];
			istart = match('start',names(motif));
			if (istart == 2) {
				names(motif)[(istart-1):(istart+3)] = c("chr","st","ed", "dir","value");
				direction = motifdir;
				direction[direction==1] = '+';
				direction[direction==-1] = '-';
				motif = motif[is.element(motif$dir, direction),];
			} else {
				names(motif)=c("chr","st","ed", "value","dir");
				motif = motif[is.element(motif$dir, motifdir),];
			}

			motifwidth = motif$ed[1]-motif$st[1]+1;
			result=NULL;
			if (is.na(S$DHS) || is.null(S$DHS) || (S$DHS=="")) {
				motif_chip_dhs = motif;									
			} else {
				dhs = readHotspotRegions(S$DHS);
				if (S$DHS == S$CHIP) {
					chip = dhs;
				} else {
					chip =data.frame(data.table::fread(S$CHIP));
					if (ncol(chip)>7) {
						names(chip)=  c("ID","chr","st","ed", "MaxD","AveD","Zscore","pvalue");
					} else if (ncol(chip)==3) {
						names(chip)=  c("chr","st","ed");
					}
				}

				if (nrow(motif)>0) {
				#				dlog(">>>>>\n");
					motif_chip= CompareHotspotsFaster(motif, chip, mc.cores=MCCORES);
				#				dlog("CompareHotspots>>>>>\n");
					if (S$DHS == S$CHIP) {
						motif_chip_dhs = motif_chip;
					} else {
						motif_chip_dhs= CompareHotspotsFaster(motif_chip, dhs, mc.cores=MCCORES);
					}
				}

				if (nrow(motif_chip_dhs)> MAX_MOTIF_SITES) {
					set.seed(1);
					motif_chip_dhs = motif_chip_dhs[sample(nrow(motif_chip_dhs), MAX_MOTIF_SITES),];
				}
			}
    #browser();
#		if (nrow(motif_chip_dhs) < 1) {
#			browser();
#		} 
		result=StamMotifAggPlot(motif_chip_dhs, S,output=output, halfwidth=halfwidth,movingaverage=movingaverage,yrange=yrange,plot1yrange=plot1yrange, title=title,cache=cache,outputdir=outputdir,strandwisePlot=strandwisePlot,scalefactor=scalefactor, graphout=graphout);
		save(result, file=savefile);
		#load('/home/baeks/Data/project/saoriATAC/analysis_LTED2/bagfoot_by_clusters/T47D_8clusters/cl1/CACHE_T47D_Control/ZNF740_10bp_6369db133f60569f7f1c39ad4005bec7_motifplotresult.dat');
		#cat(sprintf("#sites = %d\n", nrow(motif_chip_dhs)));
	} else {
		cat(sprintf('loading %s\n', savefile));
		load(savefile);
	}
	#cat('line 941\n');
	result;
}   ###### DrawStamMotifAggPlotOnChIPBoundRegions


calcNormalizedExpectedObservedCounts <- function (sites, S, cache, basehalfwidth, halfwidth, scalefactor, motifbuffer=2,halfflankingregionwidth=50, heatmapout='') {

	getSubCount <- function(from, to) {
		subrange_st = motifcenter+((basehalfwidth)+from);
		subrange_ed = motifcenter+((basehalfwidth)+to);

		isForward = counts$sites$dir==1 | counts$sites$dir=='+';
		st = ifelse(isForward ,counts$region$st + subrange_st - 1,counts$region$ed - subrange_ed + 1)
		ed = ifelse(isForward,counts$region$st + subrange_ed - 1,counts$region$ed - subrange_st + 1)

		subobserved = counts$observed[,(subrange_st:subrange_ed)-1];

		#browser();
		if (length(isForward)==1) {
		  if (isForward) {
		    subobserved = counts$observed[,(subrange_st:subrange_ed)+1];
		  }
		} else {
		  subobserved[isForward,] = counts$observed[isForward,(subrange_st:subrange_ed)+1];
		}
		subexpected = counts$expected[,subrange_st:subrange_ed];
		subcounts=list(sites=counts$sites, dir=counts$sites$dir, region =  data.frame(chr = counts$region$chr, st = st,  ed = ed), observed=subobserved*scalefactor, expected=subexpected*scalefactor);
	}

	average_aggregation<-function(heatmap){

		sum1 = apply(heatmap, 1,sum);
		notnull = !is.na(sum1);
		num= sum(notnull);
		apply(heatmap[notnull,], 2, sum)/num;
	}

	calc_logratio<-function(cc) {
	  if (is.null(dim(cc$expected)) && length(cc$expected)>0) {
      mexp = cc$expected;
	  } else {
	    rowNA1 = is.na(apply(cc$expected,1,sum));
	    mexp = average_aggregation(cc$expected[!rowNA1,]);
	  }
	  if (is.null(dim(cc$observed)) && length(cc$observed)>0) {
      mobs = cc$observed;
	  } else {
	    rowNA2 = is.na(apply(cc$observed,1,sum));
	    mobs = average_aggregation(cc$observed[!rowNA2,]);
	  }
		mexp[mexp<0.0001] = 0.0001;
		mobs[mobs<0.0001] = 0.0001;
		logratio = log(mobs/mexp,2);
		list(mobs=mobs, mexp=mexp,logratio= logratio);
	}
	

	cutcountfile=S$cutcount;
	motifname= S$motifinfo$name;
	motifcenter=S$motifinfo$motifcenter;

	motiflogo= S$motifinfo$motiflogo;
	motifwidth=nchar(motiflogo);

	cutcountname = S$cutcountname;

	ftable= S$ftable;
	ftablefile=S$ftablefile;


	nuccodes= S$nuccodes;

	nsite = nrow(sites);

	if (motifcenter==0) {
			motifcenter = floor(nchar(motiflogo)/2)+1;
	}

	savefile = file.path(cache,sprintf('%s_%s_2.dat',make.names(S$motifinfo$name), digest::digest(c(sites,quicksummaryListOrString(cutcountfile) , S$dhsinfo$name, halfwidth, motifname,ftablefile))));

	if (!file.exists( savefile) || heatmapout != '' || BFOOT_DEBUG) {
		if (is.list(cutcountfile)) {
			cutcount = cutcountfile;
		} else {
				cutcount = readcutcount3(cutcountfile);
		}
	  Rfwd={};
	  Rbwd={};
    if (nsite > 0) {
      counts=CalcObservedExpectedCuts(S, sites, cutcount,motifcenter=1, halfspan = basehalfwidth, nuccodes=nuccodes,motifwidth=motifwidth, ftable=ftable);
      cutcount = {};
      #	browser();
      subcounts <- getSubCount(-halfwidth, halfwidth);
      FRcounts <- getSubCount(-motifbuffer-(halfflankingregionwidth-1)-motifcenter, motifwidth+motifbuffer+halfflankingregionwidth-motifcenter);
      flankingregionwidth= motifwidth+2*(motifbuffer+halfflankingregionwidth)
      flankingcolumn=c(1:halfflankingregionwidth, (flankingregionwidth-halfflankingregionwidth+1):flankingregionwidth);
      R2<- calc_logratio(FRcounts);
      baselogratio = mean(R2$logratio[flankingcolumn], trim=0.1);
      motiflogratio = mean(R2$logratio[-flankingcolumn], trim=0.1);
      R=calc_logratio(subcounts);
    } else {
      baselogratio = 0;
      motiflogratio = 0;
      R=list(mobs = rep(0, 2*basehalfwidth +1),
             mexp = rep(0, 2*basehalfwidth + 1),
             logratio = rep(0,2 * basehalfwidth+1));
    }
		result=list(Rboth=R, Rfwd=Rfwd, Rbwd=Rbwd, baselogratio=baselogratio, motiflogratio=motiflogratio);
		if (heatmapout!='') {
			heatmapfile1 = file.path(sprintf('%s.expected.csv',heatmapout));
			heatmapfile2 = file.path(sprintf('%s.observed.csv',heatmapout));
			cat(sprintf('saving %s\nsaving %s\n', heatmapfile1, heatmapfile2));
			write.table(subcounts$expected, file= heatmapfile1, sep=',', row.names=F, col.names=F);
			write.table(subcounts$observed, file= heatmapfile2, sep=',',  row.names=F, col.names=F);
		}
		save(result, file=savefile);
	} else {
		cat(sprintf('loading %s\n', savefile));
		load(savefile);
	}
	result;
}


drawAggregatedObservedExpectedCutcountPlot <- function (mobs, mexp, motifcenter, halfwidth, title,  motifwidth, linea, lineb,motiflogo,  plot1yrange) {

  colorcode=c('red','darkolivegreen','darkblue');

  miny=min(mobs,mexp);
  maxy=max(mobs,mexp);

  xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
  if (length(plot1yrange)<2) {
	yrange = c(miny,maxy);
  } else {
	yrange= plot1yrange;
  }


  hgt = 0;

  op <- par(mar=c(5,6,4,2) + 0.1);

  plot((-halfwidth:halfwidth), mobs, "n", main = title, col = colorcode[1], ylim=c(miny,maxy),
	   xlab = '(BP)',  ylab='', lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
	   xaxs='i',yaxs='i', las=1);
  mtext('Observed vs. Expected\nCounts', side=2, line=3, cex=1.5);


	text(xr, rep(maxy- (maxy-miny)*0.05, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

  #		if (nplot>1) {


  lines((-halfwidth:halfwidth), mobs, col= colorcode[1],lwd=2.5,lty=1);    #observed
  lines((-halfwidth:halfwidth), mexp, col= colorcode[2],lwd=2.5,lty=1);    #expected

  if (motifwidth >0) {
	lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
	lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);
  }

  legend(-halfwidth+1, y = yrange[2] - (yrange[2]-yrange[1])*0.05, bty='n',c("Observed", "Expected"), lwd=3, col = colorcode,cex=1.3);

}


drawLogratioPlot <- function (logratio, yrange, halfwidth, title, motiflogratio,baselogratio, linea, lineb, motiflogo, motifcenter) {

	miny=min(logratio);
	maxy=max(logratio);
	if (length(yrange)<2) {
	  yrange = c(miny,maxy);
	}
	op <- par(mar=c(5,6,4,2) + 0.1);

	linewidth= 3;
	plot((-halfwidth:halfwidth), logratio,"n", main = title,
		 xlab = '(BP)', ylab='', ylim=yrange, lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
		 xaxs='i',yaxs='i',las=1);

	# baseline + mean line
	lines(c(-halfwidth, halfwidth), c(motiflogratio,motiflogratio), col= 'darkred', lwd=2, lty=2);
	lines(c(-halfwidth, halfwidth), c(baselogratio,baselogratio), col= 'black', lwd=2, lty=3);
	#lines(c(-halfwidth, halfwidth), c(0,0), col= 'black', lwd=1.5, lty=2);
	#}

	lines((-halfwidth:halfwidth), logratio,"l", col = 'darkblue',lwd=linewidth);

	# y label
	mtext('Log Ratio',side=2, line=3,cex=1.5);
	par(op);

	lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
	lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);

	 xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;

	 text(xr, rep(maxy, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

}


StamMotifAggPlot<- function(sites, S,output='', halfwidth=20,movingaverage=F, yrange=0,plot1yrange=0, title='', motifoffset=2, trim=0, strandwisePlot=F,basehalfwidth=200, cache='.',outputdir='.',scalefactor=1,ttest=F, graphout=F) {


	motifcenter=S$motifinfo$motifcenter;
	motiflogo= S$motifinfo$motiflogo;
	motifwidth=nchar(motiflogo);


	if (!file.exists(outputdir) && outputdir != '') {
		dir.create(outputdir, recursive=T);
	}

	datafile = file.path(outputdir,sprintf("aggregation_%s.dat",output));
#
	#browser();;;

	if (!filesizenotzero(datafile) || BFOOT_DEBUG) {

		RR <- calcNormalizedExpectedObservedCounts(sites, S, cache, basehalfwidth, halfwidth, scalefactor)

		mobs = RR$Rboth$mobs;
		mexp = RR$Rboth$mexp;
		logratio = RR$Rboth$logratio;
		baselogratio=RR$baselogratio;
		motiflogratio=RR$motiflogratio;
		footprintdepth=motiflogratio-baselogratio;

		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		x1 = min(xr); x2 = max(xr);

		linea = x1 - motifoffset;
		lineb = x2 + motifoffset;

		halfpoint= floor(length(logratio)/2)+1;
		motifregion=seq(from=halfpoint+x1-motifoffset, to=halfpoint+x2+motifoffset);
		motifinside = rep(F, length(logratio));
		motifinside[seq(from=halfpoint+x1-motifoffset, to=halfpoint+x2+motifoffset)] = T;



		if (ttest) {
		res<-wilcox.test(logratio~motifinside,  alternative="greater");

		mu = 0;
		xbar = mean(logratio[motifinside],trim=0.1);
		nn = sum(motifinside);
		ssd= sd(logratio);
		ttt= (xbar-mu) / (ssd/sqrt(nn));

		alpha = 0.05;
		t.alpha = qt(1-alpha, df=nn-1);

		ttestoutputfile = file.path(outputdir,sprintf('%s_ttest.txt', output));
		sink(ttestoutputfile);
		cat(sprintf('t=%f\n', ttt));
		cat(sprintf('critical value at %f = %f\n', alpha, t.alpha));
		cat(sprintf('Univariate t-test pvalue = %f\n', pt(ttt, df=nn-1)));
	#	print(res);
		sink();
		}

		trimmean=mean(logratio[motifregion],trim=trim);
		save(mobs, mexp, logratio,baselogratio,motiflogratio, footprintdepth,linea,lineb, file=datafile);

	} else {
		cat(sprintf("loading datafile:%s\n", datafile));
		load(datafile);
	}
	if (graphout && nrow(sites)>0) {
		if (length(yrange)<2) {
			rangestring = "";
		} else {
			rangestring = paste(yrange, collapse="to");
		}
		jpgfile1= file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts%s.jpg",output,rangestring));
		jpgfile2= file.path(outputdir,sprintf("log_ratio_%s_observed_expected_cuts%s.jpg",output,rangestring));
		pdffilename1 =  file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts%s.pdf",output,rangestring));
		pdffilename2 =  file.path(outputdir,sprintf("log_ratio_%s_observed_expected_cuts%s.pdf",output,rangestring));

		jpeg(jpgfile1,width = 600, height = 600);

		drawAggregatedObservedExpectedCutcountPlot(mobs, mexp, motifcenter, halfwidth, title,  motifwidth, linea, lineb,motiflogo,  plot1yrange);


		dev.off();


		Cairo::CairoPDF(file = pdffilename1, width = 8, height = 8,   bg = "white");
		drawAggregatedObservedExpectedCutcountPlot(mobs, mexp, motifcenter, halfwidth, title,  motifwidth, linea, lineb,motiflogo,  plot1yrange);

		dev.off();

		jpeg(jpgfile2,width = 600, height = 600);
		drawLogratioPlot(logratio, yrange, halfwidth, title, motiflogratio,baselogratio,  linea, lineb, motiflogo, motifcenter);
		dev.off();


		Cairo::CairoPDF(file = pdffilename2, width = 8, height = 8,   bg = "white");
		drawLogratioPlot(logratio, yrange, halfwidth, title, motiflogratio,baselogratio,  linea, lineb, motiflogo, motifcenter);
		dev.off();
	}
	result = list(name= title, numSites=nrow(sites), meanlog=motiflogratio, baselog=baselogratio, footprintdepth=footprintdepth, graphdatafile=datafile);
	result;
}
#CalcObservedExpectedCuts(S, sites, cutcount,motifcenter=1, halfspan = basehalfwidth, nuccodes=nuccodes,motifwidth=motifwidth, ftable=ftable);
CalcObservedExpectedCuts <- function(S, motif, tdensity,motifcenter=0, halfspan = 50, nuccodes=NA,motifwidth=0, ftable=NA) {
	  if (motifcenter==0) {
		  midpoint= ifelse(motif$dir==1,motif$st + floor(motifwidth/2), motif$ed-floor(motifwidth/2));
	  } else {
		  midpoint= ifelse(motif$dir==1, motif$st + motifcenter -1, motif$ed-(motifcenter-1));
	  }
	  tregion = data.frame(chr= motif$chr, st = midpoint - halfspan + ifelse(motif$dir==-1,-1,0), ed = midpoint + halfspan+ ifelse(motif$dir==-1,-1,0));  #AP1

	#  browser();
	  # Remove out-of-range regions from sites;
	  chroms = names(tdensity);
	  isValid = vector("logical", length=nrow(tregion));
	#  browser();
	  for (chr in chroms) {
		maxcoord = max(tdensity[[chr]]$ed);
		isValid[tregion$chr==chr] = tregion$ed[tregion$chr==chr] <= maxcoord;
	  };

	 if (nrow(motif)==45) {
	#	  browser();
	  }	 
	  
	  motif = motif[isValid,];
	  tregion = tregion[isValid,];
	  
	  
	  ObservedCuts = CalcHeatmapFastMatrixOutputC(tregion, tdensity, dir=motif$dir);
	  
	#  browser();
	  observedsum= apply(ObservedCuts, 1, sum);

	  ExpectedRates = CalcExpectedRates(tregion, dir=motif$dir,nuccodes=nuccodes,ftable=ftable);

	  ExpectedCuts = matrix(0, nrow=nrow(ExpectedRates), ncol = ncol(ExpectedRates));

	  for (ii in 1:nrow(ExpectedCuts)) {
	  	ExpectedCuts[ii,] = observedsum[ii]* ExpectedRates[ii,];
	  }
	  #printf('bfoot_calc.R: CalcObservedExpectedCuts()'); browser();
	  list(sites=motif, region=tregion, observed=ObservedCuts, expected=ExpectedCuts);
}  #CalcAggregationMotif

readNucleotideCodeForChromosomeForCuts<-function(nuccodefile,np) {

	nuccode= readRDS(nuccodefile);
	nuccodeAtCuts = shiftarray(nuccode$code, np/2);
	nuccodeAtCuts[1:np]=4^np+1;
	nuccode={};

	as.integer(nuccodeAtCuts);
}

CalcExpectedRates<- function(tregion, dir=NA,nuccodes=NA,ftable=NA) {

  locchrom = as.character(unique(sort(tregion$chr)));
  numchrom = length(locchrom);
  Plist = {};
  nregion = nrow(tregion);
  range=seq(from=1, to=nregion);
  chrs=unique(as.character(tregion$chr));


  calcPerchrom<-function(ch) {
   #  print(sprintf('CalcExpectedRates: %s',ch));
      nuccode = nuccodes[[ch]];
      region = subset(tregion, chr==ch);
      #rangechr = subset(range, tregion$chr==ch);
      dirchr = subset(dir, tregion$chr==ch);
      numregion = nrow(region);
      R={};
      inRange= region$ed < length(nuccode);
      outRange = !inRange;
      regionInRange = region[inRange,];
      calcPj<-function(ii) {
        site= regionInRange[ii,];
        Aj_CleavageRate=ftable$Aj[nuccode[site$st:site$ed]];
         Pj=Aj_CleavageRate / sum(Aj_CleavageRate);
         if (dirchr[ii]==1) {
         	returnvalue=Pj;
         } else {
         	returnvalue = rev(Pj);
         }
         returnvalue;
      }
      if (MCCORES==1) {
		QQ<-lapply(1:nrow(regionInRange) ,calcPj);
	  } else {
		QQ<-parallel::mclapply(1:nrow(regionInRange) ,calcPj, mc.cores=MCCORES);
	  }
      Qj= do.call("rbind",QQ)
      QC =array(0,dim=c(length(inRange), ncol(Qj)));
	  QC[inRange,] = Qj;
      QC;
  }

  Pj = parallel::mclapply(locchrom, calcPerchrom, mc.cores=MCCORES);

  filterNANresult <- function(SS) {
	  if (is.numeric(SS)) {
		  nullrow = is.nan(apply(SS,1,sum));
		  SS[nullrow,] = array(0, dim=c(sum(nullrow), ncol(SS)));
		  result = SS;
	  } else {
		  result = NA;
	  }
	  result;
  }

  Pjj = lapply(Pj, filterNANresult);
  ExpectedCuts=do.call(rbind, Pjj);
  ExpectedCuts;
}

MultipleAggregationPlot <- function(dataname1="", dataname2="", motiflistfile1="", motiflistfile2="",
					outputfile1="",outputfile2="",title = "",threshold=0,datadir1='',datadir2='',
					outputdir=' ',scalefactor= c(1,1),range=NA, logyrange=0, graphout=F) {

 	if (!file.exists(outputdir) && outputdir != '') {
		dir.create(outputdir, recursive=T);
	}

	motiflist1=read.table(motiflistfile1,stringsAsFactors=F);
	motiflist2=read.table(motiflistfile2,stringsAsFactors=F);
	#browser();

	mfs1 = as.character(motiflist1$motif);
	mfs2 = as.character(motiflist2$motif);

	commonmotifs = intersect(mfs1,mfs2);

	common1= motiflist1$motif %in% commonmotifs;
	common2= motiflist2$motif %in% commonmotifs;

	motiflist1 = subset(motiflist1,common1);  # 200 -> 186
	motiflist2 = subset(motiflist2,common2);  # 199 -> 186

	printf("reading %s..\n", outputfile1);
	printf("reading %s..\n", outputfile2);

	output1 = read.table(outputfile1,sep='\t');  #200
	output2 = read.table(outputfile2,sep='\t');  #199

	output1 = output1[common1,]; #186
	output2 = output2[common2,]; #186

	if (is.na(range[1])) {
		selection= (output1$numSites>=threshold) | (output2$numSites>=threshold);
		output1 = output1[selection,];
		output2 = output2[selection,];
		motiflist1 = motiflist1[selection,];
		motiflist2 = motiflist2[selection,];
		motiflist = motiflist1;
	} else {
		motiflist1 = motiflist1[range,];
		motiflist2 = motiflist2[range,];
		motiflist = motiflist1;
	}

	drawDoubleAggregatedObservedExpectedCutcountPlot <- function (mobs1, mobs2,motifcenter,
				title, motiflogo,  plot1yrange=0) {

		halfwidth= (length(mobs1)-1)/2;
		motifwidth = nchar(motiflogo);
		colorcode=c('darkolivegreen','purple');   # some other colors required

		miny=min(mobs1,mobs2);
		miny=0;

		maxy=max(mobs1,mobs2);

		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		if (length(plot1yrange)<2) {
			yrange = c(miny,maxy);
		} else {
			yrange= plot1yrange;
		};

		hgt = 0;

		op <- par(mar=c(5,6,4,2) + 0.1);

		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		x1 = min(xr); x2 = max(xr);

		linea = x1 - motifoffset;
		lineb = x2 + motifoffset;

		plot((-halfwidth:halfwidth), mobs1,"n", main = title, col = colorcode[1], ylim=c(miny,maxy),
		   xlab = '(BP)',  ylab='', lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
		   xaxs='i',yaxs='i', las=1);
		mtext('Normalized Cut Counts', side=2, line=3, cex=1.5);

		#browser();
		#text(xr, rep(maxy- (maxy-miny)*0.05, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

		if (length(mobs1)!=length(mobs2)) {
			print("length of mobs1 and mobs2 differ.");
		  stop(sprintf('Error in line number 1426'));
	#		browser();
		}
		lines((-halfwidth:halfwidth), mobs1, col= colorcode[1],lwd=2.5,lty=1);    #observed 1
		lines((-halfwidth:halfwidth), mobs2, col= colorcode[2],lwd=2.5,lty=1);    #observed 2


		if (motifwidth >0) {
		lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
		lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);
		}

		legend(-halfwidth+1, y = yrange[2] - (yrange[2]-yrange[1])*0.05, bty='n',c(dataname1, dataname2), lwd=3, col = colorcode,cex=1.3);
	}


	drawDoubleLogratioPlot <- function (logratio1, logratio2, meanlog1, meanlog2, motifcenter, title, motiflogo, yrange=0) {

		motifwidth = nchar(motiflogo);
		colorcode=c('darkolivegreen','purple');   # some other colors required

		miny=min(min(logratio1),min(logratio2));
		maxy=max(max(logratio1),max(logratio2));

		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		x1 = min(xr); x2 = max(xr);

		linea = x1 - motifoffset;
		lineb = x2 + motifoffset;


		if (length(yrange)<2) {
		  yrange = c(miny,maxy);
		}

		op <- par(mar=c(5,6,4,2) + 0.1);

		linewidth= 3;
		plot((-halfwidth:halfwidth), logratio1,"n", main = title,
			 xlab = '(BP)', ylab='', ylim=yrange, lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
			 xaxs='i',yaxs='i',las=1);

		# baseline + mean line
		lines(c(-halfwidth, halfwidth), c(meanlog1,meanlog1), col= colorcode[1], lwd=2, lty=2);
		lines(c(-halfwidth, halfwidth), c(meanlog2,meanlog2), col= colorcode[2], lwd=2, lty=2);

		lines(c(-halfwidth, halfwidth), c(0,0), col= 'black', lwd=1.5, lty=1);
		#}

		lines((-halfwidth:halfwidth), logratio1,"l", col = colorcode[1],lwd=linewidth);
		lines((-halfwidth:halfwidth), logratio2,"l", col = colorcode[2],lwd=linewidth);


		# y label
		mtext('Footprinting Depth',side=2, line=3,cex=1.5);
		par(op);


		lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
		lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);

		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		text(xr, rep(maxy*0.95, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

		legend(-halfwidth+1, y = yrange[1] + (yrange[2]-yrange[1])*0.15, bty='n',c(dataname1, dataname2), lwd=3, col = colorcode,cex=1.3);


	}

	if (is.na(range[1])) {
		range = 1:nrow(motiflist);
	}
	#browser();

	for (ii in 1:length(range)) {
		motif = motiflist[ii,];
		name = as.character(motif$motif);
	#	print(name);

		sites1 = output1[ii,];
		sites2 = output2[ii,];

		numsite1 = sites1$numSites;
		numsite2 = sites2$numSites;

		name = as.character(motif$motif);
		motiflogo = as.character(motif$logo);
		motifcenter =motif$center;

		load(as.character(output1$graphdatafile[ii]));
		mobs1 = mobs;
		mexp1 = mexp;
		logratio1 =  logratio - baselogratio;
		meanlog1 = motiflogratio - baselogratio;

		load(as.character(output2$graphdatafile[ii]));
		mobs2 = mobs;
		mexp2 = mexp;
		logratio2 =  logratio - baselogratio;
		meanlog2 = motiflogratio - baselogratio;

		motifoffset=2;
		plot1yrange=0;


		halfwidth= (length(logratio1)-1)/2;
		filenamebase = make.names(name);

		if (length(logyrange)<2) {
			rangestring = "";
		} else {
			rangestring = paste(logyrange, collapse="to");
		}

		jpgfile1=file.path(outputdir,sprintf("aggregation_%s_%s_vs_%s_%dbp%s.jpg",filenamebase,dataname1, dataname2,halfwidth,rangestring));
		jpgfile2=file.path(outputdir,sprintf("logratio_%s_%s_vs_%s_%dbp%s.jpg",filenamebase,dataname1, dataname2,halfwidth,rangestring));

		pdffilename1 =  file.path(outputdir,sprintf("aggregation_%s_%s_vs_%s_%dbp%s.pdf",filenamebase,dataname1, dataname2,halfwidth,rangestring));
		pdffilename2 =  file.path(outputdir,sprintf("logratio_%s_%s_vs_%s_%dbp%s.pdf",filenamebase,dataname1, dataname2,halfwidth,rangestring));

		if (graphout) {
			jpeg(jpgfile1,width = 600, height = 600);
			drawDoubleAggregatedObservedExpectedCutcountPlot(mobs1, mobs2,motifcenter,name, motiflogo);
			dev.off();

			Cairo::CairoPDF(file = pdffilename1, width = 8, height = 8,   bg = "white");
			drawDoubleAggregatedObservedExpectedCutcountPlot(mobs1, mobs2,motifcenter,name, motiflogo);
			dev.off();
			print(jpgfile1);
			jpeg(jpgfile2,width = 600, height = 600);
			drawDoubleLogratioPlot(logratio1, logratio2, meanlog1, meanlog2, motifcenter, title, motiflogo, yrange=logyrange);
			dev.off();
			Cairo::CairoPDF(file = pdffilename2, width = 8, height = 8,   bg = "white");
			drawDoubleLogratioPlot(logratio1, logratio2, meanlog1, meanlog2, motifcenter, title, motiflogo, yrange=logyrange);
			dev.off();
		}

	}
}



#------------------ FastAggregation_FixedWidth.R-----------------------------------------

#' @useDynLib bagfoot CalcAggregationPerChromFixedWidthInC
CalcAggregationPerChromFastFixedWidthInC<- function(P) {
	oldw <- getOption("warn");
	options(warn = -1);

	region =  P$region;     # Array of regions
	density= P$density;      # Array of Tag Density

	x1 = as.integer(density$st);
	x2 = as.integer(density$ed);
	y1 = as.integer(region$st);
	y2 = as.integer(region$ed);

	if (all(is.na(P$dir))) {
	   dir = seq(from=1,to=1, length=length(y1));
	} else {
	    dir = P$dir;
	}

	width=y2[1]-y1[1]+1;
	Result=.Call("CalcAggregationPerChromFixedWidthInC", x1,x2,y1,y2,as.integer(dir), density$value);
	options(warn = oldw);
	array(data=Result, dim=c(length(y1),width));
}

CalcHeatmapFastFixedWidthC<- function(tregion, tdensity, dir=NULL) {
#browser();
  chrom = levels(factor(tregion$chr ));
  locchrom = as.character(tregion$chr);

  if (is.null(dir)) {
	dir = rep(NA, nrow(tregion));
  }
  numchrom = length(chrom);
  Plist = {};
  nregion = nrow(tregion);
  range=seq(from=1, to=nregion);

  runCalcAggregationPerChromFastPerChrom<- function(ch) {
	P = {};
	#print(ch);
	P$density = tdensity[[ch]];
	maxcoord = max(P$density$ed);
	P$region = subset(tregion, chr==ch ); # & ed < maxcoord
#	browser();
	P$dir = subset(dir, tregion$chr==ch);
	withinRange = P$region$ed < maxcoord;
	outRange = !withinRange;
	P$region = P$region[withinRange,];
	P$dir = P$dir[withinRange];
	RC =array(0,dim=c(length(withinRange), tregion$ed[1]-tregion$st[1]+1));
	if (nrow(P$region)>0) {
		R = CalcAggregationPerChromFastFixedWidthInC(P);
		RC[withinRange,] = R;
	} 	
	RC[outRange,] = array(0, dim=c(sum(outRange),tregion$ed[1]-tregion$st[1]+1));
	#browser(text="line#1623");
	RC;
  }
  chrs=unique(as.character(tregion$chr));
  result = lapply(chrs, runCalcAggregationPerChromFastPerChrom);
  do.call("rbind", result);
}


CalcHeatmapFastMatrixOutputC <- function(tregion, tdensity,dir=NULL) {
    out = CalcHeatmapFastFixedWidthC(tregion, tdensity, dir=dir);
    out;
}

#------------------ bfoot_util.R -----------------------------------------------------------
#' @useDynLib bagfoot intersect_intervals
intersect_intervals_fastInC <- function(X1,X2,Y1,Y2)
{
	Isorted1= order(X1);
	Isorted2=  order(Y1);
	.Call("intersect_intervals", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}

#' @useDynLib bagfoot intersect_intervals_index
intersect_intervals_fast_indexInC <- function(X1,X2,Y1,Y2) {
	Isorted1= order(X1);
	Isorted2=  order(Y1);
	.Call("intersect_intervals_index", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}

intersect_intervals_fast_index <- intersect_intervals_fast_indexInC;
intersect_intervals_fast <- intersect_intervals_fastInC;
intersect_intervals_really_fast = intersect_intervals_fastInC;


compare_intervals <- function(x1,x2,y1,y2) ifelse(y1>x2, -1, ifelse(y2<x1,1,0));

intersect_intervals_fast_parallel <- function(P) {

    R={};
    if (length(P$Ast) ==0  || length(P$Bst)==0) {
      return(R);
    }
    R=intersect_intervals_fast(P$Ast, P$Aed, P$Bst,P$Bed);
    R;
}

intersect_intervals_really_fast_parallel <- function(P) {
    R={};
    if (length(P$Ast) ==0  || length(P$Bst)==0) {
      return(R);
    }
    R=intersect_intervals_really_fast(P$Ast, P$Aed, P$Bst,P$Bed);
    R;
}


csv2bed <- function(csvfilename, bedfilename) {
	dat<- readcsv(filename);
	D=data.frame(chr=as.character(dat$chr), st=dat$st-1, ed=dat$ed);
	writebed(bedfilename, D, header=F);
}


readbed <- function(filename)
{

    readbed= readLines(filename, n=5);
    i=1;
    while (i < 5)
    {
        if (substr(readbed[[i]],1,3)=="chr") break;
        i=i+1;
    }
    if (i >= 5)
    {
        readbed={};
        sprintf("%s is not a readable bed file.", A);
        return;
    } else
    {
        readbed=read.csv(filename,sep="", skip=i-1, header=FALSE, colClasses=c("character","numeric","numeric"));
        names(readbed) <- c("chr","st","ed");
        readbed$st = readbed$st+1;
    }
    readbed;
}

readbgr <- function(filename)
{
    if (filename=="")
    {
        R=NULL;
        R;
    } else {
        # browser();
        R= readLines(filename, n=7);
        i=1;
        while (i < 7)
        {
            if (substr(R[[i]],1,3)=="chr") break;
            i=i+1;
        }
        if (i >= 7)
        {
            R={};
            sprintf("%s is not a readable bed file.", A);
            return;
        } else
        {
            R=read.csv(filename, sep="", header=FALSE,comment.char="", skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE);
            names(R) <- c("chr","st","ed","value");
            R$st = R$st+1;
            R;
        }
    }
}

bgr2rdensity <- function(inputfilename, outputfilename=NULL, binsize=20) {

    cat(sprintf("bgr2rdensity v.1.0.0    May 18, 2011\n"));
    # tagdensityfiles = c("","");
    #inputfilename= "/home/sjbaek/Data/data/Dalal/DHS/bgr/DHS_HeLa_r1.bgr"
    #outputfilename = "DHS_HeLa_r1.rdensity";

    if (is.null(outputfilename)) {
        outputfilename= sub('.bgr','.rdensity', inputfilename);
    }

    b1= readbgr(inputfilename);
    chrlist = levels(b1$chr);

    rd = list(density=vector("list", length(chrlist)),binsize=binsize);
    rd$read <- function(chr, addr) { rd$density[[chr]][1 + floor(addr/ binsize)]; };
    names(rd$density) <- chrlist;

    getaddr <- function(x) { r= 1+floor((x-1)/binsize) };
    invaddr <- function(r) { r= 1+(r-1)*binsize };

    for (ch in chrlist) {

        sel = (b1$chr == ch);
        maxaddr= max(b1$ed[sel]);
        maxx = 1 + floor(maxaddr/ binsize);

        st = b1$st[sel];
        ed = b1$ed[sel];
        dd = b1$value[sel];
        rrst = getaddr(st);
        rred = getaddr(ed);
        mst = invaddr(rrst);
        med = invaddr(rred);
        bb=length(rrst);
        #   for (i in 1:length(rrst)) {
        #  for (bb in seq(from=250, to=700771, by=50000))
        #{
        density.raw = vector("numeric", length=maxx);
        for (i in 1:bb) {
            if (rred[i] > rrst[i]) {
                density.raw[rrst[i]] = density.raw[rrst[i]] + dd[i] * (mst[i] + binsize - st[i]);
                density.raw[rred[i]] = density.raw[rred[i]] + dd[i] * (ed[i] - med[i]+1);
                if (rred[i] > rrst[i]+1) {
                    density.raw[(rrst[i]+1):(rred[i]-1)] = density.raw[(rrst[i]+1):(rred[i]-1)] + dd[i]*binsize;
                }
            } else {
                density.raw[rrst[i]] = density.raw[rrst[i]] + dd[i] * (ed[i] - st[i]+1);
            }
        }
        cat(paste(ch,"-",bb,":",sum(density.raw)," ", sum(dd[1:bb] * (ed[1:bb]-st[1:bb]+1)),"\n"));
        rd$density[[ch]]= as.integer(round(density.raw/20));
    }
    save(rd, file=outputfilename, compress=TRUE);

}
#

specialgrep <- function(x,y,...){
  grep(
    paste("^",
          gsub("([].^+?|[#\\-])","\\\\\\1",x)
          ,sep=""),
    y,...)
}

readcsv <- function(filename)
{
    R=read.table(filename, sep=",", header=TRUE);
  #  header=c("chr","st","ed","MaxD","AveD","Zscore");
	idx=sapply(c('chr','st'),function(x) {specialgrep(x, names(R)); });
	names(R)[idx[["chr"]]] = "chr";
	names(R)[idx[["st"]]] = "st";
	names(R)[idx[["st"]]+1] = "ed";
   # browser();
 #   nc = min(5,ncol(R)+1);
#    names(R)[2:nc] = header[1:(nc-1)];
    #readcsv$st = readcsv$st+1;
    R;
}


readcsv2 <- function(filename)
{
	library(tools)
	filebase = basename(filename);
	ext = tolower(file_ext(filebase));
	if (ext == "bed") {
		R=readbed(filename);
	} else {
		R=readcsv(filename);
	}
    R;
}


readhotspot<-function(filename)
{
    readcsv(filename);
}

readannot<-function(filename)
{
    read.csv(filename, sep=",", header=TRUE);
}

readbed_as_csv <- function(filename)
{
    hs=readbed(filename);
    lnth= length(hs$st);
    chs= data.frame(chr = hs$chr, st=hs$st+1, ed=hs$ed, MaxD=vector("numeric",length=lnth), AveD=vector("numeric",length=lnth), Zscore=vector("numeric",length=lnth));
    chs;
}

filterHotspot <- function(hotspot, category, threshold) {
    test= hotspot[[category]]>=threshold;
    R=hotspot[test,];
    R;
}


extendHotspot <- function(hotspot, ext) {
	if (ext >=0) {
		hotspot$st = hotspot$st - ext;
		hotspot$ed = hotspot$ed + ext;
		st = sapply(hotspot$st, function(x) { ifelse(x>0,x,1);});
		hotspot$st = st;
	}
	hotspot;
}


writebed <- function(filename, D, header=T)
{
    trackname=basename(filename);
    description= basename(filename);
    fid = file(filename, 'w');
    oheader= sprintf('browser hide all\ntrack name="%s" description="%s"\n', trackname, description);

	if (header) {
	    cat(paste(oheader),file=fid);
	}
    for (i in 1: length(D$st))
    {
        cat(sprintf('%s\t%d\t%d\n',as.character(D$chr[i]), D$st[i]-1, D$ed[i]),file=fid);
    }
    close(fid);
}


writeBGRHeader <- function(fp, dataname, datadescription)
{


}

writeBGRPerChrom <- function(fp, chr, count,  binsize = 20, thr=2)
{

    span= length(count);
    sum = 0;
    c=0;
    cc=0;
    prevc=0;
    meanc=0
    prevj=0;
    sz=0;
    MINTHR = 1.0;
    ii = 1;     jj = 0;
    while (ii <= span) {
        c = count[ii];
        if (c > thr) {
            jj = ii + 1;
            while (count[jj] == c && jj <= span) {
                jj=jj+1;
            }
            if (c > 0 && jj <= span) {
                if (c<30.0) {
                    cat(sprintf("%s %d %d %.1f\n", chr, ii-1, jj-1, c), file=fp);
                } else {
                    cat(sprintf("%s %d %d %.0f\n", chr, ii-1, jj-1, c), file=fp);
                }
            }
            ii = jj;
        } else {
            jj = ii;
            cc = 0;
            prevj = jj;
            prevc = 0;
            while (jj <= span) {
                if (count[jj] > thr)
                    break;
                sum = 0;
                sz = 1;
                while (jj <= span) {
                    if (count[jj] > thr || sz > binsize)
                        break;
                    sum = sum+count[jj];
                    jj=jj+1;
                    sz=sz+1;
                }
                cc=cc+1;
                meanc =  sum / (sz - 1);
                if (cc > 1 &&  (round(meanc) !=  round(prevc))) {
                    if (prevc > MINTHR && prevj <= span)
                    {
                        if (prevc<30.0) {
                            cat(sprintf("%s %d %d %.1f\n", chr, ii-1, prevj-1, prevc),file=fp);
                        } else {
                            cat(sprintf("%s %d %d %.0f\n", chr, ii-1, prevj-1, prevc), file=fp);
                        }
                    }
                    ii = prevj;
                }
                prevj = jj;
                prevc = meanc;
            }

            if (meanc > MINTHR && jj <= span) {
                if (meanc<30.0) {
                    cat(sprintf("%s %d %d %.1f\n", chr, ii-1, jj-1, meanc), file=fp);
                } else {
                    cat(sprintf("%s %d %d %.0f\n", chr, ii-1, jj-1, meanc),file = fp);
                }
            }
            ii = jj;
        }
    }
}


#' @useDynLib bagfoot readcutcount
readCutCount <- function(filename, chr) {
    maxaddr = 250000000;
    out <- .C("readcutcount", filename = as.character(datafilepath),chr= as.character(chr), span=as.integer(0), count= as.integer(vector("numeric",length=maxaddr)));
    cutcount <- out$count[1:out$span];
    cutcount;
}



file_check <- function(filelist)  {
    res = file.exists(filelist);
    ox<-ifelse(res,"O","X");
    for (i in 1:length(ox)) {
        cat(sprintf('[%s] %s\n', ox[i], filelist[[i]]));
    }
    return (all(res));
}

################ bfoot.bagplot ###############################################

#' @export
gen_bagplot<-function(dat, dataname1=dataname1, dataname2=dataname2, factor=3, pdf=T,  noBag=F, noFence=F, noMedianPoint=F) {
#	factor = 3;
	x=dat$lrtotal;
	y=dat$fdepth1-dat$fdepth2;

	trimbp <- function(sss) {
		doit <- function(sss) {
		ttt= gregexpr('_[0-9]+bp', sss)[[1]][1];
		ifelse(ttt<1, sss, substring(sss, 1,ttt-1));
		}
		sapply(sss, doit);
	}

	name = trimbp(as.character(dat$motif));

	datan = data.frame(x=x, y=y, name=as.character(name),stringsAsFactors=F);
	ss <- aplpack::compute.bagplot(datan$x, datan$y, factor=factor, approx.limit = nrow(datan));

	outlier = ss$pxy.outlier;
	outer = ss$pxy.outer;
	bag = ss$pxy.bag;
	hdepths = ss$hdepths;
	dat = datan;

	findOutliersIndex<-function(data, outlier) {
		cols <- colnames(outlier);
		dat2 <- data[,cols];
		nr = nrow(outlier);
	    rg = 1:nrow(dat2);
		matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,2] == dat2[,2]) & (outlier[ii,1]==dat2[,1])]))));
		if (length(matchingIdx) != nr) {
			stop('cannot locate outliers in the data.');
		}
		matchingIdx;
	}

	nametoshow = rep("", nrow(datan));
	nametoshowInGray = rep("", nrow(datan));
	if (!is.null(outlier)) {
		outlieridx = findOutliersIndex(datan, outlier);
		nametoshow[outlieridx] = as.character(datan$name[outlieridx]);
	}

	outeridx = findOutliersIndex(datan, outer);
	inneridx = findOutliersIndex(datan, bag);
	nametoshowInGray[outeridx] = as.character(datan$name[outeridx]);

	category = vector("character", length=nrow(datan));
	category[inneridx] = "bag";
	category[outeridx] = "fence";
	category[outlieridx ] = "outlier";

	mx = mean(datan$x);
	my = mean(datan$y);

	udist = sqrt((datan$x-mx)^2 + (datan$y-my)^2);
	bagplottable= data.frame(name=datan$name, deltahyp=datan$x, deltafd = datan$y, category=category, hdepth= hdepths, udist=udist);

	outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s",dataname2,dataname1);
	write.csv(bagplottable, file=sprintf('%s_bagplot_output.csv', outputfilename));

	labx=sprintf("Normalized Cutcount Diff. (%s-%s)",dataname2,dataname1);
	laby=sprintf("-Footprinting Depth Diff. (%s-%s)",dataname2,dataname1);

	do_draw<-function() {
		ss <- aplpack::bagplot(datan$x, datan$y, factor=factor,show.looppoints = FALSE,show.bagpoints=FALSE,
show.baghull = !noBag,  show.loophull= !noFence,  show.whiskers=F,xlab=labx,ylab=laby, col.looppoint = '#000000',cex=2,pch=15, cex.lab=1.5,  cex.axis=1.5);
		points(datan$x[outeridx], datan$y[outeridx], pch=16, col= "#332280");      # outer
		points(datan$x[inneridx], datan$y[inneridx], pch=18, col= "#000000");   # bag
		text(datan$x, datan$y, labels=nametoshow,adj=1, cex=1.5);
		text(datan$x, datan$y, labels=nametoshowInGray,col="#000000A0", adj=1, cex=0.9);
		plotrange <- par("usr");
		xmin = plotrange[1];
		xmax = plotrange[2];
		ymin = plotrange[3];
		ymax = plotrange[4];
		lines(c(xmin, xmax),c(0,0));
		lines(c(0,0),c(ymin, ymax));
	}
	do_draw();
	outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s",dataname2,dataname1);
	tiff(filename= sprintf("%s.tiff",outputfilename),compression="zip", width=2000,height=1500, units="px", pointsize=20);
	do_draw();
	dev.off();

	if (pdf) {
		pdf(file= sprintf("%s.pdf",outputfilename), width=16,height=12, pointsize=10,useDingbats=FALSE);
		do_draw();
		dev.off();
	}

}

#' @export
gen_bagplot_chisq<-function(dat, dataname1=dataname1, dataname2=dataname2, qvaluethreshold=0.05,factor=3, pdf=T,  noBag=F, noFence=F, noMedianPoint=F,pdfwidth=16, pdfheight=12, textsizefactor=1.0, printmode=F) {
#	factor = 3;
	x=dat$lrtotal;
	y=dat$fdepth1-dat$fdepth2;

	trimbp <- function(sss) {
		doit <- function(sss) {
		ttt= gregexpr('_[0-9]+bp', sss)[[1]][1];
		ifelse(ttt<1, sss, substring(sss, 1,ttt-1));
		}
		sapply(sss, doit);
	}

	mat <- data.frame(dhp=dat$lrtotal, dfd=dat$fdepth1-dat$fdepth2);
	name = trimbp(as.character(dat$motif));

	datan = data.frame(x=x, y=y, name=as.character(name),stringsAsFactors=F);
	ss <- aplpack::compute.bagplot(datan$x, datan$y, factor=factor, approx.limit = nrow(datan));

	outlier = ss$pxy.outlier;
	outer = ss$pxy.outer;
	bag = ss$pxy.bag;
	hdepths = ss$hdepths;
	dat = datan;

	findOutliersIndex<-function(data, outlier) {
		cols <- colnames(outlier);
		dat2 <- data[,cols];
		nr = nrow(outlier);
	    rg = 1:nrow(dat2);
		matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,2] == dat2[,2]) & (outlier[ii,1]==dat2[,1])]))));
		if (length(matchingIdx) != nr) {
			stop('cannot locate outliers in the data.');
		}
		matchingIdx;
	}

	nametoshow = rep("", nrow(datan));
	nametoshowInGray = rep("", nrow(datan));
	if (!is.null(outlier)) {
		outlieridx = findOutliersIndex(datan, outlier);
		nametoshow[outlieridx] = as.character(datan$name[outlieridx]);
	}

	outeridx = findOutliersIndex(datan, outer);
	inneridx = findOutliersIndex(datan, bag);
	nametoshowInGray[outeridx] = as.character(datan$name[outeridx]);

	category = vector("character", length=nrow(datan));
	category[inneridx] = "bag";
	category[outeridx] = "fence";
	category[outlieridx ] = "outlier";

	mx = mean(datan$x);
	my = mean(datan$y);

	Sx <- cov(mat);
	#D2 <- mahalanobis(mat, colMeans(mat), Sx);
	D2 <- mahalanobis(mat,apply(mat, 2, median), Sx);
	pvalue =   1-pchisq(D2, ncol(mat));
	qvalue = p.adjust(pvalue, method="BH");  # FDR adjusted p-value (BH)
	#browser();
	#browser();

	udist = sqrt((datan$x-mx)^2 + (datan$y-my)^2);
	bagplottable= data.frame(name=datan$name, deltahyp=datan$x, deltafd = datan$y, category=category, hdepth= hdepths, udist=udist,pvalue=pvalue,qvalue=qvalue);

	outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s_qvalue",dataname2,dataname1);
	write.csv(bagplottable, file=sprintf('%s_bagplot_output.csv', outputfilename));

	labx=sprintf("Normalized Cutcount Diff. (%s-%s)",dataname2,dataname1);
	laby=sprintf("-Footprinting Depth Diff. (%s-%s)",dataname2,dataname1);

	do_draw<-function() {
		ss <- aplpack::bagplot(datan$x, datan$y, factor=factor,show.looppoints = FALSE,show.bagpoints=FALSE,
show.baghull = !noBag,  show.loophull= !noFence,  show.whiskers=F,xlab=labx,ylab=laby, col.looppoint = '#000000',cex=2,pch=0, cex.lab=1.5,  cex.axis=1.5, lwd = 3);
		points(datan$x[outeridx], datan$y[outeridx], pch=16, col= "#332280");      # outer
		points(datan$x[inneridx], datan$y[inneridx], pch=18, col= "#000000");   # bag

		#PVALUETHRESHOLD = 0.05;
		points(datan$x[qvalue < qvaluethreshold], datan$y[qvalue < qvaluethreshold],cex=2, pch=15, col= "#FF0000");   #  p-value < 0.05;q

		text(datan$x, datan$y, labels=nametoshow,adj=1, cex=1.5);
		text(datan$x, datan$y, labels=nametoshowInGray,col="#000000A0", adj=1, cex=0.9);
		plotrange <- par("usr");
		xmin = plotrange[1];
		xmax = plotrange[2];
		ymin = plotrange[3];
		ymax = plotrange[4];
		lines(c(xmin, xmax),c(0,0));
		lines(c(0,0),c(ymin, ymax));
	}
	do_draw();

	outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s_qvalue%g",dataname2,dataname1,qvaluethreshold);
	tiff(filename= sprintf("%s.tiff",outputfilename),compression="zip", width=2000,height=1500, units="px", pointsize=20);
	do_draw();
	dev.off();

	if (pdf) {
		pdf(file= sprintf("%s.pdf",outputfilename), width=16,height=12, pointsize=10,useDingbats=FALSE);
		do_draw();
		dev.off();
	}

}


############################# bfoot_prep.R ####################################################


isBigEnough<-function(filename) {
	allowablesize = 300;
	if (file.exists(filename)) {
		ifelse(file.info(filename)$size > allowablesize , T,F);
	} else {
		F;
	}

}

ds239487_ptm <- proc.time();  # initial time
dlog <- function(output, ...) {
	etime = proc.time()-ds239487_ptm;
	cat(sprintf("[%g sec] %s", etime[["elapsed"]], output, ...));
	ds239487_ptm <- proc.time();
}

basenameWithoutExt<-function(filename)  {
	sub("^([^.]*).*", "\\1",basename(filename))
}

readBAMIndex<-function(bamFile, refgenome="mm9") {

	if (system("samtools", intern=F)!=1) {
		stop('Cannot run "samtools".  Install "samtools" first!');
	};

	if (refgenome=="mm9") {
		chroms = paste("chr", c(1:19,"X","Y"), sep="");
	}

	if (refgenome =="hg19") {
		chroms = paste("chr", c(1:22,"X","Y"), sep="");
	}

	baiFile = sprintf('%s.bai',bamFile);
	if (!isBigEnough(baiFile)) {
		dlog(sprintf("Generating the index file of %s\n", bamFile));
		system(sprintf('samtools index %s', bamFile));
	}
	idxout= system(sprintf('samtools idxstats %s', bamFile), intern=T);
	con <- textConnection(idxout);
	dat <- read.csv(con,sep='\t',header=F);
	names(dat) <- c("chr","maxloc","num","etc");
	dat[dat$chr %in% chroms,];
}

getMaxLocBAM<-function(bamFile) {   # returns the maximum locations of each chromosome
	dat = readBAMIndex(bamFile);
	ll = list();
	lapply(1:nrow(dat), function(ii) { ll[[as.character(dat$chr[ii])]] <<-  dat$maxloc[ii]; });
	ll;
}

#' @export
countReadsBAM<-function(bamFile, refgenome = "mm9") {
	#bamFile = fastedbam;
	dat = readBAMIndex(bamFile, refgenome= refgenome);
	ncount = sum(dat$num);
	ncount;
}

#' Make a cutcount data file from a BAM file
#'
#' This function reads a bam file and convert it to a BedGraph file which
#' contains counts of cuts (DNase) or insertions on the genomic locations.
#'
#' @param bamfile Path to the input bam file
#' @param name Name of the output bedgraph file
#' @param refgenome  Reference genome, eg) "mm9", "hg19"
#' @param atac Logical  TRUE for ATAC-seq
#' @return Bedgraph filename
#' @export
makeCutCountBAMWithName<- function(bamFile, name="", refgenome = "mm10",atac=TRUE) {  #modified Pedro
	#OutputDir =  dirname(bamFile);
	OutputDir = getwd();
	bgrfilename = file.path(OutputDir, sprintf('%s_AC_cutcount.bgr',name));
	bgrgzfilename = paste(bgrfilename,'.gz', sep='');

	bamindex = readBAMIndex(bamFile);  # Make sure that the index file is generated before creating cutcounts
	# browser();
	if (!isBigEnough(bgrgzfilename)) {
		BuildCutCountProfileBAM(bamfile = bamFile, bgrfilename=bgrfilename, ftablefile = '',  refgenome=refgenome, compress=T,atac=atac);
	};
	bgrgzfilename;
}

#' Make a cutcount data file from a BAM file
#'
#' A shorter version of makeCutCountBAMWithName
#'
#' @param bamfile Path to the input bam file
#' @param refgenome  Reference genome, eg) "mm9", "hg19"
#' @return Bedgraph filename
#' @export
makeCutCountBAM<- function(bamFile, refgenome = "mm9") {
	makeCutCountBAMWithName(bamFile, name=basenameWithoutExt(bamFile), refgenome=refgenome);
}

#' Combine two hotspot (peak) files by
#'
#' This function combines two hotspot files in the CSV format into
#' taking unions of hotspots
#'
#' @param bamfile Path to the input bam file
#' @param refgenome  Reference genome, eg) "mm9", "hg19"
#' @return Bedgraph filename
#' @export
combineTwoHotspots <- function(csvfile1, csvfile2, name1, name2) {
  # Last Edited: 2/2/2018

  Outputdir = getwd();
  pooledhotspotfilename = file.path(Outputdir, sprintf('pooled_%s_%s_hotspot.csv',name1,name2));
  hotspotfiles = list(csvfile1, csvfile2);
  hotspotlist = lapply(hotspotfiles, readcsv);
  phs = pool_hotspots3(hotspotlist);
  if (!"ID" %in% names(phs)) {
    phs = data.frame(ID=1:length(phs$chr),phs);
  }

  write.table(phs, file=pooledhotspotfilename, sep=",", row.names=F);
  pooledhotspotfilename;
}

######################################################
pool_hotspots3<- function(hotspotdata) {

  chrs = as.character(sort(unique(unlist(lapply(hotspotdata, function(x) {unique(x$chr);})))));  numData = length(hotspotdata);
  ihotspots = vector("list", length=numData);

  R={};
  if (length(hotspotdata)<2)
  {
     R=hotspotdata[[1]];
     R;
     return;
  }
  ihotspots = vector("list", length=numData);
  maxnumhotspot = 0;
  for (j in 1:length(hotspotdata))
  {
      ihotspots[[j]]= 1:length(hotspotdata[[j]]$st);
      maxnumhotspot = maxnumhotspot + length(hotspotdata[[j]]$st);
  }

  R$chr = vector("character", length=maxnumhotspot);
  R$st = vector("numeric", length=maxnumhotspot);
  R$ed = vector("numeric", length=maxnumhotspot);

  range_st= 1;
  for (i in 1:length(chrs))
  {
      jj=1;
      hotspots=vector("list");
      range_hotspots =  vector("list", length= numData);
      for (j in 1:numData) {
		  hotspots[[j]]=vector("list");
		  range_hotspots[[j]] = ihotspots[[j]][hotspotdata[[j]]$chr== chrs[i]];
		  hotspots[[j]]$st= hotspotdata[[j]]$st[range_hotspots[[j]]];
		  hotspots[[j]]$ed= hotspotdata[[j]]$ed[range_hotspots[[j]]];
      }
      st=hotspots[[1]]$st;
      ed=hotspots[[1]]$ed;

      for (j in 2:numData) {
		st=c(st,hotspots[[j]]$st);
		ed=c(ed,hotspots[[j]]$ed);
      }
	  if (length(st)>0)
	  {
		  ord= order(st);
		  st1 = st[ord];  ed1 = ed[ord];
		  leng = length(st);
		  sst = st1[1];  cmp = ed1[1];

		  chr2= vector("character", length = leng);
		  st2= vector("numeric", length = leng);
		  ed2= vector("numeric", length = leng);

		  ii=2;
		  while (ii<=leng)
		  {
				if (cmp < st1[ii])
				{
					  chr2[jj]= chrs[i];  st2[jj] = sst;     ed2[jj] = cmp;

					  jj=jj+1;
					  cmp = ed1[ii];    sst = st1[ii];
				} else {
					  cmp = max(cmp,ed1[ii]);
				}
				if (ii==leng)
				{
				  chr2[jj] = chrs[i];
				  st2[jj] = sst;
				  ed2[jj] = cmp;
				}
			ii=ii+1;
		}

		range_ed = range_st+jj-1;
		R$chr[range_st:range_ed] = chr2[1:jj];
		R$st[range_st:range_ed] = st2[1:jj];
		R$ed[range_st:range_ed] = ed2[1:jj];
		range_st = range_st+jj;
	  }
  } #i
    R$chr=R$chr[1:range_ed];
    R$st=R$st[1:range_ed];
    R$ed=R$ed[1:range_ed];

  hhs =  as.data.frame(matrix(0, nrow=range_ed, ncol=ncol(hotspotdata[[1]])));
  names(hhs) <- names(hotspotdata[[1]]);
  if ("ID" %in% names(hhs)) {
	  hhs$ID = 1:range_ed;
  }
  hhs$chr = R$chr;
  hhs$st = R$st;
  hhs$ed = R$ed;
  hhs;
}



##########################################calcNucBiasCorrectedCutcounts.R #######


readCutSitesPerChromFromBAM <- function(bamfile, chrom,shifts=c(0,0)) {

#  cat('readCutSitesPerChromFromBAM\n');

  TAGLENGTHLIMIT=1000;

  param <- Rsamtools::ScanBamParam(what=c('rname', 'pos','qwidth','mpos','isize','strand'),flag=Rsamtools::scanBamFlag(isSecondMateRead=FALSE ),
                       which = GenomicRanges::GRanges(chrom, IRanges::IRanges(1,536870912)));   # T

  bam <- data.frame(Rsamtools::scanBam(bamfile,param=param)[[1]]);
  bam=bam[abs(bam$isize) < TAGLENGTHLIMIT,];

  paired = any(bam$isize>0);
  if (paired) { # Paired-end sequencing
    st = bam$pos + ifelse(bam$strand=='+', 0, bam$qwidth + bam$isize);
    ed = bam$pos+  ifelse(bam$strand=='+', bam$isize-1, bam$qwidth-1);
  } else { # Single-end sequencing
    st = bam$pos;
    ed = bam$pos+bam$qwidth-1;

  }

  if (paired) {
    cuts = sort(c(st-1+shifts[1],ed+shifts[2]));
  } else {
    cuts= sort(ifelse(bam$strand == '+', st-1+shifts[1], ed+shifts[2]));

  }

  cuts;
}

getNucleotideString<-function(np) {
	nucleotide= c('A','C','G','T');
	NN = nucleotide;
	if (np>=2) {
	  for (i in 2:np) {
		NT = c();
		for (j in 1:length(NN)) {
		  NT = c(NT,paste(NN[j],nucleotide,sep=''))
		}
		NN=NT;
	  }
	}
	NN;
}
shiftarray<-function(map, bp) {
	bp = as.integer(bp);
	ln = length(map);
	if (bp>0) {
		newmap = c(rep(0,bp), map[1:(ln-bp)]);
	} else if (bp<0) {
		newmap = c(map[(-bp+1):ln], rep(0,-bp));
	} else {
		newmap = map;
	}
	newmap;
}

pickUnmappaleBasesByMappability<- function(mapfiledir, chr, map_seqlength=35) {

	mapfile = sprintf('%s/%sb.out',  mapfiledir, chr);
	cat(sprintf('reading mappability file: %s...\n', mapfile));
	fid1=file(mapfile,'rb');
	fileSize <- file.info(mapfile)$size;
	mappability=readBin(fid1,integer(),n=fileSize,size=1, signed = FALSE,endian='little');
	close(fid1);
	forwardcutmappability = shiftarray(mappability,-1);
	backwardcutmappability = shiftarray(mappability,map_seqlength);
	unmappable = (forwardcutmappability != 1) & (backwardcutmappability != 1);
	mappability={};
	forwardcutmappability={};
	backwardcutmappability={};
	unmappable;
}

readNucleotideCodeForChromosome<-function(chr, nmer=4, genome=NA, nuccodefileDir='.', mapdir='') {

  unmappable = NA;

  np =nmer;
  nuccodefile = file.path(nuccodefileDir,sprintf('nuccode_%s_%gmer_%s.dat',GenomeInfoDb::providerVersion(genome), np, chr));

  if (file.exists(nuccodefile)) {
	    cat(sprintf('readNucleotideCodeForChromosome: loading saved data for %s: refgenome=%s\n', chr, GenomeInfoDb::providerVersion(genome)));
        print(nuccodefile);
        nuccode = readRDS(nuccodefile);
  } else {

	  if (mapdir != '') {
  		unmappable = pickUnmappaleBasesByMappability(mapdir, chr);
	  }

      seq=as.character(genome[[chr]])  # convert to a sequence string

      nseq = nchar(seq);
      seqarray=vector("integer", nseq);
      nucleotide= c('A','C','G','T');

  #    cat(sprintf('readNucleotideCodeForChromosome:l79:%s\n', chr));
      sst <- strsplit(seq, "")[[1]]
  #    cat(sprintf('readNucleotideCodeForChromosome:l81:%s\n', chr));
      seqcode=vector(mode="integer",length=nseq)

      for (ii in 1:length(nucleotide)) {
        seqcode[sst == nucleotide[ii]]=ii;
      }
      seqcode[seqcode==0] = NA;

      if (length(unmappable)>1) {
      		cat(sprintf('masking  unmappable bases.'));
			lenmap= length(unmappable);
			lenseq= nchar(seq);

			if (lenmap > lenseq) {
			stop('The mappability file is too long.');
			}
			if (lenmap < lenseq) {
				unmappable= c(unmappable, rep(FALSE, length=lenseq-lenmap));
			}
		    seqcode[unmappable] = NA;
	  }

    #  cat(sprintf('readNucleotideCodeForChromosome:l89:%s\n', chr));
      #np = 4;
      maxcode = 4^np+1;
      bases = 4^((np-1):0);  # 1024,256,64,16,4,1

      idx=seq(from=-np/2+1, to=np/2);

      shiftsequence<-function(shift) {
        rg=(1:nseq)+shift;
        rg[rg<1] = NA;
        rg[rg>nseq] = NA;
        seqcode2 = seqcode[rg];
        seqcode2;
      };

      codes = vector(mode="integer",length=nseq);

      for (ii in 1:length(idx)) {
        codes = codes + (shiftsequence(ii)-1) * bases[ii];
      }

      codes = codes + 1;
      codes[is.na(codes)]= maxcode;


      freq = vector(mode="integer", length=maxcode);
      rl = rle(sort(codes));

      label = getNucleotideString(np);
      label = c(label,'other');
      freq[rl$values] = rl$length;
      dat=data.frame(seq=label, count=freq)

      nuccode=list(code=codes, freqtable = dat);

      saveRDS(nuccode, file = nuccodefile);
  }
  nuccode;
}

calcFrequencyTableBAMChromosome<-function(bamfile, chr, genome=NA, np=6,mapdir='',shifts=c(0,0)) {

  tempdata = sprintf('temp_data_%s.dat',digest::digest(list(bamfile, chr, GenomeInfoDb::providerVersion(genome),np, shifts, mapdir)));
  if (file.exists(tempdata)) {
#  	  cat(sprintf('calcFrequencyTableBAMChromosome: loading saved data for %s: refgenome=%s\n', chr, providerVersion(genome)));
    freqtable = readRDS(tempdata);
  } else {

#~ 	  if (mappability) {
#~ 		 mapdir = getMapDirectory(genome);
#~ 	  } else {
#~ 		  mapdir ="";
#~ 	  }

	  #unmappable = pickUnmappaleBasesByMappability(mapdir, chr);

	  cuts=readCutSitesPerChromFromBAM(bamfile, chr, shifts=shifts);  # ok

	 # cat(sprintf('calcFrequencyTableBAMChromosome: finished reading cuts for %s\n', chr));
	  nuccode<-readNucleotideCodeForChromosome(chr,nmer=np,genome=genome,mapdir=mapdir); # ok


	  cutfreq= rle(cuts);
	  tab=cbind(loc=cutfreq$values, count = cutfreq$lengths);
	  nucfreq = vector("numeric", length=length(nuccode$freqtable$count));
	  codetochange=nuccode$code[cutfreq$values-floor(np/2)];
	  for (x in 1:length(codetochange)) {
		nucfreq[codetochange[x]] = nucfreq[codetochange[x]] + cutfreq$lengths[x];
	  }

      freqtable = data.frame(count=nucfreq,refseqcount= nuccode$freqtable$count, row.names=nuccode$freqtable$seq);
      cat(sprintf('calcFrequencyTableBAMChromosome: saving the result for %s\n', chr));
	  saveRDS(freqtable, file = tempdata);
	 # browser()
  }
  freqtable;
}

loadReferenceGenome<- function(refgenome) {
	genome = NULL;
	if (refgenome=='mm9') {
      genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9;
    } else if (refgenome=='hg19') {
      genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
    } else if (refgenome=='mm10') {
	  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10;
	} else if (refgenome=='hg38') {
	  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38;
	} else if (refgenome=='rn5') {
	  genome <- BSgenome.Rnorvegicus.UCSC.rn5::BSgenome.Rnorvegicus.UCSC.rn5;
	} else if (refgenome=='rn6') {
	  genome <- BSgenome.Rnorvegicus.UCSC.rn6::BSgenome.Rnorvegicus.UCSC.rn6;
	} else {
	   cat('Supported Ref. genomes : mm9, mm10, hg19, hg38, rn5, rn6\n');
	   stop(sprintf('not supported ref. genome %s\n', refgenome));
	}
    genome;
}


calcFreqeuncyTableBAM <- function ( bamfile='',refgenome='', np=6,mapdir='',shifts=shifts) {

  tempdata = sprintf('temp_main_data_%s.dat',digest::digest(list(bamfile, refgenome,np, shifts, mapdir)));
  if (file.exists(tempdata)) {
    freqtable = readRDS(tempdata);
  } else {

    genome <- NA;
    chrominfo = loadChromosomeRange(refgenome);
    chroms = names(chrominfo);
    genome <- loadReferenceGenome(refgenome);

   # print(chroms);

   freqs=parallel::mclapply(chroms, function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts), mc.cores=2);
  #freqs=lapply(c('chr1','chr2'), function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts));

    countsum = freqs[[1]]$count;
    refseqcountsum = freqs[[1]]$refseqcount;
    if (length(freqs)>1) {
      for (ll in 2:length(freqs)) {
        countsum = countsum + freqs[[ll]]$count;
        refseqcountsum = refseqcountsum + freqs[[ll]]$refseqcount;
      }
    }

    freqtable = data.frame(count=countsum,refseqcount= refseqcountsum, row.names=row.names(freqs[[1]]));
    saveRDS(freqtable, file = tempdata);
  }
  freqtable;
}


#' Generate n-mer bias tables from a given bam file
#'
#' This function generate n-mer bias tables from a given BAM file.
#' If the mappability files directory is given, mappability is used to
#' generate the bias table.  The Reference genome must be given to run this function.
#'
#' @param bamfile Path to the input bam file
#' @param outfile  Output Filename
#' @param refgenome  Reference genome, eg) "mm9", "hg19"
#' @param np  Number of base pairs, 2 = dimer, 4 = tetramer, 6 = hexamer (default)
#' @param mapdir  Directory to the mappability files
#' @param atac  Logical, T for ATAC-seq data
#' @return Output Filename
#' @export
MakeBiasCorrectionTableBAM<- function(bamfile="", outfile="", refgenome="", np=6, mapdir='',atac=F) {
#examples:
#tab=MakeBiasCorrectionTableBAM(bamfile="/mnt/Data1/MA/sjbaek/project/ghostfootprint/data/nakedDNA/NakedDNA_SRR769954_hg19_sorted.bam", outfile="Hexamer_FT_NakedDNA_hg19_withMap.txt", refgenome="hg19", np=6, mappability=T);
#tab2=MakeBiasCorrectionTableBAM(bamfile="/mnt/Data1/MA/sjbaek/project/ghostfootprint/data/nakedDNA/NakedDNA_SRR769954_hg19_sorted.bam", outfile="Hexamer_FT_NakedDNA_hg19_withoutMap.txt", refgenome="hg19", np=6, mappability=F)
#tab3=MakeBiasCorrectionTableBAM(bamfile= "/home/sjbaek/Data/data/stamseq/DS8497/DS8497_sorted_merged.bam", outfile="Hexamer_FT_DS8497_mm9_withMap.txt", refgenome="mm9", np=6, mappability=T);

	if (!file.exists(outfile)) {
		if (!atac) {
			shifts = c(0,0);
		} else {
			shifts = c(5,-5);
			cat(sprintf("The reads are shifted by (%d,%d) base-pairs for ATAC-seq data.\n", shifts[1],shifts[2]));
		}

		freq=calcFreqeuncyTableBAM(bamfile =bamfile,  refgenome=refgenome,  np=np,mapdir=mapdir,shifts=shifts);
		#	browser();
		colnames(freq) <- c("ObCuts", "GenomicPositionCount");
		gpsum = sum(as.numeric(freq$GenomicPositionCount[-nrow(freq)]));
		totalcuts =sum(freq$ObCuts[-nrow(freq)]);   # Both sum exclude the last element (others)
		ObCutRatio = freq$ObCuts / totalcuts;
		GPRatio = freq$GenomicPositionCount / gpsum;
		#AjRatio = freq$ObCuts / freq$GenomicPositionCount;
		ObCutRatio[nrow(freq)] = NA;
		GPRatio[nrow(freq)] = NA;
		CorrectionFactor= GPRatio/ObCutRatio;
		CorrectionFactor[nrow(freq)] = 1;
		tab=cbind(freq, ObCutRatio, GPRatio, CorrectionFactor);
		write.table(tab, file=outfile);
	}
	outfile;
}


convertBiasCorrectionTable<- function(ftablefile, np) {

#examples:
#tab4=convertBiasCorrectionTable('Hexamer_FT_DS8497_mm9_withMap.txt', 4);
#tab2=convertBiasCorrectionTable('Hexamer_FT_DS8497_mm9_withMap.txt', 2);

	ftable=read.table(ftablefile);
	#if (np==2 || np==4) {
	subseq= getNucleotideString(np);
	repeatN = function(x) { paste(rep('N',x),collapse='');};
	seqNseqN = paste(repeatN((6-np)/2),subseq,repeatN((6-np)/2), sep='');

	tab=data.frame(t(sapply(1:length(seqNseqN), function(ii)  {
			   seq=seqNseqN[ii];
			   sel=CalcSeqCodeForSeq(seq)
			   subtable = ftable[sel,];
			   w=apply(subtable,2,sum);
		})));
	row.names(tab) = subseq;
	tab$CorrectionFactor = tab$GPRatio / tab$ObCutRatio;
	tab = rbind(tab, ftable[nrow(ftable),]);

	if (np == 4) {
		outputfile = sub('Hexamer','Tetramer', ftablefile);
	} else if (np == 2) {
		outputfile = sub('Hexamer','Dimer', ftablefile);
	}
	if (outputfile == ftablefile) {
		outputfile = paste('NucBias_',np,'_',ftablefile, sep='');
	}
	write.table(tab, file=outputfile);
	tab;
}


FindSeqCodeInBase2ContainingSeqInBase1<-function(code, base1=4, base2=6) {  # default
	n=base1; m=base2;

	if (code <1 | code > 4^n) {
		cat(sprintf('%d is out of range.\n', code));
	}

	diff = m-n;
	hd=diff/2;
	loop=seq(from=0, to=4^hd-1);
	as.vector(sapply(sapply(loop, function(x) x*4^(m-hd)+(code-1)*4^(hd) ), function(y) y+loop)+1);

}

CalcSeqCodeForSeq<- function(seq) {

	np = nchar(toupper(seq));
	seqchar <- strsplit(seq, "")[[1]]
	nucleotide= c('A','C','G','T');

	code = c(0);
	mult = 1;
	for (ii in np:1) {
		ch = substring(seq,ii,ii);
		if (ch == 'A')  {
#			code = code * 4;  #no change
		} else if (ch =='C') {
			code = code  + 1 * mult;
		} else if (ch =='G') {
			code = code  + 2 * mult;
		} else if (ch =='T') {
			code = code  + 3 * mult;
		} else if (ch =='N') {
			code =  c(code, code+mult,code+2*mult,code+3*mult);
		} else {
			code = NA
		}
		mult = mult*4;
	}
	code = code+1;
	code;
}

#' @export
merge_two_frequency_tables<- function(file1, file2, outfile) {
	t1<-read.table(file1);
	t2<-read.table(file2);

	if (all(t1$GenomicPositionCount==t2$GenomicPositionCount)) {
		nr= nrow(t1);
		ObCuts <- t1$ObCuts + t2$ObCuts;
		GenomicPositionCount= t1$GenomicPositionCount
		GPRatio <- c(t1$GenomicPositionCount[-nr] / sum(as.double(t1$GenomicPositionCount[-nr])), NA);
		ObCutRatio <- c(ObCuts[-nr] / sum(as.double(ObCuts[-nr])),NA);
		CorrectionFactor=c(GPRatio[-nr]/ObCutRatio[-nr], 1);

		t<-data.frame(ObCuts=ObCuts,GenomicPositionCount=GenomicPositionCount , ObCutRatio=ObCutRatio,GPRatio=GPRatio,CorrectionFactor=CorrectionFactor);

		row.names(t) <- row.names(t1);
		write.table(t, file=outfile);
	} else {
		stop('Different ref. genomes were used for Table 1 and Table 2.');
	}
}

merged_two_frequency_table = merge_two_frequency_tables;


################################ CutCountProfiles.R ########################################

loadChromosomeRange<-function(refgenome) {
# last edited: 2/2/2018

  allowableChrs = c(paste('chr',1:22,sep=''),"chrX","chrY");
  genome=loadReferenceGenome(refgenome);
  seqnames = genome@seqinfo@seqnames;
  filter= seqnames %in% allowableChrs;
  seqnames = seqnames[filter];
  seqlengths= genome@seqinfo@seqlengths[filter];
  ll <- lapply(seq_along(seqnames), function(x) c(as.integer(1), as.integer(seqlengths[x])));
  names(ll) = seqnames;
  ll;
}

getDateStamp <- function() {
    format(Sys.time(), "%m-%d-%y");
}

col2ColorCode<-function(col) {
    cc=col2rgb(col)
    sprintf('%d,%d,%d', cc[1],cc[2],cc[3]);
}

writeBGRheader<-function(filename, name, desc, option, color) {

    filep = file(filename,"w+");
    str=sprintf( "track type=bedGraph name=\"%s\" description=\"%s\" %s color=%s\n",
                 name, desc, option, color);
    cat(str, file=filep);
    close(filep);
}

#' @useDynLib bagfoot writeBGRperChromosomeInt
writeBGRperChromosomeR<-function(bgrfilename, chr, count, binsize=20, threshold=2) {
	countsize= length(count);

    out <- .C("writeBGRperChromosomeInt", filename = as.character(bgrfilename),chrom=as.character(chr),
              count=as.integer(count), span=as.integer(countsize), binsize=as.integer(binsize),
              thr=as.numeric(threshold));
}

### Input
BuildCutCountProfileBAM <- function (bamfile='', dataname='', bgrfilename='', ftablefile='', refgenome='hg19', chroms=NA, compress=T, atac=F) {

	if (!atac) {
		shifts = c(0,0);
	} else {
		shifts = c(5,-5);
		cat(sprintf("The reads are being shifted by (%d,%d) base-pairs for ATAC-seq data.\n", shifts[1],shifts[2]));
	}

	if (dataname=='') {
		dataname= sub('.bam','',basename(bamfile),ignore.case=T);
	}

	if (bgrfilename=='') {
		if (any(shifts==0)) {
			bgrfilename= file.path(sprintf("%s_adjusted_by_%d_%d_cutcount%s.bgr",dataname, shifts[1],shifts[2], ifelse(ftablefile=='','',sprintf('_%s', basename(ftablefile)))));
		} else {
			bgrfilename= file.path(sprintf("%s_cutcount%s.bgr",dataname, ifelse(ftablefile=='','',sprintf('_%s', basename(ftablefile)))));
		}
	}

	np=NA;

	if (ftablefile != '') {
		biascorrection=T;
		ftable=read.table(ftablefile);
		np= nchar(rownames(ftable)[1]);
	} else {
		biascorrection=F;
	}


   heightPixel=c(128,80,0);
   binsize=20;

   trackcolor='black';


	bgrdir = dirname(bgrfilename);
	bgrbase = basename(bgrfilename);
	bgrbasenoext = sub('.bgr','', bgrbase, ignore.case = T);
	filetype = '';


	gboption = sprintf("visibility=full maxHeightPixels=%d:%d:%d autoScale=on",
	                   heightPixel[1],heightPixel[2],heightPixel[3]);
	description = sprintf('%s %s',
	                      dataname, getDateStamp());

	chromrange = loadChromosomeRange(refgenome);

	if (is.na(chroms)) {
	  	chroms = names(chromrange);
	}

  	tempfiles=paste(file.path(bgrdir, bgrbasenoext),'_',chroms,'.bgr', sep='');
	bgrtempfiles = hash::hash(chroms,tempfiles);

	writeBGRheader(bgrfilename, dataname, description, gboption, col2ColorCode(trackcolor));

	genome <- loadReferenceGenome(refgenome);

	MakeCutCountProfilesForChromosome<-function(chr)   {
		cuts<-readCutSitesPerChromFromBAM(bamfile, chr, shifts);

		cuts <- cuts[cuts>0]; # delete negative indices

	#	cat(sprintf('Reads were shifted by %d and %d bps.\n', shifts[1],shifts[2]));
#		browser();
		cutdensity = rep(0, chromrange[[chr]][2]);
		rl = rle(cuts);  # cuts is already sorted.
		if (biascorrection==T) {
			nuccode=readNucleotideCodeForChromosome(chr, nmer=np, genome=genome, nuccodefileDir='.', mappability=T);
			nuccodeAtCuts = shiftarray(nuccode$code, np/2);
			nuccodeAtCuts[1:np]=4^np+1;
			nuccode={};
			bias = ftable$CorrectionFactor;
			rl$lengths  = rl$lengths[rl$values != 0];
			rl$values=rl$values[rl$values != 0];
			cutdensity[rl$values] = rl$lengths * bias[nuccodeAtCuts[rl$values]];
		} else {
			cutdensity[rl$values] = rl$lengths;
		}

		if (any(is.na(cutdensity))) {
			cutdensity[is.na(cutdensity)] = 0;
		}
	#
		writeBGRperChromosomeR(bgrtempfiles[[chr]], chr, cutdensity, binsize=20, threshold=0);
	#	cat(sprintf('creating cut count profiles for %s...\n', chr));
		cutdensity={};
	}

#	l<-parallel::mclapply(sample(chroms), MakeCutCountProfilesForChromosome, mc.cores=6);
	l<-lapply(sample(chroms), MakeCutCountProfilesForChromosome);

	for (chr in chroms) {
		tempfile = bgrtempfiles[[chr]];
		if (file.exists(tempfile)) {
	    	cmd1 = sprintf('cat %s >>%s',tempfile, bgrfilename);
		    cmd2 = sprintf('rm %s',tempfile);
		    system(cmd1, intern=TRUE);
		    system(cmd2, intern=TRUE);
	    }
	}

	if (compress==TRUE) {
	    system(sprintf("gzip -f %s", bgrfilename), wait=FALSE);
	}
}


#' @export				
GFootSingleRun=function(cutcount=NA, sitefile=NA, motifDB=NA, gfootOption=NA, name=cutcount@name, range=NA, graphout=F, yrange=c(-2.2,0.5), mc.cores= getOption("mc.cores", 2L)) {
	
				#obj <- loadCutcount(obj);	
	sitefilename = extractBaseName(basename(sitefile));
	outputcsvfile = sprintf("%s_On_%s_footprint_depth_table.csv", name, sitefilename);
	if (!file.exists(outputcsvfile)) {
		cutcountonsite= readCutCountSites(cutcount, site_file=sitefile);   #read cutcount
		cat(sprintf("Num. of CPU cores to use: %d\n",mc.cores));
		
		MCCORES= mc.cores;
		halfwidths=c(50);
		normalization="default";
		outputdir = file.path(".", sprintf("OUTPUT_%s", name));
		cachedir = file.path(".", sprintf("CACHE_%s", name));
			
		drawMotifAggPlotOnMotifSetsForMultipleRangesForSingleRun(
			sitefiledir=motifDB@directory,  
			cache= cachedir,                
			outputdir=outputdir,            
			dataname = name,
			cutcountdata = cutcountonsite,
			ncutcounts = c(cutcount@count),
			normalization = normalization,
			sitefile = sitefile,
			motiflist = motifDB@motiflistfile,
			ftablefile = gfootOption@biasfile,
			nuccodefilePattern = gfootOption@hexamer_pattern,
			halfwidthlist = halfwidths,
			range=range,
			graphout=graphout,
			yrange=yrange
		);

		outputfile = scatterDataOutGroupsSingleRun(
			outputdir=  outputdir ,
			dataname = name,
			motiflist = motifDB@motiflistfile,
			sitefile = sitefile,
			halfwidthlist = halfwidths,threshold = 0,
			range=range);
		return(outputfile);
	} else {
		cat(sprintf("%s exists: skipping...\n",outputcsvfile)	);
		return(outputcsvfile);
	}	
}
	


drawMotifAggPlotOnMotifSetsForMultipleRangesForSingleRun <- function(
	sitefiledir='',
	cache='',
	outputdir='',
	dataname = '',
	cutcountdata = '',
	normalization="default",
	ncutcounts = NA,
	sitefile = '',
	motiflist = '',
	ftablefile = Fed_mm9_withMap,
	nuccodefilePattern = nuccode_hexamer_withMappability_35mer_mm9,
	halfwidthlist = c(50,100),
	range= NA,
	yrange=c(-2.2,0.5),
	graphout=F) {

	threshold = 0;
	outputdir_base = outputdir;

 	if (!exists("ftablefile")) { # read once,
 	  freqtable<-readFrequencyTable(ftablefile);np = round(log(nrow(freqtable)-1)/log(4));
 	} else {
	  np = 6;   # by default
	}

  	chrs = paste("chr", c(1:22, "X","Y"),sep='');
	nuccodefiles= sapply(chrs,function(x) gsub("{chr}",x,nuccodefilePattern,fixed=T));
	doexist=file.exists(nuccodefiles);
	chrs= chrs[doexist];
	nuccodefiles=nuccodefiles[doexist];
	nuccodes= parallel::mclapply(nuccodefiles, function(nuccodefile) {
			nuccode=readNucleotideCodeForChromosomeForCuts(nuccodefile,np);
			nuccode; }, mc.cores=MCCORES,mc.preschedule=F);

	if (normalization=="default") {
	#	browser();
		oldncutcounts= ncutcounts;
		ncutcounts =  countCutcountOnSites(cutcountdata@cutcount,sitefile);
	}

	if (is.na(ncutcounts[1])) {
		stop("Error in drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons: ncutcounts should not be zero.")
	}

	for (halfwidth in halfwidthlist) {
			out=getOutputDirAndFilename(outputdir_base, dataname, sitefile, halfwidth);
			if (!file.exists(out$filename)) {
				run<-drawMotifAggPlotOnMotifSets(dataname = dataname,
					cutcount = cutcountdata@cutcount,
					ncutcount = ncutcounts,
					hotspotfile = sitefile,
					motiflist = read.table(motiflist,stringsAsFactors=F) ,
					freqtablefile=ftablefile,
					nuccodes = nuccodes,
					sitefiledir=sitefiledir,
					cache=cache,
					outputdir=out$dir,
					halfwidth=halfwidth,
					logyrange=yrange,
					range=range,
					graphout=graphout);
				gc();
			} else {
				cat(sprintf('%s exists.', out$filename));
			}
			
	}  # for (half

}



scatterDataOutGroupsSingleRun<-function(
	outputdir='',
	dataname = '',
	motiflist = '',
	sitefile = '',
	halfwidthlist = c(50,100),threshold = 0,
	range=NA) {

    outputfiles=c();
	
	for (halfwidth in halfwidthlist) {
	
		sitefilename = extractBaseName(basename(sitefile));
		motiflistfile =  motiflist
		
		out=getOutputDirAndFilename(outputdir, dataname, sitefile, halfwidth);
		outputfile=out$filename;
		
		title = sprintf('Footprinting Aggregation Plot (%s)', dataname);

		datadir =  sprintf("%s_%dbp",dataname, halfwidth);
		MAoutputdir = file.path(outputdir, 'comparison');

		outputcsvfile = sprintf("%s_On_%s_footprint_depth_table.csv", dataname, sitefilename);
		if (!file.exists(outputcsvfile)) {
			outfilename=ScatterDataOutSingleRun(dataname=dataname,
				 motiflistfile=motiflistfile,
				outputfile=outputfile,
				title = title,threshold=threshold,datadir=datadir, outputcsv=outputcsvfile,
				range=range	);
		}		
		outputfiles=c(outputfiles, file.path(getwd(),basename(outputcsvfile)));
	}	
	outputfiles;
}

ScatterDataOutSingleRun <- function(dataname="",  motiflistfile="",
					outputfile="",title = "",threshold=0,datadir='', outputcsv=' ',range=NA) {

	motiflist0=read.table(motiflistfile, stringsAsFactors=F);

	calcObsMeans <- function (mobs, motifcenter, motiflogo, meanlog, numsite) {
		motifoffset=2;
		halfwidth= (length(mobs)-1)/2;
		motifwidth = nchar(motiflogo);
		xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
		hgt = 0;
		x1 = min(xr); x2 = max(xr);
		linea = x1 - motifoffset;
		lineb = x2 + motifoffset;
		motifregion = linea:lineb  + halfwidth + 1;
		flanking = setdiff(1:length(mobs), motifregion);
		c1=log(sum(mobs[motifregion]),2);
		c2=log(sum(mobs[flanking]),2);
		c3=log(sum(mobs),2);
		meancc = mean(mobs);
	 	data.frame(logmotif=c1, logflanking=c2, logtotal=c3, fdepth= meanlog,meancc= meancc, numsite);
	}

	calcprog<-function(ii) {
		motif = motiflist[ii,];
		name = as.character(motif$motif);
		print(name);
		#rank = difforder[ii];
	#	browser();
		sites = output[ii,];
				
		numsite = sites$numSites;
		
		motif = motiflist[ii,];
		name = as.character(motif$motif);
		motiflogo = as.character(motif$logo);
		motifcenter =motif$center;

		load(as.character(output$graphdatafile[ii]));
#		mobs1 = mobs;
#		mexp1 = mexp;
#		logratio1 =  logratio - baselogratio;
		meanlog = motiflogratio - baselogratio;
		
#		plot1yrange=0;
#		halfwidth= (length(logratio)-1)/2;
#		filenamebase = make.names(name);
		#if (numsite < 1) {
	#		data.frame(logmotif=0, logflanking=0, logtotal=0, fdepth= 0,meancc= 0, numsite);
	#	} else {
	#	 browser();
			return(calcObsMeans(mobs, motifcenter,motiflogo, meanlog, numsite));
	#	}
	}

#	if (!is.na(range[1])) {
#		motiflist <- motiflist[range,];
#	}
#	browser();
	cat(sprintf("outputfile = %s\n", outputfile));
	
	#browser();	
	output = read.table(outputfile,sep='\t',stringsAsFactors=F);  #200
	
	selection= (output$numSites>=threshold & output$numSites>0);

	output = output[selection,];
	if (is.na(range)) {
		range = 1:length(selection);
	}
	motiflist = motiflist0[range[selection],];
	#browser();
	cal<-do.call("rbind", lapply(1:nrow(motiflist), calcprog));
	tabl = data.frame(motiflist, cal);
	write.table(tabl, file=outputcsv)
	outputcsv;
}
			
#' @export
RunPairwiseComparison <- function(file1 = '', file2= '',name1= '', name2= '', sitename = '', bagplot=T) {
		#sitefilename = extractBaseName(basename(sitefile));
		outputcsvfile = sprintf("%s_vs_%s_On_%s_footprint_depth_table.csv", name1, name2, sitename);

		dat1 =  read.table(file1);
		dat2 =  read.table(file2);
		if (nrow(dat1) != nrow(dat2)) {
			stop(sprintf('%s and %s have different numbers of rows.',file1,file2 ))
		}
		#lrmotif" "lrflanking" "lrtotal" "fdepth1" "fdepth2" "meancc1" "meancc2" "numsite"
		lrmotif = dat2$logmotif - dat1$logmotif;
		lrflanking = dat2$logflanking - dat1$logflanking;
		lrtotal = dat2$logtotal - dat1$logtotal;
		fdepth1 = dat1$fdepth;
		fdepth2 = dat2$fdepth;
		meancc1 = dat1$meancc;
		meancc2 = dat2$meancc;
		if (any(dat1$numsite != dat2$numsite)) {
			stop(sprintf('%s and %s have different motif sites.',file1,file2 ));
		}
		numsite = dat1$numsite;

		newtab = data.frame(dat1[,1:7], lrmotif, lrflanking, lrtotal, fdepth1, fdepth2, meancc1, meancc2, numsite);
		
		if (bagplot) {
				gen_bagplot_chisq(newtab, dataname1=name1, dataname2=name2);	
		}
		write.table(newtab, file=outputcsvfile)

}

