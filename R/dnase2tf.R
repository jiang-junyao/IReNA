Get_hotspot<-function(peak_counts){
  ATAC2<-ATAC1[1:(nrow(peak_counts)-7), 1:3]
  return(ATAC2)
}


dnase2tf <- function(datafilepath, filename, mapfiledir='',outputfilepath, biascorrection="none", assemseqdir='', dftfilename='', maxw=30,
                     minw=6,z_threshold =0, numworker = 2, FDRs=NA, paired=F) {

  if (is.na(FDRs)) {
    FDRs=c(0, 0.001, 0.005, 0.01, 0.05, 0.1 ,0.5, 1);
  }

  ##% Default Values
  numrep = 1;  # Number of repetition for random footprint generation
  max_fp_chr = 1600000;  # max footprint per chrom
  max_fp_hotspot = 10000; #max footprint per hotspot

  readBAMchr <- function(bamfile, chrom) {

    TAGLENGTHLIMIT=1000;
    param <- ScanBamParam(what=c('flag','rname', 'pos','qwidth','mpos','isize','strand','seq'),
                          flag=scanBamFlag(isFirstMateRead=F),
                          which = GRanges(chrom, IRanges(1,536870912)));
    #which = GRanges(chroms, IRanges(1,.Machine$integer.max)));
    baifile = paste(bamfile,'.bai',sep='');
    if (!file.exists(baifile)) {
      stop(sprintf('failed to load BAM index : %s', baifile));
    }
    bam <- data.frame(scanBam(bamfile,param=param)[[1]]);
    bam=bam[bam$isize < TAGLENGTHLIMIT,];

    if (any(bam$isize>0)) { # Paired-end sequencing
      st = bam$pos + ifelse(bam$strand=='+', 0, bam$qwidth + bam$isize);
      ed = bam$pos+  ifelse(bam$strand=='+', bam$isize-1, bam$qwidth-1);
    } else { # Single-end sequencing
      st = bam$pos
      ed = bam$pos+bam$qwidth-1;
    }
    B=cbind(st,ed);
    tab=data.frame(chr=rep(as.character(chrom),nrow(bam)), st=apply(B,1,min), ed=apply(B,1,max), dir=ifelse(bam$strand=='+','F','R'));
    tab;
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
      sprintf("%s is not a readable bed file.", filename);
      return;
    } else
    {
      readbed=read.csv(filename,sep="", skip=i-1, header=FALSE, colClasses=c("character","numeric","numeric"));
      names(readbed)[1:3] <- c("chr","st","ed");
      readbed$st = readbed$st+1;
    }
    readbed;
  }


  readAssembly<-function(filepath, chr) {
    asmfile=file.path(filepath, paste(chr, ".fa",sep=""));
    cat(sprintf('reading : %s...\n', asmfile));
    seq=toupper(paste(scan(file=asmfile, what=character(), skip=1, strip.white=T),collapse=""));
    seq;
  }

  readDFTfile <- function(filename) {
    sss<-read.csv(filename,sep='\t');
    colnames(sss)[1]="nuc";
    sss[,1] = as.character(sss[,1]);
    ditable = sss[1:10,];
    tettable = sss[11:nrow(sss),];

    calcBiasCorrectionRatio<-function(tab)  {
      correction=hash();
      for (kk in 1:nrow(tab)) {
        nc = tab$nuc[kk];
        w=strsplit(nc,"/")[[1]];
        correctionratio=tab$expected[kk]/tab$observed[kk];
        if (length(w)<2) {
          correction[[nc]]= correctionratio;
        } else {
          correction[[w[1]]]= correctionratio;
          correction[[w[2]]]= correctionratio;
        }
      }
      correction;
    }

    diBCR= calcBiasCorrectionRatio(ditable);
    tetraBCR=calcBiasCorrectionRatio(tettable);
    list(diBCR=diBCR, tetraBCR=tetraBCR);
  }


  splitfunction <- function(ccount,cmap, rstart, rend,prob, minw, maxw, pv_threshold)
  {
    W1 = vector("numeric", length=max_fp_hotspot);
    W2 = vector("numeric", length=max_fp_hotspot);
    W3 = vector("numeric", length=max_fp_hotspot);
    W4 = vector("numeric", length=max_fp_hotspot);
    nfootprint=0;

    #  print("splitfunction()\n"); browser();
    R<- .C("splitregions4c",
           as.numeric(ccount),
           as.integer(cmap),
           as.integer(rstart),
           as.integer(rend),
           as.numeric(prob),
           as.integer(minw),
           as.integer(maxw),
           as.numeric(pv_threshold),
           W1=as.numeric(W1),
           W2=as.numeric(W2),
           W3=as.numeric(W3),
           W4=as.numeric(W4),
           n=as.integer(nfootprint));
    #  R;
    if (R$n>0) {
      T = cbind(R$W1[1:R$n],R$W2[1:R$n],R$W3[1:R$n],R$W4[1:R$n] );
    } else {
      T = {};
    };
    T;
  }

  readDinucR<-function(filename) {
    fileSize <- file.info(filename)$size;
    out=.Call("readDinuc", normalizePath(filename),fileSize)
    out$dinuc[1:out$length];
  }

  read_map<- function(mapfiledir, chr) {
    mapfile = sprintf('%s/%sb.out',  mapfiledir, chr);
    cat(sprintf('reading mappability file: %s...\n', mapfile));
    fid1=file(mapfile,'rb');
    fileSize <- file.info(mapfile)$size;
    mappability=readBin(fid1,integer(),n=fileSize,size=1, signed = FALSE,endian='little');
    close(fid1);
    mappability;
  }

  randomf <- function(maxv, dummy, n) {1+floor(runif(n,max=maxv))};
  ##% Arguments Handling
  ### Settings

  cat(sprintf('filename:%s\n',filename));
  cat(sprintf('data file path:%s\n',datafilepath));
  cat(sprintf('map file path:%s\n',mapfiledir));
  cat(sprintf('output file path:%s\n',outputfilepath));
  cat(sprintf('assembly file path:%s\n',assemseqdir));
  cat(sprintf('dft filename:%s\n',dftfilename));
  cat(sprintf('max. footprint width:%d\n',maxw));
  cat(sprintf('min. footprint width:%d\n',minw));
  cat(sprintf('z-score threshold:%f\n',z_threshold));

  dinuc_bias_correction=F;
  tetranuc_bias_correction=F;

  diBCR={};
  tetraBCR={};

  #browser();

  ### Initialization
  if (assemseqdir!='' && dftfilename!='') {
    R=readDFTfile(dftfilename);
    diBCR= R$diBCR;
    tetraBCR= R$tetraBCR;
    if (biascorrection=="dimer")  {
      dinuc_bias_correction=T;
      tetranuc_bias_correction=F;
      print("Dinucleotide bias correction is enabled.\n");
    } else if (biascorrection=="tetramer") {
      dinuc_bias_correction=F;
      tetranuc_bias_correction=T;
      print("Tetranucleotide bias correction is enabled.\n");
    }
  }
  #splitfunction=@splitRegions4;

  bgwsize= maxw*3;
  prob = maxw/bgwsize;

  ##% reads a hotspot file
  endingWith <- function(str, pattern) {
    return (substring(str, nchar(str)-nchar(pattern)+1, nchar(str)) == pattern);
  }

  if (endingWith(tolower(filename),".bed") || endingWith(tolower(filename),".bgr")) {
    hotspot_all =readbed(filename);
  } else if (endingWith(tolower(filename),".csv")) {
    hotspot_all =read.csv(filename); names(hotspot_all)[2:4] = c('chr','st','ed');
  } else {
    stop(sprintf("Unsupported file format :%s", filename));
  }

  chroms=levels(as.factor(hotspot_all$chr));
  numchrom = length(chroms);

  #browser();

  callFootprintsPerChrom<-function(ch) {
    chr = chroms[ch];
    rstart = 0; rend =0;
    dcount =0;
    W2={};
    footprintno=1;
    fplist= array(0, dim=c(max_fp_chr ,5));
    fprandom = array(0, dim=c(max_fp_chr ,5));
    chr=chroms[ch];

    dinuc={};

    ##  reads the sequence assembly of the chromosome
    if (dinuc_bias_correction || tetranuc_bias_correction) {
      seqassembly=readAssembly(assemseqdir,chr);
    }

    BCR2 <- function(rstart,rend) {
      subseq= substr(seqassembly,rstart-1,rend);
      dinucs=sapply(1:(nchar(subseq)-1), function(x) {substring(subseq,x,x+1);});
      bcr2=sapply(dinucs,function(dinuc) {if (has.key(dinuc,diBCR)) R=diBCR[[dinuc]] else R=1; R; });
      bcr2;
    }

    BCR4<- function(rstart,rend) {
      subseq= substr(seqassembly,rstart-2,rend+1);
      tetranucs=sapply(1:(nchar(subseq)-3), function(x) {substring(subseq,x,x+3);});
      bcr4=sapply(tetranucs,function(tetranuc) {if (has.key(tetranuc,tetraBCR)) R=tetraBCR[[tetranuc]] else R=1; R; });
      bcr4;
    }

    hotspotarray={};

    subset_chr= hotspot_all$chr==chr;
    hotspotarray$st = hotspot_all$st[subset_chr];
    hotspotarray$ed = hotspot_all$ed[subset_chr];
    hotspotarray$chr = hotspot_all$chr[subset_chr];
    hotspotarray$length= length(hotspot_all$st);

    if (length(hotspotarray$st)==0)
      next;

    if (substring(tolower(datafilepath), nchar(datafilepath)-3,nchar(datafilepath))=='.bam') {
      datafile = datafilepath;
      CC = readBAMchr(datafile, chr);
      cc_st = CC[,2];
      cc_ed = CC[,3];
      cc_dir= CC[,4];
      ntags = length(cc_dir);
    } else {
      datafile = sprintf('%s_%s.txt', datafilepath,  chr);
      cat(sprintf('reading data file: %s...\n', datafile))
      CC = fread(datafile,verbose=F);   # Much faster than read.csv
      cc_st = CC$V2;
      cc_ed = CC$V3;
      cc_dir= CC$V4;
      ntags = length(cc_dir);
    }
    rm(CC);

    maxcoord=max(max(cc_ed),max(hotspotarray$ed))+500;
    data= vector("numeric",length= maxcoord);
    range=1:ntags;

    if (paired) { # if  paired-end reads
      for (jj in range) {
        hst=cc_st[jj];
        hed=cc_ed[jj]+1;
        data[hst] = data[hst] + 1;
        data[hed] = data[hed] + 1;
      }
    } else {
      for (jj in range[cc_dir == 'F']) {
        hst=cc_st[jj];
        data[hst] = data[hst] + 1;
      }
      for (jj in range[cc_dir == 'R']) {
        hed=cc_ed[jj]+1;
        data[hed] = data[hed] + 1;
      }
    }

    numhotspot = length(hotspotarray$st);

    #numhotspot=20;  #### DEBUG

    regions = array(0,dim=c(numhotspot, 2));

    ##%  Searches footprint candidates from each hotspot regions
    cdata_list = vector("list",length= numhotspot);

    for (iii in 1:numhotspot) {
      #browser();
      #print(iii);
      regstart=hotspotarray$st[iii];
      regend = hotspotarray$ed[iii];

      if ((regstart <=1) || (regend - regstart < 1))
        next;

      pp=data[regstart:regend];
      nonzerocuts = which(pp>0);
      if (length(nonzerocuts)==0) {
        rstart = hotspotarray$st[iii];
        rend = hotspotarray$ed[iii];
      } else {
        aa = nonzerocuts[1];
        bb = nonzerocuts[length(nonzerocuts)];
        rstart = max(hotspotarray$st[iii]+aa-1, hotspotarray$st[iii]);
        rend = min(hotspotarray$ed[iii], hotspotarray$st[iii]+bb-1);
      }


      if (dinuc_bias_correction) {        # With dinucleotide correction

        #	browser();
        data[rstart:rend]= data[rstart:rend]*BCR2(rstart,rend);
      }

      if (tetranuc_bias_correction) {
        data[rstart:rend]= data[rstart:rend]*BCR4(rstart,rend);
      }

      regions[iii,1] = rstart;
      regions[iii,2] = rend;
      cdata_list[[iii]] = cumsum(data[rstart:rend]);
    }
    rm(data);

    ##% reads a mappability file
    if (mapfiledir!='') {
      cat(sprintf('reading the mappability file...%s\n', chr));
      mappability = read_map(mapfiledir, chr);
      mappa = mappability==1;
      maxaddr=length(mappa);
      rm(mappability);
    } else {
      maxaddr= max(hotspotarray$ed)+maxw+1;
      mappa = rep(1,maxaddr);
    }
    cmap = cumsum(mappa);

    #	 	browser();

    mapfp <- function(a,b) {cmap[b]-cmap[a-1];}; # Sum of cutcounts between a and b
    cat(sprintf('chromosome :%s\n',chr));

    cutcountsum = vector("numeric", length= numhotspot);

    for (iii in 1:numhotspot) {
      rstart = regions[iii,1];
      rend = regions[iii,2];
      if (rend> maxaddr || (rend-rstart+1<maxw ) )
        next;

      #print(iii);
      #browser();
      cutcountsum[iii]= cdata_list[[iii]][length(cdata_list[[iii]])]
      cutcountsum[is.na(cutcountsum)]=0;

      cat(sprintf('Hotspot no. %d: %s hotspot=[%d,%d] mapsum=%d ccsum=%f\n', iii, chr,
                  rstart,rend, mapfp(rstart,rend),cutcountsum[iii]));


      if (cutcountsum[iii]<20) {
        W={};
      } else {
        W=splitfunction(cdata_list[[iii]],cmap[rstart:rend], rstart,rend,prob,minw,maxw,z_threshold);
      }
      sz=0;

      if (!(is.null(W))) {
        if (nrow(W)>0) {
          selected = W[,4]< z_threshold;
          if (any(selected)) {
            sz = sum(selected);
            #fplist = rbind(fplist,cbind(subset(W,selected), rep(1,sz)  ));
            fplist[footprintno:(footprintno+sz-1),] = cbind(subset(W,selected), rep(ch,sz)  );
            footprintno=footprintno+sz;
          }
        }
      }
      cat(sprintf('Hotspot no. %d [%s]: %d footprint candididates \n', iii, chr, sz ));

    }
    rm(cdata_list);

    randomfootprintno = 1;
    #for rep=1:numrep;
    crandom_list = list();
    for (iii in 1:numhotspot) {
      regstart=regions[iii,1];
      regend = regions[iii,2];
      if (regend > maxaddr)
        next;

      if ((regstart <=1) || (regend - regstart < 1))
        next;

      numfp = floor(cutcountsum[iii]);

      if (numfp<1)
        next;

      mapregion = mappa[regstart:regend];
      reg = 1:(regend-regstart+1);
      mapidx =  reg[mapregion];

      rcutcount = vector("numeric",length=regend-regstart+1);

      #browser();

      if (length(mapidx)<1) {
        randomcutcount=randomf(regend-regstart+1,1,numfp);
        for (ii in 1:length(randomcutcount)){
          tmp=randomcutcount[ii];
          rcutcount[tmp] = rcutcount[tmp]+1;
        }
      } else {
        randomcutcount=randomf(length(mapidx),1,numfp);
        for (ii  in 1:length(randomcutcount)) {
          tmp=mapidx[randomcutcount[ii]];
          rcutcount[tmp] = rcutcount[tmp]+1;
        }
      }

      if (dinuc_bias_correction) {        # With dinucleotide correction
        rcutcount= rcutcount* BCR2(regstart,regend);
      }

      if (tetranuc_bias_correction) {        # With tetranucleotide correction
        rcutcount= rcutcount* BCR4(regstart,regend);
      }

      crandom_list[[iii]] = cumsum(rcutcount);
    } #iii

    for (iii in 1:numhotspot) {
      regstart=regions[iii,1];
      regend = regions[iii,2];
      if (regend > maxaddr)
        next;

      if ((rend-rstart+1<maxw ) || cutcountsum[iii]<20) {
        W2={};
      } else {
        W2= splitfunction(crandom_list[[iii]],cmap[regstart:regend], regstart,regend,prob,minw,maxw,z_threshold);
      }

      if (!( is.null(W2) )) {
        if (nrow(W2)>0) {
          selected_random = W2[,4]< z_threshold;
          if (any(selected_random)) {
            sz = sum(selected_random);
            fprandom[randomfootprintno:(randomfootprintno+sz-1),] = cbind(subset(W2,selected_random), rep(ch,sz)  );
            randomfootprintno=randomfootprintno+sz;
          }
        }
      }
    } #iii

    list(fplistResult = fplist[1:(footprintno-1),],fprandomResult= fprandom[1:(randomfootprintno-1),]);
  } # end of callFootprintsPerChrom

  #	browser();

  Result= mclapply(1:numchrom, callFootprintsPerChrom, mc.cores=numworker);
  #   Result= lapply(1:numchrom, callFootprintsPerChrom);

  #Result=callFootprintsPerChrom(4);
  fplist={};
  fprandom={};

  ### Print footprints thresholded by FDRs
  #	browser();
  for (ch in 1:length(chroms)) {
    fplist = rbind(fplist, Result[[ch]]$fplistResult);
    fprandom = rbind(fprandom, Result[[ch]]$fprandomResult);
  }

  for (iii in 1:length(FDRs)) {
    fdr=FDRs[iii];
    fid1 = file(sprintf('%s_fdr%f.bgr',outputfilepath,fdr),'w');
    fid2 = file(sprintf('%s_fdr%f.bed',outputfilepath,fdr),'w');
    selection = fplist[,4] <= quantile(fprandom[,4],fdr);
    cat(sprintf('track type=bedGraph name="%s_FDR%f" description="%s_FDR%f"  visibility=full color=92,0,161 altColor=0,100,200\n',
                outputfilepath,fdr,outputfilepath,fdr ),file=fid1);
    cat(sprintf('chr\tstart\tstop\tp-value\n'), file=fid2);
    range = 1:nrow(fplist);
    for (i in range[selection]) {
      cat(sprintf('%s\t%d\t%d\t%f\n',chroms[fplist[i,5]], fplist[i,1]-1, fplist[i,2],fplist[i,4]), file=fid1);#bgr zero-based
      cat(sprintf('%s\t%d\t%d\t%f\n',chroms[fplist[i,5]], fplist[i,1]-1, fplist[i,2],fplist[i,4]), file=fid2);#bgr zero-based
    }
    close(fid1);
    close(fid2);
  }
}
