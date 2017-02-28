###########################################################################################
# 20111116
#            produce a pdf script that summarizes a dataset produced by MaxQuant
#
#
############################################################################################
# required functions:
#
# summaryPDF           - the main function that does the whole job and produces
#                        the latex code for the pdf
# summaryPDF.Sites     - function to make some slides on modification sites
# siteOverlap          - helper function for 'summaryPDF.Sites'; calculates the overlap
#                        of detected/quantifies modification sites accross different
#                        experiments
# rmConRev             - removes contaminants and reverse hits from a MaxQuant table
#                        like 'proteinGroups', 'peptides', ...
#                      - estimates Q-values
# fancyBarplot         - make a nice barplot
# fancyDensPlot        - make a nice density plot
# sigRatioPlot.ph      - plots log(ratios, 2) vs. log(intensies, 10) of protein groups
#                        or modification sites
# checkIncorporation   - produces some plots illustrating the incorporation level of
#                        SILAC AAs
# barplotInt           - helper function for 'checkIncorporation'
# q.value              - estimated Q-values
# getProteinState      - determine detected/quantified protein groups in certain
#                        experiments defined in the 'experimentalDesign' file
# significanceA        - see MQ paper
# significanceB        - see MQ paper
# erfc                 - helper function; implements complementary error
#                        function of a gaussian distribution
# getNRSites           - removes redundancy in modification site tables
# compareExperiments   - determines overlap of detected protein groups and peptides
#                      - does it separately for light/medium/heavy versions of the peptides/protein groups
#                      - detemines how many peptides/protein groups were detected in how many experiments
# splitMultipleProteinIDsInEvidence - helper function for 'compareExperiments'
# checkSeparation      - determines the number of non-redundant peptide sequecences in each LC-MS run
# venndiagram          - function to plot non-proportional venn diagrams
# olReport             - helper function for 'venndiagram'
# my.col2rgb           - helper function
#
############################################################################################
#
# changelog: 20111116
#            20120314 preparations for building the R package
#            20120329 function 'capwords' to fix column names
#
#            20120625 major update for compatibility to 1.2.7.5
#            20131217 MQ 1.4.1.2 compatibility
#            20140625 MQ 1.5.0.0 compatibility
#            20140822 MQ 1.5.0.30 compatibility
#            20140922 MQ 1.5.1.0 compatibility
#
#
############################################################################################

## load dll if working with Windows
if(Sys.info()["sysname"] == "Windows") {

    if(Sys.info()["machine"] == "x86-64") {
        ##dyn.load("H:/Projects/MaxQuantParser/C/maxquantparser_x86_64.dll")
        ##dyn.load("T:/User/Karsten/MaxQuantParser/C/maxquantparser_x86_64.dll")
	dyn.load("C/maxquantparser_x86_64.dll")
    } else
    ##dyn.load("H:/Projects/MaxQuantParser/C/maxquantparser.dll")
    ##dyn.load("T:/User/Karsten/MaxQuantParser/C/maxquantparser.dll")
    dyn.load("C/maxquantparser.dll")
}


#############################################################################################################
#
#                                          summaryPDF
#
# files needed:  proteinGroups.txt
#                peptides.txt
#                evidence.txt
#                (msms.txt)
#                modifiedPeptides.txt
#                parameters.txt
#                summary.txt
#                (msScans.txt)
#                (experimentalDesign.txt)
#
# arguments
#    title        - character, the title for the pdf. Avois special character since this string is also used
#                   to generate file names
#    min.peptides - numeric, minimal number of peptides that have to be assigned to a protein
#                   group in each experiment to treat it as detected
#    mc           - logical, if true a barchart depicting the number of missed cleavages is plotted
#    separation   - logical, checks the separation, e.g. off gel, based on the 'Experiment'
#    rf           - logical, barchart of identified peptides per raw file
#    ma           - logical, density of calibrated mass error based on evidences
#    pep          - logical, distribution of PEP values of several tables
#    ol           - logical, if true the overlap between the specified experiments (experimentalDesign-file)
#                   is calculated
#    norm.ratio   - logical, if true the normalized ratios are used for the 'ratio' vs. 'intenesity' figure
#    msms         - logical, if false the 'msms.txt' is not loaded
#                 - sometimes necessary when working with 32bit systems
#    ic           - logical, if true the function tries to estimate an incorporation rate H/L
#    clean.up     - logical, if true all temporal files will be deleted after the pdf is produced
#    rm.mod       - character, specifies the modification that should not be included to calculate overlaps, ....
#    psig         - numeric, p-value indicating significantly regulated pg/ph
#    padjust      - character, method to adjust p-values (pg / ph)
#
#
# changed: 20090508  - implementation
#          200907    - added mass accuracy
#          20090904  - column 'Peptide.Counts' changed -> updated
#          20091013  - added parameter 'separation'
#          20091202  - added parameter 'rf', 'ma', 'pep', 'min.peptides'
#                    - added several slides, e.g. pep distribution, FDRs, ...
#          20091203  - some information about a SILAC experiment áre returned
#          20091204  - added support for version 1.14.3
#          20100121  - modified peptides table is now splitted if there are too many entries
#          20100211  - added some phospho stuff
#          20100212  - parameters phLscore, phPEP
#          20100408  - site redundancy
#                    - parameter psig, padjust
#          20100420  - bugfix: if the order of rawfiles in experimental design and the summary.txt were
#                           not the same, the raw file plot was coloured incorrectly
#          20100421  - renamed to 'summaryPDF'
#          20100621  - silac.type is set to FALSE if MQ before 1.0.13.x is used
#          20100915  - switched to function sigRatioPlot.ph also for protein groups
#          20101021  - added some new slides:
#                        - overlap protein groups (based on unique evidences)
#                        - peptides/protein groups per experiment
#                        - phosphosites (nr) per raw file
#                    - same x-scale in all ratio vs. intensity plots
#          201011    - # nr phospho sites per rawfile
#          20101124  - distribution of missed cleavages
#          20101129  - fixed bug in MC plot
#          201101    - overlap of peptide/protein groups split into light/medium/heavy if SILAC
#          20110310  - parameter 'ic'
#                    - incorporation check ( only one experiment/file at the moment)
#          20110318  - bugfix: ratio vs intensity of ph-sites - significant sites are reported again...
#          20110331  - same order of raw files in barplots of detected peptides and sites
#                    - added support for all kinds of modification tables
#          20110401  - parameter 'norm.ratio' now also affects sites
#                    - changed default parameters for 'padjust' and 'psig' to 'none' and '0.01'
#          20110404  - added parameter 'clean.up': removes all temporary files
#                    - removed parameter 'phospho'
#          20110405  - support for IC with multiple experiments
#          20110629  - bugfix: number of sites per raw file is now correct ( normal table, without CON/REV )
#          20111010  - added page numbers
#          20111118  - parameter 'low.level.qc'
#                    - additional slides on
#                             - MS and MSMS scans
#                             - ion injection times
#                             - cylce times
#                      per experiment and raw file
#          20111122  - added slide on SILAC mixing error
#                    - 'experimentalDesign' file is copied in the 'txt' folder automatically
#          20120208  - added slide on number of exclusively detected modification sites per raw file
#          20120210  - added slide on percent identified MS/MS per raw file
#                    - SILAC mixing error is now reported on raw scale (not log...)
#          20120213  - bugfix: cycle times are reported correctly if no experimetal design file was specified
#          20120329  - column names of tables are fixed, e.g.
#                      'Raw.file' -> 'Raw.File'
#          20120621  - mass deviations of fragment ions (png)
#                    - PEP boxplot as .'png'
#                    - msms table: only used columns are imported
#          20120625  - compatibility with 1.2.7.5, e.g. case of column names doesn't matter anymore
#                    - modified peptides table: only certain columns are imported
#          20120706  - sites with loc prob of zero are filtered out in the non-redundnat table
#          20120709  - several minor bugfixes
#          20120710  - fragment mass error: ppm vs. Da
#                    - for versions without having the ppm error in the tables, it is calculated
#          20131210  - compatible with 1.4.1.2
#                    - version check for loading modificationSpecificPeptides had to be adapted
#          20140822  - slide on enrichment efficiency of variable modifications
#          20140922  - compatible with 1.5.1.0 (no PEP in protein groups table)
#          20150112  - added protein sequence coverage
#
#
#############################################################################################################
summaryPDF <- function( title="test", min.peptides=1, mc=T, separation=F, rf=T, ma=T, pep=T, ol=T, norm.ratio=T, msms=T, phLscore=.75, phPEP=0.01, low.level.qc=T, load.files=T, psig=0.01, padjust="none", ic=F, clean.up=T, rm.mod = "Oxidation" )
{

    ##############################
    # current version
    ##############################
    vers = "0.9.6"

    #############################
    ## legend paramaters for  raw file plots
    leg.cex=1
    leg.row=4


    ##############################
    # start time
    ##############################
    time.start <- Sys.time()

    dir="."
    outDir="./summary"


    ########################################################
    # store the current working directory
    #########################################################
    d <- getwd()
    setwd(dir)

    ########################################################
    # files in the folder
    ########################################################
    filesInFolder <- list.files()

    ########################################################
    # - check wether the current directory contains any MQ files
    # - it tests for 'proteinGroups.txt'
    ########################################################
    check.folder <- length( grep('^proteinGroups.txt$', filesInFolder))

    while( check.folder == 0 ){

        dd <- choose.dir()
        if(is.na(dd)){
            setwd(d)

            stop("canceled\n")
        }

        setwd(dd)
        filesInFolder <- list.files()
        check.folder <- length( grep('^proteinGroups.txt$', filesInFolder))
    }


    ########################################################
    # get all site tables (variable modifications)
    ########################################################
    siteTabs <- filesInFolder[ grep(".*Sites\\.txt$", filesInFolder)  ]

    if(length(siteTabs) > 0){

        # do not consider Oxidation
        if(!is.null(rm.mod)){
            for(m in rm.mod){
                mod.idx <- grep(paste("^", m, sep=""), siteTabs)
                if(length(mod.idx) == 1 )
                    siteTabs <- siteTabs[-mod.idx]
            }
        rm(mod.idx)
        }

        names(siteTabs) <- sub(" .*", "", sub("Sites\\.txt", "", siteTabs))

        # name of the modification like in the header, e.g. Phospho..STY.
        mod.header.names <-  gsub("( |\\(|\\))", ".", sub("Sites\\.txt$","", siteTabs))
        names(mod.header.names) <- names(siteTabs)

    }


    #####################################
    #  paramters table
    #####################################
    param <- read.delim("./parameters.txt")


    # determine the number of pages for the list of paramaters,
    # each page holds 20 parameters
    nParamPages <- ceiling(dim(param)[1] / 20)

    ######################################
    #  determine MaxQuant version !
    ######################################
    MQversion <- as.character(param[ grep("^Version", param[, 1]) , 2])

    # representation as numerical value
    MQversion.num <- gsub("\\.", "", MQversion)
    # in order to asure that 1.13.9 > 1.12.35, i.e. 1139 > 11235
    if(nchar(MQversion.num) < 5) MQversion.num <- as.numeric(MQversion.num)*10

    rm(param)

    #######################################################
    # check whether the 'experimentalDesign' file is located
    #   in the current folder
    # - for newer MQ versions, this file is localated on folder
    #   up the hierarchy
    #######################################################
    if(MQversion.num > 12000){

        # if the file is not in the current folder
        if( length( grep("^experimentalDesignTemplate", filesInFolder)  ) == 0 ){

            # check whether the file can be found on folder above
            if(length( grep("^experimentalDesignTemplate", dir(".."))  ) > 0){
                dummy <- file.copy( paste( "../", dir("..")[ grep("^experimentalDesignTemplate", dir("..")) ], sep=""), "." )
            } else {
                warning("\nCannot find the 'experimentalDesign'-file!!\nTrying to proceed without...\n")
            }
        }

    }

    ########################################################
    # check whether there are experiments defined
    ########################################################
    ed.filename <- dir()[ grep("^experimentalDesignTemplate", dir())][1]
    experiments = NULL

    if(!is.na(ed.filename)){
        ed <- read.delim(ed.filename)

        if(sum( unlist(lapply( ed[, "Experiment"], is.na  ))   ) == 0  ){
            experiments <- unique( as.character(ed[, "Experiment"]) )

          #####################################################
          # fix the experiments written in the experimental design
          #  1. R converts special characters like '-' and '+' (but not the underscore '_')
          #     to dots -> in order to access the table columns
          #     one has to use 'experiments.dot'
          # 2. Latex doesn't like underscores. Thus for printing the experiment names
          #    to the tex file all '_' are replaced by '-'
          #####################################################

          # replace any special characters by "." -> in order to access column names
          experiments.dot <- gsub("-|\\+", "\\.", experiments)
          names(experiments.dot) <- experiments

          experiments.tex <- gsub("_", "-", experiments)
          names(experiments.tex) <- experiments

         }
    }

    ######################################
    #  check whether its SILAC or not
    ######################################

    # import protein groups header and check whether a column start with 'Ratio'
    # if so, its SILAC ...
    pg.head <- scan("./proteinGroups.txt", nlines=1, sep="\t", what = "character", quiet=T)
    quant = ifelse( length( grep("^Ratio.H.L", pg.head, ignore.case=T) ) > 0, T, F  )

    ###########################################
    # check whether its double or triple SILAC
    ###########################################
    if( quant){

        if( length( grep("^Ratio.M.L", pg.head, ignore.case=T) ) > 0){
            silac.type = "triplets"
        } else {
            silac.type = "doublets"
        }
    }
    else {
        silac.type = "singlets"
    }

    ########################################################
    # incorporation check: only if SILAC doublets
    ########################################################
    if(ic){
        if(silac.type != "doublets")
            ic = F
        if(ic)
            norm.ratio=F
    }
    ########################################################
    # check fragmentation type(s)
    ########################################################
    msms.head <-  make.names(scan("./msms.txt", nlines=1, sep="\t", what="character", quiet=T))

    frag.idx <- grep("^Fragmentation", msms.head, ignore.case=T)
    rev.idx.tmp <- grep("^Reverse$", msms.head, ignore.case=T)

    if(!is.null(frag.idx)){
       colclass.tmp <- rep("NULL", length(msms.head))
       colclass.tmp[c( frag.idx, rev.idx.tmp)] <- "character"
        # only read the column containing the type of fragmentation
        frag.types <- read.delim("./msms.txt", sep="\t", colClasses=colclass.tmp, header=T, comment.char="" )
        frag.types <- rmConRev(frag.types, con=F, q=F)
        frag.types <- table(frag.types[, grep("^Fragmentation", colnames(frag.types), ignore.case=F)])

   } else {
        frag.types=0
        names(frag.types)="Mmh!! Could not determine the type of fragmentation...\n"
   }

    #####################################
    # check whether LFQ
    #####################################
    lfq=F
    #if( (length(grep("^LFQ", pg.head)) > 0) |  (length(grep("^Norm.Intensity", pg.head)) > 0) )
    #    lfq=T


    ###############################################################################
    #
    #                         generate the Sweave file
    #
    ###############################################################################

    # latex header
    SweaveFile <- paste( "\\documentclass{beamer}
      %\\usepackage{beamerthemesingapore}
      \\usepackage{beamerthememarburg}
      %\\usepackage{beamerthemeberkeley}
      %\\usepackage{beamerthemebergen}
      %\\usepackage{beamerthemeannarbor}
      %\\usepackage{beamerthemepaloalto}


      %\\usepackage{longtable}
      %\\usepackage{pdfpages}

      %  page numbers
      %\\setbeamertemplate{footline}[frame number]

      % remove navigation panel
      \\setbeamertemplate{navigation symbols}{}

      % prefix string for figures
      \\SweaveOpts{prefix.string=pic/plot}

      %###################################################
      \\setbeamertemplate{footline}
      {%
       \\leavevmode%
       \\hbox{%
       \\begin{beamercolorbox}[wd=.333333\\paperwidth,ht=2.25ex,dp=1ex,center]{author in head/foot}%
       \\usebeamerfont{institute in head/foot}\\insertshortinstitute%~~(\\insertshortinstitute)
       \\end{beamercolorbox}%
       \\begin{beamercolorbox}[wd=.333333\\paperwidth,ht=2.25ex,dp=1ex,center]{title in head/foot}%
       \\usebeamerfont{title in head/foot}\\insertshorttitle
       \\end{beamercolorbox}%
       \\begin{beamercolorbox}[wd=.333333\\paperwidth,ht=2.25ex,dp=1ex,right]{date in head/foot}%
       \\usebeamerfont{date in head/foot}\\insertshortdate{}\\hspace*{2em}
       \\insertframenumber{} / \\inserttotalframenumber\\hspace*{2ex}
       \\end{beamercolorbox}}%
       \\vskip0pt%
      }

      %########################################################
      %
      %        preambel
      %
      %########################################################
      \\title[]{", title," }
      \\institute[MQparse v",vers,"]{}
      %\\institute{Proteome Center Tuebingen}


      \\begin{document}
      \\begin{frame}
       \\titlepage
      \\end{frame}
      \n\n", sep="")


    # load some packages
    SweaveFile <- paste( SweaveFile, "<<load_packages, echo=F>>=\n
         require(xtable)\n
         require(gplots)\n
         require(car)\n
         require(xtable)\n
    \n@\n", sep="")

    # initialize workspace
    if(load.files){
    SweaveFile <- paste(SweaveFile, "<<init, echo=F>>=\n
         norm.ratio=",norm.ratio,"\n
         # setwd(\"..\")\n@
         \n\n", sep="")


    #####################################################################################################
    #
    #                                  import the MaxQuant tables
    #
    # - determine some general numbers
    #
    #####################################################################################################

    ########################################
    #            protein groups
    ########################################
    SweaveFile <- paste(SweaveFile, "<<load_proteinGroups, echo=F>>=

      proteinGroups.head <-  scan('../proteinGroups.txt', nlines=1, sep='\t', what = 'character', quiet=T)
      proteinGroups=read.delim('../proteinGroups.txt', row.names=grep('^id$', proteinGroups.head), stringsAsFactors=F)

      ## compatibility to MQ 1.5.0.30
      colnames(proteinGroups) <- sub('Potential.contaminant', 'Contaminant', colnames(proteinGroups))

      nPG <- nrow(proteinGroups)

      # contaminants (proteins)
      nCon <- sum(proteinGroups[, grep('^Contaminant$', colnames(proteinGroups), ignore.case=T)] == '+')
      if(is.na(nCon))
             nCon = 0

      # reverse hits (proteins)
      nRev <- sum(proteinGroups[, grep('^Reverse$', colnames(proteinGroups), ignore.case=T)] == '+')
      if(is.na(nRev))
             nRev = 0

      # only identified by site
      nPgOnlyBySite <- sum(proteinGroups[, grep('^Only.identified.by.site$$', colnames(proteinGroups), ignore.case=T)] == '+')
      if(is.na(nPgOnlyBySite))
                  nPgOnlyBySite = 0

      ########################
      # remove CON/REV
      ########################
      proteinGroups <- rmConRev( proteinGroups, rev.col=grep('^Reverse$', colnames(proteinGroups), ignore.case=T), con.col=grep('^Contaminant$', colnames(proteinGroups), ignore.case=T), oib=colnames(proteinGroups)[ grep('^Only.identified.by.site$', colnames(proteinGroups), ignore.case=T) ], q=F  )

    \n@\n",sep="")

    ##############################################
    #             parameters
    ##############################################
    SweaveFile <- paste(SweaveFile, "<<load_parameters, echo=F>>=
      param <- read.delim('../parameters.txt')\n@\n",sep="")

    ##############################################
    #                evidence
    ##############################################
    SweaveFile <- paste(SweaveFile, "<<load_evidence, echo=F>>=
      evidence.head <-  scan('../evidence.txt', nlines=1, sep='\t', what = 'character', quiet=T)
      evidence <- read.delim('../evidence.txt', row.names=grep('^id$', evidence.head), stringsAsFactors=F)

      ## compatibility to MQ 1.5.0.30
      colnames(evidence) <- sub('Potential.contaminant', 'Contaminant', colnames(evidence))


      nEv <- nrow(evidence)
      nEvRev <- sum( evidence[, grep('^Reverse$', colnames(evidence), ignore.case=T )] == '+')
      if(is.na(nEvRev))
            nEvRev = 0

    ##########################
    #  remove Con/REV
    ##########################
    evidence <- rmConRev(evidence, rev.col=grep('^Reverse$', colnames(evidence), ignore.case=T), con.col=grep('^Contaminant$', colnames(evidence), ignore.case=T))


    ##########################
    # used columns
    ##########################

    \n@\n",sep="")

    ###############################################
    #            modified peptides
    ###############################################
    if( length( grep("^1.0.14", MQversion) ) > 0 | length( grep("^[1-9]\\.[1-9]", MQversion) ) > 0){
        SweaveFile <- paste(SweaveFile, "<<load_modifiedPeptides, echo=F>>=

                modifiedPep.head <-  scan('../modificationSpecificPeptides.txt', nlines=1, sep='\t', what = 'character', quiet=T)

                modifiedPep.col.idx <- match(c('id',
                                        'PEP',
                                        'Modifications',
                                        ifelse(length(grep('Potential.contaminant', modifiedPep.head)) > 0,'Potential.contaminant', 'Contaminant'),
                                        'Reverse'), modifiedPep.head)

                modifiedPep.col.classes <- rep('NULL', length(modifiedPep.head))
                modifiedPep.col.classes[modifiedPep.col.idx] <- 'character'

                modifiedPep <- read.delim('../modificationSpecificPeptides.txt', sep='\t', colClasses=modifiedPep.col.classes, header=T, comment.char='', row.names='id')

                ## compatibility to MQ 1.5.0.30
                colnames(modifiedPep) <- sub('Potential.contaminant', 'Contaminant', colnames(modifiedPep))


                nModPep <- nrow(modifiedPep)
                nModPepRev <- sum( modifiedPep[, grep('^Reverse$', colnames(modifiedPep), ignore.case=T)] == '+')
                if(is.na(nModPepRev))
                      nModPepRev=0

                ###########################
                # remove CON/REV
                ###########################
                modifiedPep <- rmConRev(modifiedPep, rev.col=grep('^Reverse$', colnames(modifiedPep), ignore.case=T), con.col=grep('^Contaminant$', colnames(modifiedPep), ignore.case=T))

                \n@\n",sep="")

    } else{
        SweaveFile <- paste(SweaveFile, "<<load_modifiedPeptides, echo=F>>=
                modifiedPep.head <-  scan('../modifiedPeptides.txt', nlines=1, sep='\t', what = 'character', quiet=T)

                modifiedPep.col.idx <- match(c('id',
                                        'PEP',
                                        'Modifications',
                                         ifelse(length(grep('Potential.contaminant', modifiedPep.head)) > 0,'Potential.contaminant', 'Contaminant'),
                                        'Reverse'),  modifiedPep.head)

                modifiedPep.col.classes <- rep('NULL', length(modifiedPep.head))
                modifiedPep.col.classes[modifiedPep.col.idx] <- 'character'

                modifiedPep <- read.delim('../modifiedPeptides.txt', sep='\t', colClasses=modifiedPep.col.classes, header=T, comment.char='', row.names='id')


                ## compatibility to MQ 1.5.0.30
                colnames(modifiedPep) <- sub('Potential.contaminant', 'Contaminant', colnames(modifiedPep))


                nModPep <- nrow(modifiedPep)
                nModPepRev <- sum( modifiedPep[, grep('^Reverse$', colnames(modifiedPep), ignore.case=T)] == '+')
                if(is.na(nModPepRev))
                      nModPepRev=0

                ###########################
                # remove CON/REV
                ###########################
                modifiedPep <- rmConRev(modifiedPep, rev.col=grep('^Reverse$', colnames(modifiedPep), ignore.case=T), con.col=grep('^Contaminant$', colnames(modifiedPep), ignore.case=T))

               \n@\n",sep="")
    }


    ################################################
    #                     msms
    ################################################
    if(msms){
        SweaveFile <- paste(SweaveFile, "<<load_msms, echo=F>>=
                # identify columsn to be loaded
                msms.head <- make.names( scan('../msms.txt', nlines=1, sep='\t', what='character', quiet=T))
                msms.firstline <- read.delim('../msms.txt', row.names=1, nrow=1, stringsAsFactors=F)


                # check whether there are fragment mass deviations in the table
                me.col.idx <- msms.head[grep('^Mass.Deviations', msms.head, ignore.case=T)]


                msms.col.idx <- match(c('id',
                                        'PEP',
                                        'm.z',
                                        me.col.idx,
                                        'Masses',
                                        'Fragmentation',
                                        'Reverse'), msms.head)

                msms.col.classes <- rep('NULL', length(msms.head))
                msms.col.classes[msms.col.idx] <- 'character'

                msms <- read.delim('../msms.txt', sep='\t', colClasses=msms.col.classes, header=T, comment.char='', row.names=grep('^id$', msms.head, ignore.case=T), stringsAsFactors=F)

                nMsms <- nrow(msms)
                nMsmsRev <- sum( msms[, grep( '^Reverse$', colnames(msms), ignore.case=T )] == '+')
                if(is.na(nMsmsRev))
                      nMsmsRev=0

                ###########################
                #  remove CON/REV
                ###########################
                msms <- rmConRev(msms, con=F, rev.col=grep('^Reverse$', colnames(msms), ignore.case=T))

                # extract PEP
                msmsPEP <- as.numeric(msms[, grep('^PEP$', colnames(msms), ignore.case=T)] )
                # mass to charge
                msms.mz <- as.numeric( msms[, grep('^m.z$', colnames(msms), ignore.case=T)])
                names(msmsPEP) <- names(msms.mz) <- rownames(msms)

                ##############################################
                # if Da and ppm are present
                ##############################################
                if(length(me.col.idx) == 2 ){
                     # msms mass error in ppm
                     msms.me.ppm <- msms[ , me.col.idx[ grep('ppm', me.col.idx) ] ]
                     # msms error in Da
                     msms.me.Da <- msms[ , me.col.idx[ grep('Da', me.col.idx) ] ]
                     names(msms.me.ppm) <- names(msms.me.Da) <- rownames(msms)
                }
                #############################################
                # if only Da is present
                #############################################
                if( (length(me.col.idx) == 1) & (length( grep('CID', msms.firstline[grep('Fragmentation', names(msms.firstline), ignore.case=T)])) == 1) ){

                     # msms error in Da
                     msms.me.Da <- msms[ , me.col.idx ]
                     names(msms.me.Da) <- rownames(msms)

                     # fragment masses
                     msms.masses <- msms[, grep('^Masses$', colnames(msms), ignore.case=T)]
                     names(msms.masses) <- row.names(msms)



                    ## msms.me.ppm <- unlist(lapply( names(msms.me.Da), function(x)  paste( as.numeric( unlist( strsplit( as.character(msms.me.Da[x]), ';'))) / (as.numeric( unlist( strsplit( as.character(msms.masses[x]), ';')))*1e-6) , collapse=';')   ))

                }
                ###############################################
                # fragment mass error
                ###############################################
                if(length(me.col.idx) < 1 | length(me.col.idx) > 2){
                     msms.me.Da <- rep(0, nrow(msms))
                     msms.me.ppm <- msms.me.Da
                     names(msms.me.Da) <- names(msms.me.ppm) <- rownames(msms)
                }

                # fragmentation type
                msms.frag <- msms[ , grep( '^Fragmentation$', colnames(msms), ignore.case=T)]
                names(msms.frag) <- rownames(msms)

                # fragment masses: needed to calculate ppm error, if neccessary
                msms.masses <- msms[, grep('^Masses$', colnames(msms), ignore.case=T)]
                names(msms.masses) <- rownames(msms)

                rm(msms)
\n@\n",sep="")

    } else{
        SweaveFile <- paste(SweaveFile, "<<load_msms, echo=F>>=
                nMsms <- nMsmsRev <- msms.mz <- 0

                \n@\n",sep="")
    }

    ###################################################
    #                    peptides
    ###################################################
    SweaveFile <- paste(SweaveFile, "<<load_peptides, echo=F>>=

                peptides.head <- make.names( scan('../peptides.txt', nlines=1, sep='\t', what='character', quiet=T))
                peptides <- read.delim('../peptides.txt', row.names=grep('^id$', peptides.head), stringsAsFactors=F)


                ## compatibility to MQ 1.5.0.30
                colnames(peptides) <- sub('Potential.contaminant', 'Contaminant', colnames(peptides))


                # contaminants
                nPepCon = NA
                if(length(grep('^Contaminant$', colnames(peptides), ignore.case=T) > 0))
                        nPepCon <- sum(peptides[, grep('^Contaminant$', colnames(peptides), ignore.case=T)] == '+')
                if(is.na(nPepCon))
                     nPepCon = 0
                # reverse
                nPepRev <- sum(peptides[, grep( '^Reverse$', colnames(peptides), ignore.case=T ) ] == '+')
                if(is.na(nPepRev))
                     nPepRev = 0

                ###############################
                # remove CON/REV
                ###############################
                peptides <- rmConRev(peptides, rev.col=grep('^Reverse$', colnames(peptides), ignore.case=T), con.col=grep('^Contaminant$', colnames(peptides), ignore.case=T))


    \n@\n",sep="")

    ###################################################
    #                    site tables
    ###################################################
    for(s in names(siteTabs)){
        SweaveFile <- paste(SweaveFile, "<<load_",s,", echo=F>>=
                 \n",s,"Sites.head <-  make.names( scan('../", siteTabs[s],"', nlines=1, sep='\t', what='character', quiet=T))
                 \n",s,"Sites <- read.delim('../",siteTabs[s],"', row.names=grep('^id$', ",s,"Sites.head), stringsAsFactors=F)


                ## compatibility to MQ 1.5.0.30
                colnames(",s,"Sites) <- sub('Potential.contaminant', 'Contaminant', colnames(",s,"Sites))

                 #colnames(", s,"Sites) <- unlist(lapply(colnames(",s,"Sites), function(x) paste(capwords(unlist(strsplit(x, '\\\\.'))), collapse='.') ))
                 \n@\n",sep="")
    }

    ###################################################
    #                    summary
    ###################################################
    SweaveFile <- paste(SweaveFile, "<<load_summary, echo=F>>=
                 summaryTxt <- read.delim('../summary.txt')
                 #colnames(summaryTxt) <- fixColumnNames(colnames(summaryTxt))
    \n@\n",sep="")

    ###################################################
    #               experimental design
    ###################################################
    if(!is.na(ed.filename)){
         SweaveFile <- paste(SweaveFile, "<<load_experimentalDesign, echo=F>>=
                 ed <- read.delim('../",ed.filename,"')\n@\n",sep="")

    } else{ # build  an own experimentalDesign table

          SweaveFile <- paste(SweaveFile, "<<build_experimentalDesign, echo=F>>=
                   rf <- unique(as.character(evidence[, grep('^Raw.File$', colnames(evidence), ignore.case=T)]))
                   exp.dummy <- rep('', length(rf))
                   ed <- data.frame(Name=rf, Experiment=exp.dummy)\n@\n",sep="")

    }

    if(low.level.qc){
        ###################################################
        #                    msScans
        ###################################################
        SweaveFile <- paste(SweaveFile, "<<load_msScans, echo=F>>=
                   msScans <- read.delim('../msScans.txt', stringsAsFactors=F)
                   #colnames(msScans) <- unlist(lapply(colnames(msScans), function(x) paste(capwords(unlist(strsplit(x, '\\\\.'))), collapse='.') ))

                   # fill the first column with names of the raw files
                   rf.idx <- which( nchar(msScans[ , grep('^Raw.File$', colnames(msScans), ignore.case=T) ] ) > 0)
                   if(length(rf.idx) < dim(msScans)[1]){
                       rf.idx <- c(rf.idx, dim(msScans)[1] + 1)

                       for(tmp in 1:(length(rf.idx)-1) ){
                            # get the file name
                            rf.tmp <- msScans[ rf.idx[tmp], grep('^Raw.File$', colnames(msScans), ignore.case=T)]

                            msScans[(rf.idx[tmp]+1):(rf.idx[tmp+1]-1) , grep('^Raw.File$', colnames(msScans), ignore.case=T)] <- msScans[rf.idx[tmp], grep('^Raw.File$', colnames(msScans), ignore.case=T)]
                       }
                   }
        \n@\n",sep="")


        ######################################################
        #             allPeptides
        ######################################################
        # identify columsn to be loaded
        #ap.head <- make.names( scan('../allPeptides.txt', nlines=1, sep='\t', what='character', quiet=T))
        #ap.firstline <- read.delim('../allPeptides.txt', row.names=1, nrow=1, stringsAsFactors=F)


        # columns that are loaded
        #ap.rf.col.idx <- ap.head[ grep('^Raw.File', ap.head, ignore.case=T) ]
        #ap.int.col.idx <- ap.head[ grep('^Intensity', ap.head, ignore.case=T) ]


        #ap.col.idx <- match(c(ap.rf.col.idx,
        #                      ap.int.col.idx,
        #                      ), ap.head)

        #ap.col.classes <- rep('NULL', length(ap.head))
        #ap.col.classes[ap.col.idx] <- 'character'

        #ap <- read.delim('../allPeptides.txt', sep='\t', colClasses=msms.col.classes, header=T, comment.char='', row.names=grep('^id$', msms.head, ignore.case=T), stringsAsFactors=F)


    } # end if low.level.qc


    } # end if load.files

    ###################################################################
    #
    #                  some general numbers
    #
    ###################################################################
    SweaveFile <- paste(SweaveFile, "<<someNumbers, echo=F>>=\n

      #############################################
      # check whether there are experiments defined
      ##############################################
      experiments = NULL
      if(sum(grepl('^Experiment$', colnames(evidence), ignore.case=T) ) == 1){
          experiments <- unique( as.character(evidence[, grep('^Experiment$', colnames(evidence), ignore.case=T)]) )

          experiments.dot <- gsub(\"-|\\\\+\", \".\", experiments)
          names(experiments.dot) <- experiments
      }

      #########################
      # protein groups
      #########################

      ####################
      # identified by site
      #nPgOnlyBySite = 0
      # check whether the column is present ( depends on MQ version  )
      #if( length( grep( 'Only.identified.by.site', colnames(proteinGroups) ,ignore.case=T) ) > 0){
      #    nPgOnlyBySite = sum( proteinGroups[, grep( 'Only.identified.by.site', colnames(proteinGroups) ,ignore.case=T)] == '+'  )
      #    if(is.na(nPgOnlyBySite))
      #            nPgOnlyBySite = 0
      #}

      #####################
      # single peptide pg
      # find the column name for the peptide counts. Since version .13. some additional columns have been added
      # for unique peptide counts and so on...
      pepCountColname <- colnames(proteinGroups)[grep('^Peptide.Counts', colnames(proteinGroups), ignore.case=T) ]
      if(length(pepCountColname) > 1){
          pepCountColname <- pepCountColname[ grep('^Peptide.Counts..all.$', pepCountColname, ignore.case=T) ]
      }
      nSinglePepProt <- sum( as.character(proteinGroups[, pepCountColname ]) == '1'  )
      # check if there are pg encompassing several proteins that all where identified by a single peptide,
      # e.g. 'Peptide.Counts' = 1;1;1
      tmp <- unlist(lapply( as.character(proteinGroups[, pepCountColname]), function(x){ xx <-  unlist( strsplit(x,';')  ); if(length( xx  ) > 1  ){ ifelse( sum(as.numeric(xx)) == length(xx), return(1),  return(0)   )   };return(0)  } ) )
      nSinglePepProt <- nSinglePepProt + sum(tmp)


      #########################
      #   peptides
      #########################
      # non-redundant sequences
      nPepNR <- length( unique( as.character(evidence[, grep('^Sequence$', colnames(evidence), ignore.case=T )])) )


      ###############################
      # unique peptides: depends again on the version...
      ################################
      nPepUn=NA
      #if( 'Unique' %in% colnames(peptides)  ){
      if(length( grep('^Unique$', colnames(peptides), ignore.case=T )  )>0){
         nPepUn <- sum(peptides[, grep('^Unique$', colnames(peptides), ignore.case=T ) ] == 'yes')
      }
      #if( 'Unique..Groups.' %in% colnames(peptides)) {
      if(length( grep('^Unique..Groups.$', colnames(peptides), ignore.case=T )  )>0){
         nPepUn <- sum(peptides[, grep('^Unique..Groups.$', colnames(peptides), ignore.case=T )] == 'yes')
      }
      #########################
      # modified peptides
      #########################
      # number of modified peptides
      nModPep.modified <- sum( modifiedPep[, grep('^Modifications$', colnames(modifiedPep), ignore.case=T)] != 'Unmodified', na.rm=T )
      # the modifications
      modsNumb <- data.frame(table(modifiedPep[, grep('^Modifications$', colnames(modifiedPep), ignore.case=T)]), row.names=1)


      #########################
      #          FDRs
      #########################
      fdr.table <- matrix(0, ncol=5, nrow=5, dimnames=list( c('protein groups', 'peptides', 'modified peptides', 'evidences', 'msms'), c('total', 'forward', 'reverse', 'FDR1', 'FDR2') ))
      fdr.table['protein groups', ] <- c( nPG, nPG-nRev, nRev, round((nRev/(nPG-nRev))*100,2), round((2*nRev/(nPG))*100,2)  )
      fdr.table['peptides', ] <- c( nPepNR, nPepNR-nPepRev, nPepRev, round((nPepRev/(nPepNR-nPepRev))*100,2), round((2*nPepRev/(nPepNR))*100,2)  )
      fdr.table['modified peptides', ] <- c( nModPep, nModPep-nModPepRev, nModPepRev, round((nModPepRev/(nModPep-nModPepRev))*100,2), round((2*nModPepRev/(nModPep))*100,2)  )
      fdr.table['evidences', ] <- c( nEv, nEv-nEvRev, nEvRev, round((nEvRev/(nEv-nEvRev))*100,2),  round((2*nEvRev/(nEv))*100,2)  )
      fdr.table['msms', ] <- c( nMsms, nMsms-nMsmsRev, nMsmsRev, round((nMsmsRev/(nMsms-nMsmsRev))*100,2), round((2*nMsmsRev/(nMsms))*100,2)  )


      ######################################
      #          Raw Files
      ######################################
      rf <- as.character( ed[, 'Name']  )

@\n", sep="")

    ##########################################
    #           Site FDRs
    ##########################################
    if(length(siteTabs) > 0){

        for(s in names(siteTabs)){

            ################################################
            #              FDR table
            ##################################################
             SweaveFile <- paste(SweaveFile, "\n<<",s,"FDR, echo=F>>=\n
                         n",s,"Site = dim(",s,"Sites)[1]
                         n",s,"SiteRev <- length(grep( \"^\\\\+$\", ",s,"Sites[, 'Reverse'] ))
                         n",s,"SiteCon <- length(grep( \"^\\\\+$\", ",s,"Sites[, 'Contaminant'] ))

                         fdr.table <- rbind(fdr.table, c( n",s,"Site, n",s,"Site-n",s,"SiteRev, n",s,"SiteRev, round( (n",s,"SiteRev/(n",s,"Site-n",s,"SiteRev))*100,2),round( 2*(n",s,"SiteRev/(n",s,"Site))*100,2) )  )
                         rownames(fdr.table)[dim(fdr.table)[1]] <- \"",s," sites\"
                         \n\n@\n\n",sep="")

            ####################################################
            # filter the  sites:
            #  1. remove reverse hits and contaminants
            #  2. L > 0.75
            #  3. L > 0.75 & PEP < 0.01
            ####################################################
            SweaveFile <- paste(SweaveFile, "\n<<",s,"Filter, echo=F>>=\n
                         # remove reverse hits and contaminants\n
                         ",s,"Sites <- rmConRev(",s,"Sites)

                         #if(n",s,"SiteRev > 0) ",s,"Sites <- ",s,"Sites[ -grep( '^\\\\+$', ",s,"Sites[, 'Reverse'] ), ]
                         #if(length(grep('^\\\\+$', ",s,"Sites[, 'Contaminant']) ) > 0) ",s,"Sites <- ",s,"Sites[ -grep( \"^\\\\+$\", ",s,"Sites[, 'Contaminant']  ), ]
                         # get localized sites
                         ",s,"Sites.loc <- ",s,"Sites[ which(",s,"Sites[, grep('^Localization.Prob$', colnames(",s,"Sites ), ignore.case=T)] >= ",phLscore,"), ]

                         # localized sites + PEP filter
                         ",s,"Sites.loc.PEP.0.01 <- ",s,"Sites.loc[which(",s,"Sites.loc[, grep( '^PEP$', colnames(",s,"Sites), ignore.case=T)] <= ",phPEP," ), ]
                         \n\n@\n\n",sep="")

        }

    }

    ###############################################################################################################################
    #
    #                                         the first slide ...
    #
    ###############################################################################################################################

    ##########################################
    # identification: protein groups
    ##########################################
    SweaveFile <- paste(SweaveFile, "\\section{Identification}
         \\begin{frame}
         \\frametitle{Protein Groups}
          \\begin{description}
            \\item[Protein groups:] \\Sexpr{nPG}     \\
            \\item[Single peptide PG:] \\Sexpr{nSinglePepProt} (\\Sexpr{round((nSinglePepProt/nPG)*100, 3)} \\%)
            \\item[Sequence coverage (median):] \\Sexpr{median( as.numeric(proteinGroups[, grep(paste('^Sequence.coverage....$'), colnames(proteinGroups),ignore.case=T,value=T)] ))} \\% \\
            \\item[Only identified by site:] \\Sexpr{nPgOnlyBySite} (\\Sexpr{round((nPgOnlyBySite/nPG)*100, 3)} \\%)
            \\item[Contaminants:] \\Sexpr{nCon} (\\Sexpr{round((nCon/nPG)*100, 3)} \\%)   \\
            \\item[Reverse hits:] \\Sexpr{nRev}  \\
          \\end{description}
         \\end{frame}\n", sep="")

    ##########################################
    # identification: peptides
    ##########################################
        SweaveFile <- paste(SweaveFile, " \\begin{frame}
         \\frametitle{Peptides}
          \\begin{description}
            \\item[Identified evidences:] \\Sexpr{nEv}     \\\
            \\item[Reverse hits:] \\Sexpr{nEvRev}    \\\
          \\end{description}
          \\vspace{4mm}
          \\begin{description}
            \\item[Identified peptides (non-redundant):] \\Sexpr{nPepNR}  \\\
            \\item[Peptides unique to a PG:] \\Sexpr{ nPepUn} (\\Sexpr{round((nPepUn/nPepNR)*100, 3)} \\%) \\\
            \\item[Contaminants:] \\Sexpr{nPepCon} (\\Sexpr{round((nPepCon/nPepNR)*100, 3)} \\%) \\\
            \\item[Reverse hits:] \\Sexpr{nPepRev} \\\
          \\end{description}
          \\vspace{4mm}
         \\begin{description}
            \\item[Spectra submitted:] \\Sexpr{summaryTxt[dim(summaryTxt)[1], 'MS.MS.Submitted']} \\\\
            \\item[Spectra identified:] \\Sexpr{summaryTxt[dim(summaryTxt)[1], 'MS.MS.Identified'] } (\\Sexpr{round((summaryTxt[dim(summaryTxt)[1], 'MS.MS.Identified']/ summaryTxt[dim(summaryTxt)[1], 'MS.MS.Submitted'])*100, 2)} \\%) \\\\
         \\end{description}
        \\end{frame}\n", sep="")

    ##########################################
    # False discovery rates
    ##########################################
         SweaveFile <- paste(SweaveFile, "\\section{False Discovery Rates}\n\\begin{frame}\n
         \\frametitle{False Disovery Rates}\n
         False discovery rates are estimated as\n

         \\vspace{2mm}\n
         \\begin{center}
            $\\hat{FDR_1}= \\frac{N_{rev}}{N_{fwd}} \\cdot 100$ and $\\hat{FDR_2}= \\frac{2 \\cdot N_{rev}}{N_{fwd} + N_{rev}} \\cdot 100$\n
         \\end{center}

         \\vspace{2mm}\n", sep="")

         SweaveFile <- paste( SweaveFile, "<<fdrTab, echo=F, results=tex>>=
          print(xtable(fdr.table, digits=c(0,0,0,0,2,2)), size='small')\n@\n", sep="")

         SweaveFile <- paste( SweaveFile, "\\end{frame}\n", sep="")


    ##########################################
    # modified peptides
    ##########################################


            SweaveFile <- paste(SweaveFile, "\\begin{frame}
                          \\frametitle{Modified Peptides}\n
                             \\begin{description}
                               \\item[Number of modified peptides:] \\Sexpr{nModPep.modified}     \\\
                             \\end{description}\n",  sep="")

            SweaveFile <- paste(SweaveFile, "
                         \\begin{center}
                         \\begin{minipage}{2in}
                         \n<<modTab1, echo=F, results=tex>>=\n

                         # ',' as separator
                         if(length(grep(',', modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T)] )) > 0){
                              modPep.all <- unlist( strsplit( modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T)], ','))

                         } else  if(length(grep(';', modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T)] )) > 0) {
                              # ';' as separator: MQ 1.5.0.0
                              modPep.all <- unlist( strsplit( modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T)], ';'))

                         } else {
                               modPep.all <- modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T)]
                         }


                         modPep.collapsed <- sub('^. ', '', modPep.all)

                         modTab.tmp <- table(modPep.all)
                         modTab.tmp.collapsed <- table(modPep.collapsed)

                         print(xtable( modTab.tmp ), size='tiny', include.colnames=F)

                         \n@\n
                         \\end{minipage}
                         \n",sep="")


    SweaveFile <- paste(SweaveFile, "
                         \\begin{minipage}{2in}

                         \n<<modTab2, echo=F, results=tex>>=\n
                         print(xtable( modTab.tmp.collapsed ), size='tiny', include.colnames=F)
                         \n@\n
                         \\end{minipage}
                         \\end{center}
                         \n",sep="")




            SweaveFile <- paste( SweaveFile, "\\end{frame}\n", sep="")
    ##########################################
    #  missed cleavages
    ##########################################
    if(mc){

        SweaveFile <- paste(SweaveFile, "\\section{Missed Cleavages}\n", sep="")

        # scan the header of peptides table
        peptides.head <- read.delim("peptides.txt", nrow=1, stringsAsFactors=F)
        missed.cleav <- colnames(peptides.head)[grep('^[M|m]issed.[C|c]leavages', colnames(peptides.head)) ]

        for(mc in 1:length(missed.cleav)){

            SweaveFile <- paste(SweaveFile, "\\begin{frame}
                          \\frametitle{Missed Cleavages}\n
                            \\begin{center}
                         \n<<MC",mc,", echo=F, fig=T, height=7, width=8>>=\n
                            mc.tab <- table(peptides[,'",missed.cleav[mc],"'])
                            mc.tab.vec <- as.vector(mc.tab)
                            names(mc.tab.vec) <- names(mc.tab)
                            fancyBarplot( mc.tab.vec, col='grey', xlab='", missed.cleav[mc] ,"')
                            legend('topright', legend=paste(names(mc.tab.vec), ' m.c.:', round(100*(mc.tab.vec/sum(mc.tab.vec)),1), ' %' ), bty='n', cex=1.5  )

                         \n@\n
                            \\end{center}
                         \n\\end{frame}\n",sep="")

        }
    }

    if(length(siteTabs) > 0){

        ##########################################################################
        #
        #                          modification sites
        #
        # - nr-sites
        # - export of nr-sites
        # - phospho site ratios
        # - distribution of localization probability
        # - boxplots PEP sites
        # - boxplots PEP unmodified vs. sites
        # - overlap: detected/quantified sites per experiment
        ###########################################################################
        for(s in names(siteTabs)){

            SweaveFile <- paste(SweaveFile, summaryPDF.Sites( s, what="general", mod.header=mod.header.names[s], phLscore=phLscore, phPEP=phPEP, quant=quant, experiments=experiments  ), sep="")

        }

    }

    #########################################
    # Raw Files and Experiments
    #########################################
    if(rf){
      SweaveFile <- paste(SweaveFile, "\\section{Raw Files}
         \\begin{frame}
         \\frametitle{Raw Files: Peptide Sequences}
            \\begin{description}
              \\item[Number of files] \\Sexpr{length(rf)}     \\\
            \\end{description}\n",sep="")


      ##############################################
      # barplot: identified peptides per raw file
      ##############################################
      SweaveFile <- paste( SweaveFile, "\\setkeys{Gin}{width=0.8\\textwidth}\n
            \\begin{center}\n<<rfPlot, echo=F, fig=T, height=10, width=15>>=
            # column name of raw files
            #rf.colum.name <- colnames(summaryTxt)[grep('^(r|R)aw\\.(F|f)ile', colnames(summaryTxt))]
            if(!is.null(experiments)){

                   cols.tmp <- rainbow(length(experiments))
                   names(cols.tmp) <- experiments

                   # remove rows that correspond to experiments, i.e. get only the numbers for each raw file
                   pepPerRF <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Peptide.Sequences.Identified$', colnames(summaryTxt), ignore.case=T)]

                   names(pepPerRF) <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]

                   pepPerRF <- pepPerRF[-c(length(pepPerRF))]

                   # determine the experiment for each rawfile
                   barCol = as.character(ed[, 'Experiment'])
                   names(barCol) <- as.character( ed[, 'Name'] )

                   for(e in experiments){
                         barCol[which(barCol == e )] <- cols.tmp[e]
                   }

                   # bring the vectors 'barCol' and 'pepPerRF' in the same order
                   pepPerRF <- pepPerRF[names(barCol)]
            }

            if(is.null(experiments)){
                   pepPerRF <- summaryTxt[, grep('^Peptide.Sequences.Identified$', colnames(summaryTxt), ignore.case=T) ]
                   names(pepPerRF) <- summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T) ]
                   barCol='grey'
                   # remove the 'total' row
                   pepPerRF <- pepPerRF[-c(length(pepPerRF))]
            }

            #################################################
            # order according to raw file name
            #################################################
            rf.order <- names(pepPerRF)[ order(names(pepPerRF))]
            pepPerRF <- pepPerRF[rf.order]
            if(!is.null(experiments))
                      barCol <- barCol[rf.order]

            par(mar = c(20, 8, 4, 2)  )
            fancyBarplot( as.vector(pepPerRF), main='Peptide sequences identified (summary.txt)', col=barCol, ylim=c(0,max(pepPerRF, na.rm=T)+max(pepPerRF, na.rm=T)*.3), cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names(pepPerRF), cex.main=3, ylab='count'  )\n
            grid(nx=NA, ny=NULL)
            # legend
            if(!is.null(experiments))
                 legend('top', legend=experiments, fill=cols.tmp, ncol=ceiling(length(experiments)/",leg.row,"), cex=",leg.cex,", bty='n'  )
                 ##legend('top', legend=experiments, fill=cols.tmp, ncol=ifelse(length(experiments) > 8, 3, 2), cex=1.5, bty=\"n\"  )\n
            par(mar=c(5,4,4,2))
         \n@\n\\end{center}", sep="")
      SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")



      ##############################################
      # barplot: percent identified spectra per raw file
      ##############################################

      SweaveFile <- paste(SweaveFile, "
        \\begin{frame}
        \\frametitle{Raw Files: Percent MS/MS Identified}
            \\begin{description}
              \\item[Number of files] \\Sexpr{length(rf)}     \\\
            \\end{description}\n",sep="")

      SweaveFile <- paste( SweaveFile, "\\setkeys{Gin}{width=0.8\\textwidth}\n
            \\begin{center}\n<<rfPlotpercent, echo=F, fig=T, height=10, width=15>>=

            if(!is.null(experiments)){

                   cols.tmp <- rainbow(length(experiments))
                   names(cols.tmp) <- experiments

                   # remove rows that correspond to experiments, i.e. get only the numbers for each raw file
                   #specIdentPerRF <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  'MS.MS.Identified....']
                   specIdentPerRF <- summaryTxt[ !( as.character( summaryTxt[ , grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^MS.MS.Identified\\\\.\\\\.\\\\.\\\\.$', colnames(summaryTxt), ignore.case=T)]
                   names(specIdentPerRF) <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]

                   specIdentPerRF <- specIdentPerRF[-c(length(specIdentPerRF))]

                   # determine the experiment for each rawfile
                   barCol = as.character( ed[ ,'Experiment'] )
                   names(barCol) <- as.character( ed[,'Name'] )

                   for(e in experiments){
                         barCol[which(barCol == e )] <- cols.tmp[e]
                   }

                   # bring the vectors 'barCol' and 'pepPerRF' in the same order
                   specIdentPerRF <- specIdentPerRF[names(barCol)]
            }

            if(is.null(experiments)){
                   specIdentPerRF <- summaryTxt[, grep('^MS.MS.Identified\\\\.\\\\.\\\\.\\\\.$', colnames(summaryTxt), ignore.case=T) ]
                   names(specIdentPerRF) <- summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T) ]
                   barCol='grey'
                   # remove the 'total' row
                   specIdentPerRF <- specIdentPerRF[-c(length(specIdentPerRF))]
            }

            #################################################
            # order according to raw file name
            #################################################
            specIdentPerRF <- specIdentPerRF[rf.order]
            if(!is.null(experiments))
                      barCol <- barCol[rf.order]

            par(mar = c(20, 4, 4, 2)  )
            fancyBarplot( as.vector(specIdentPerRF), main='MS/MS identified in percent (summary.txt)', col=barCol, ylim=c(0,max(specIdentPerRF, na.rm=T)+max(specIdentPerRF, na.rm=T)*.3), cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names(specIdentPerRF), cex.main=3, ylab='percent'  )\n
            grid(nx=NA, ny=NULL)
            # legend
            if(!is.null(experiments))
                 legend('top', legend=experiments, fill=cols.tmp, ncol=ceiling(length(experiments)/",leg.row,"), cex=",leg.cex,", bty='n'  )
                 ##legend('top', legend=experiments, fill=cols.tmp, ncol=ifelse(length(experiments) > 8, 3, 2), cex=1.5, bty=\"n\"  )\n
            par(mar=c(5,4,4,2))
         \n@\n\\end{center}", sep="")
      SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")


      ##############################################
      # barplot: average absolute mass error per raw file
      ##############################################

      SweaveFile <- paste(SweaveFile, "
        \\begin{frame}
        \\frametitle{Raw Files: Average Absolute Mass Error}
            \\begin{description}
              \\item[Number of files] \\Sexpr{length(rf)}     \\\
            \\end{description}\n",sep="")

      SweaveFile <- paste( SweaveFile, "\\setkeys{Gin}{width=0.8\\textwidth}\n
            \\begin{center}\n
            \n<<rfPlotAvME, echo=F, fig=T, height=10, width=15>>=

            if(!is.null(experiments)){

                   # remove rows that correspond to experiments, i.e. get only the numbers for each raw file
                   avMEPerRF <- summaryTxt[ !( as.character(summaryTxt[, grep( '^Raw.File$', colnames(summaryTxt), ignore.case=T )]) %in% experiments ),  grep('^Av..Absolute.Mass.Deviation($|..ppm.$)', colnames(summaryTxt), ignore.case=T)]

                   names(avMEPerRF) <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]
                   avMEPerRF <- avMEPerRF[-c(length(avMEPerRF))]

                   stMEPerRF <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Mass.Standard.Deviation($|..ppm.$)', colnames(summaryTxt), ignore.case=T)]
                   names(stMEPerRF) <- summaryTxt[ !( as.character(summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]) %in% experiments ),  grep('^Raw.File$', colnames(summaryTxt), ignore.case=T)]
                   stMEPerRF <- stMEPerRF[-c(length(stMEPerRF))]


            }

            if(is.null(experiments)){
                   # average absolute mass error
                   avMEPerRF <- summaryTxt[, grep( '^Av..Absolute.Mass.Deviation($|..ppm.$)', colnames(summaryTxt), ignore.case=T) ]
                   names(avMEPerRF) <- summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T ) ]
                   # remove the 'total' row
                   avMEPerRF <- avMEPerRF[-c(length(avMEPerRF))]

                   # standard deviation
                   stMEPerRF <- summaryTxt[, grep('^Mass.Standard.Deviation($|..ppm.)', colnames(summaryTxt), ignore.case=T ) ]
                   names(stMEPerRF) <- summaryTxt[, grep('^Raw.File$', colnames(summaryTxt), ignore.case=T) ]
                   # remove the 'total' row
                   stMEPerRF <- stMEPerRF[-c(length(stMEPerRF))]

            }

            #################################################
            # order according to raw file name
            #################################################
            avMEPerRF <-  avMEPerRF[rf.order]
            stMEPerRF <-  stMEPerRF[rf.order]

            barCol.tmp <- c('grey', 'black')

            par(mar = c(20, 4, 4, 2)  )
            fancyBarplot( rbind( as.vector( avMEPerRF),as.vector( stMEPerRF)), main='Mass error in ppm (summary.txt)',  col=barCol.tmp, ylim=c(0,max( c(avMEPerRF,stMEPerRF), na.rm=T)+max( c(avMEPerRF, stMEPerRF), na.rm=T)*.3), cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names( avMEPerRF), cex.main=3, ylab='ppm', add.numb=F  )\n
            grid(nx=NA, ny=NULL)
            # legend
            legend('top', legend=c('average absolute', 'standard deviation'), fill=barCol.tmp, cex=",leg.cex,", bty='n', ncol=2  )\n
            par(mar=c(5,4,4,2))
         \n@\n\\end{center}", sep="")
      SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

      #####################################################
      ##  Uncalibrated mass error
      #####################################################
      SweaveFile <- paste(SweaveFile, "
        \\begin{frame}
        \\frametitle{Raw Files: Uncalibrated Mass Error}
            \\begin{description}
              \\item[Number of files] \\Sexpr{length(rf)}     \\\
            \\end{description}\n",sep="")

      SweaveFile <- paste( SweaveFile, "\\setkeys{Gin}{width=0.8\\textwidth}\n
            \\begin{center}\n
            \n<<rfPlotUncalME, echo=F, fig=T, height=10, width=15>>=

            ## uncalibrated mass error per raw files
            uncalME <- tapply( evidence[, grep('^Uncalibrated.Mass.Error..ppm.$', colnames(evidence), ignore.case=T)], evidence[, grep('^Raw.file$', colnames(evidence),ignore.case=T)], function(x) x)

            #################################################
            # order according to raw file name
            #################################################
            uncalME <- uncalME[rf.order]
            ylim.tmp <- max(evidence[, grep('^Uncalibrated.Mass.Error..ppm.$', colnames(evidence), ignore.case=T)], na.rm=T)
            par(mar = c(20, 4, 4, 2)  )
            boxplot(uncalME, ylim=c(-ylim.tmp, ylim.tmp+ylim.tmp*.3), las=2, pch=16, cex=.3, main='Uncalibrated Mass Error (ppm)', cex.axis=1.5, ylab='ppm', col=barCol.tmp, cex.main=3)
            grid(nx=0, ny=NULL)

      \n@\n\\end{center}", sep="")
      SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")


      ###########################################################################################
      ##
      ##    - fraction of modified peptides
      ##    - overall number of detected site sper raw file
      ##    - number of site exclusively detected in certain raw files
      ##
      ###########################################################################################
      if(length(siteTabs) > 0){

          # loop over the modifications

          for(s in names(siteTabs)){


              ####################################################
              ##     percentage of modified evidences per raw file
              ## - use 'evidence' table
              ####################################################
              SweaveFile <- paste(SweaveFile, "
                  \\begin{frame}
                  \\frametitle{Raw Files: Enrichment Efficiency ",s,"sites}
                  \\begin{description}
                      \\item[Number of files] \\Sexpr{length(rf)}     \\\
                  \\end{description}\n

                  \\setkeys{Gin}{width=0.8\\textwidth}\n
                  \\begin{center}\n<<rfPlot_enrichment_",s,", echo=F, fig=T, height=10, width=15>>=


                         mod.enrich",s," <- vector('numeric', length(rf.order))
                         names( mod.enrich",s," ) <- rf.order

                         ## loop over raw files
                         for(rf.tmp in rf.order){

                             ev.rf.tmp <- evidence[ which( evidence[, grep('^Raw.File$', colnames(evidence), ignore.case=T) ] == rf.tmp), ]

                             mod.enrich",s,"[ rf.tmp ] <- length( grep('",s,"', ev.rf.tmp[, grep( '^Modifications$' , colnames(ev.rf.tmp), ignore.case=T) ])) / nrow(ev.rf.tmp)
                             rm(ev.rf.tmp)
                         }


                         # make the figure
                         par(mar = c(20,4,4,2)  )
                         fancyBarplot(as.vector(mod.enrich",s,"), main='Fraction of modified evidences (",s,")', ylim=c(0,1.3), ylab='Fraction', cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names(mod.enrich",s,"), cex.main=3, col=barCol, axes=F)
                         axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2))

                         abline(h=1, lwd=2, col='black')

                         # average enrichment efficiency
                         if(!is.null(experiments)){
                            av.enrich.eff <- unlist( lapply( experiments, function(x) mean( mod.enrich",s,"[ ed[ grep(paste('^',x, '$', sep=''), ed[, 3]), 1 ] ], na.rm=T) ))
                         } else {
                            av.enrich.eff <- mean( mod.enrich",s,", na.rm=T)
                         }

                        if(!is.null(experiments)){
                           #legend('top', legend=experiments, fill=cols.tmp, ncol=ifelse(length(experiments) > 8, ",leg.row,", 2), cex=",leg.cex,", bty='n'  )
                           legend('top', legend= paste(experiments, ' (', round(av.enrich.eff,3),')', sep=''), fill=cols.tmp, ncol=ceiling(length(experiments)/",leg.row,"), cex=",leg.cex,", bty='n'  )
                        } else {
                           legend('top', legend=round(av.enrich.eff,3), bty='n', cex=",leg.cex,")
                        }
                  par(mar=c(5,4,4,2), cex=1.5, bty='n')



                  \n@\n\\end{center}",sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")



              ####################################################
              #          detected  sites per raw file
              # - use 'evidence' table
              # - loop over the raw files and count the number of
              #   sites
              #
              ####################################################

              SweaveFile <- paste(SweaveFile, "
                  \\begin{frame}
                  \\frametitle{Raw Files: ",s,"sites}
                  \\begin{description}
                      \\item[Number of files] \\Sexpr{length(rf)}     \\\
                  \\end{description}\n

                  \\setkeys{Gin}{width=0.8\\textwidth}\n
                  \\begin{center}\n<<rfPlot_",s,", echo=F, fig=T, height=10, width=15>>=

                  # number of sites (evidences) per raw file
                  #sitesPerRF",s," <- evidence[unique(unlist(strsplit(peptides[ unique(unlist(strsplit(as.character(",s,"Sites[, grep('^Peptide.IDs$', colnames(",s,"Sites), ignore.case=T) ]), ';'))), 'Evidence.IDs'],';'))),  ]

                  # raw file names
                  rf.tmp <- unique(evidence[, grep('^Raw.File$', colnames(evidence), ignore.case=T)])
                  sitesPerRF",s," <- vector('numeric', length(rf.tmp))
                  names(sitesPerRF",s,") <- rf.tmp
                  rm(rf.tmp)

                  # loop over raw files
                  for(r in names(sitesPerRF",s,")){

                        # site ids of evidences of current raw file
                        ev.tmp <- unlist(strsplit(as.character(evidence[ grep(paste('^', r,'$', sep=''), evidence[,grep('^Raw.File$', colnames(evidence), ignore.case=T)]),  grep('^", mod.header.names[s],".Site.IDs$', colnames(evidence), ignore.case=T) ]), ';'))
                        # remove entries with no site id
                        ev.tmp <- ev.tmp[ which(nchar(ev.tmp) > 0) ]


                        sitesPerRF",s,"[r] <- length(unique(ev.tmp))
                        rm(ev.tmp)
                  }


                  sitesPerRF",s,"tmp <- vector('numeric', length(rf.order))
                  names(sitesPerRF",s,"tmp) <- rf.order
                  sitesPerRF",s,"tmp[ names(sitesPerRF",s,")  ] <- sitesPerRF",s,"
                  sitesPerRF",s,"tmp <- sitesPerRF",s,"tmp[rf.order]

                  # plot
                  par(mar = c(20,4,4,2)  )
                  fancyBarplot(as.vector(sitesPerRF",s,"tmp), main='Detected ",s," sites per raw file (evidence.txt)', ylim=c(0,max(sitesPerRF",s,"tmp, na.rm=T)+max(sitesPerRF",s,"tmp, na.rm=T)*.3), cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names(sitesPerRF",s,"tmp), cex.main=3, col=barCol, ylab='Counts')

                  if(!is.null(experiments))
                      legend('top', legend=experiments, fill=cols.tmp, ncol=ceiling(length(experiments)/",leg.row,"), cex=",leg.cex,", bty='n'  )
                      #legend('top', legend=experiments, fill=cols.tmp, ncol=ifelse(length(experiments) > 8, 3, 2), cex=1.5, bty='n'  )

                  par(mar=c(5,4,4,2))
                  \n@\n\\end{center}",sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")


              ####################################################
              #     exclusively detected  sites per raw file
              # - use 'evidence' table
              ####################################################
              SweaveFile <- paste(SweaveFile, "
                  \\begin{frame}
                  \\frametitle{Raw Files: Unique ",s,"sites}
                  \\begin{description}
                      \\item[Number of files] \\Sexpr{length(rf)}     \\\
                  \\end{description}\n

                  \\setkeys{Gin}{width=0.8\\textwidth}\n
                  \\begin{center}\n<<rfPlot_unique_",s,", echo=F, fig=T, height=10, width=15>>=



                         # store the raw files per detected sites
                         rf.per.site",s," <- vector('list', dim(",s,"Sites)[1])
                         names(rf.per.site",s,") <- rownames(", s,"Sites)

                         # loop over site ids -> get corresponding evidences -> get the raw file
                         for(ss in names(rf.per.site",s,"))
                                    rf.per.site",s,"[[ss]] <- unique(evidence[ grep(paste('(^|;)', ss, '($|;)', sep=''), evidence[ ,  grep('^", mod.header.names[s], ".Site.IDs$', colnames(evidence), ignore.case=T) ])  , grep('^Raw.File$', colnames(evidence), ignore.case=T)])

                         # get all sites that have been solely detected in a single raw file
                         exclusive.sites.rf.",s," <- table( unlist( rf.per.site",s,"[ which( unlist(lapply(rf.per.site",s,", length)) == 1)]  ))

                         # bring it in the same order as the figure before
                         exclusive.sites.rf.",s,".tmp <- vector('numeric', length(rf.order))
                         names(exclusive.sites.rf.",s,".tmp) <- rf.order
                         exclusive.sites.rf.",s,".tmp[names(exclusive.sites.rf.",s,")] <- exclusive.sites.rf.",s,"

                         # make the figure
                         par(mar = c(20,4,4,2)  )
                         fancyBarplot(as.vector(exclusive.sites.rf.",s,".tmp), main='Exclusively detected ",s," sites per raw file (evidence.txt)', ylim=c(0,max(exclusive.sites.rf.",s,".tmp, na.rm=T)+max(exclusive.sites.rf.",s,".tmp, na.rm=T)*.3), cex.axis=1.5, cex.lab=1.5, las=2, names.arg=names(exclusive.sites.rf.",s,".tmp), cex.main=3, col=barCol, ylab='Counts')

                        if(!is.null(experiments))
                           legend('top', legend=experiments, fill=cols.tmp, ncol=ceiling(length(experiments)/",leg.row,"), cex=",leg.cex,", bty='n'  )
                           #legend('top', legend=experiments, fill=cols.tmp, ncol=ifelse(length(experiments) > 8, 3, 2), cex=1.5, bty=\"n\"  )

                  par(mar=c(5,4,4,2))



                  \n@\n\\end{center}",sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")





          } # end for

      } # end if length(siteTbas)

    } # end rf
    #########################################
    #  overlap of peptide sequences between
    #  the experiments
    #########################################
    if(ol){
       if(!is.null(experiments) & length(experiments) > 1){

           ############################
           # calculate overlap
           ############################
           SweaveFile <- paste(SweaveFile, "
                \n<<overlapCalc, echo=F>>=
                overlap <- compareExperiments( evidence, file=NULL, experiments=NULL )
                \n@\n", sep="")

           ############################
           # ALL peptides
           ############################
           SweaveFile <- paste(SweaveFile, "\\section{Overlap}
                \\begin{frame}
                \\frametitle{Overlap between Experiments}
                \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                \n<<overlapPlot, echo=F, fig=T, height=8, width=10>>=

                overlap.peptides <- overlap[['overlap peptides']]
                # normalization
                overlap.peptides.norm <- overlap.peptides
                for(i in 1:dim(overlap.peptides)[1]){
                      overlap.peptides.norm[ i:dim(overlap.peptides)[1], i:dim(overlap.peptides)[1]  ] <- overlap.peptides.norm[ i:dim(overlap.peptides)[1], i:dim(overlap.peptides)[1] ] <- overlap.peptides[ i:dim(overlap.peptides)[1],  i:dim(overlap.peptides)[1]]/overlap.peptides[i,i]
                }
                # plot
                 heatmap.2( overlap.peptides.norm, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.peptides, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='Number of peptide sequences', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                \n@\n\\end{center}", sep="")
           SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

           ##############################
           #  SILAC
           ##############################
           if(quant & (MQversion.num <  12000)){

                ############################
                # LIGHT peptides
                ############################
                SweaveFile <- paste(SweaveFile, "\n
                \\begin{frame}
                \\frametitle{Overlap between Experiments (LIGHT peptides)}
                \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                  \n<<overlapPlot_pepL, echo=F, fig=T, height=8, width=10>>=

                  overlap.peptides.L <- overlap[['overlap peptides L']]
                  # normalization
                  overlap.peptides.norm.L <- overlap.peptides.L
                  for(i in 1:dim(overlap.peptides.L)[1]){
                      overlap.peptides.norm.L[ i:dim(overlap.peptides.L)[1], i:dim(overlap.peptides.L)[1]  ] <- overlap.peptides.norm.L[ i:dim(overlap.peptides.L)[1], i:dim(overlap.peptides.L)[1] ] <- overlap.peptides.L[ i:dim(overlap.peptides.L)[1],  i:dim(overlap.peptides.L)[1]]/overlap.peptides.L[i,i]
                  }
                  # plot
                  heatmap.2( overlap.peptides.norm.L, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.peptides.L, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='LIGHT peptide sequences', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                  \n@\n\\end{center}", sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

                ############################
                # MEDIUM peptides
                ############################
                if(silac.type == "triplets"){
                     SweaveFile <- paste(SweaveFile, "\n
                       \\begin{frame}
                        \\frametitle{Overlap between Experiments (MEDIUM peptides)}
                        \\setkeys{Gin}{width=0.75\\textwidth}\n
                        \\begin{center}
                        \n<<overlapPlot_pepM, echo=F, fig=T, height=8, width=10>>=

                          overlap.peptides.M <- overlap[['overlap peptides M']]
                          # normalization
                          overlap.peptides.norm.M <- overlap.peptides.M
                          for(i in 1:dim(overlap.peptides.M)[1]){
                                 overlap.peptides.norm.M[ i:dim(overlap.peptides.M)[1], i:dim(overlap.peptides.M)[1]  ] <- overlap.peptides.norm.M[ i:dim(overlap.peptides.M)[1], i:dim(overlap.peptides.M)[1] ] <- overlap.peptides.M[ i:dim(overlap.peptides.M)[1],  i:dim(overlap.peptides.M)[1]]/overlap.peptides.M[i,i]
                          }
                          # plot
                          heatmap.2( overlap.peptides.norm.M, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.peptides.M, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='MEDIUM peptide sequences', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                      \n@\n\\end{center}", sep="")
                      SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

                }
                ############################
                # HEAVY peptides
                ############################
                SweaveFile <- paste(SweaveFile, "\n
                \\begin{frame}
                  \\frametitle{Overlap between Experiments (HEAVY peptides)}
                  \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                  \n<<overlapPlot_pepH, echo=F, fig=T, height=8, width=10>>=

                  overlap.peptides.H <- overlap[['overlap peptides H']]
                  # normalization
                  overlap.peptides.norm.H <- overlap.peptides.H
                  for(i in 1:dim(overlap.peptides.H)[1]){
                      overlap.peptides.norm.H[ i:dim(overlap.peptides.H)[1], i:dim(overlap.peptides.H)[1]  ] <- overlap.peptides.norm.H[ i:dim(overlap.peptides.H)[1], i:dim(overlap.peptides.H)[1] ] <- overlap.peptides.H[ i:dim(overlap.peptides.H)[1],  i:dim(overlap.peptides.H)[1]]/overlap.peptides.H[i,i]
                  }
                  # plot
                  heatmap.2( overlap.peptides.norm.H, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.peptides.H, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='HEAVY peptide sequences', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                  \n@\n\\end{center}", sep="")
                SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

           }

           ######################################
           # protein groups detected with UNIQUE
           # in the particular experiment
           ######################################
           SweaveFile <- paste(SweaveFile, "
                \\begin{frame}
                \\frametitle{Overlap between Experiments}

                \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                \n<<overlapPlot_PG, echo=F, fig=T, height=8, width=10>>=
                overlap.proteins <- overlap[['overlap proteins']]
                # normalization
                overlap.proteins.norm <- overlap.proteins
                for(i in 1:dim(overlap.proteins)[1]){
                      overlap.proteins.norm[ i:dim(overlap.proteins)[1], i:dim(overlap.proteins)[1]  ] <- overlap.proteins.norm[ i:dim(overlap.proteins)[1], i:dim(overlap.proteins)[1] ] <- overlap.proteins[ i:dim(overlap.proteins)[1],  i:dim(overlap.proteins)[1]]/overlap.proteins[i,i]
                }
                # plot
                 heatmap.2( overlap.proteins.norm, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.proteins, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='Number of protein groups detected \\n by unique evidences', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                \n@\n\\end{center}", sep="")
           SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

           #################
           # SILAC
           #################
           if(quant & (MQversion.num <  12000)){
               ##################################
               # LIGHT proteins
               ##################################
               SweaveFile <- paste(SweaveFile, "
                \\begin{frame}
                \\frametitle{Overlap between Experiments (LIGHT protein groups)}

                \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                  \n<<overlapPlot_PGL, echo=F, fig=T, height=8, width=10>>=
                  overlap.proteins.L <- overlap[['overlap proteins L']]
                  # normalization
                  overlap.proteins.norm.L <- overlap.proteins.L
                  for(i in 1:dim(overlap.proteins.L)[1]){
                      overlap.proteins.norm.L[ i:dim(overlap.proteins.L)[1], i:dim(overlap.proteins.L)[1]  ] <- overlap.proteins.norm.L[ i:dim(overlap.proteins.L)[1], i:dim(overlap.proteins.L)[1] ] <- overlap.proteins.L[ i:dim(overlap.proteins.L)[1],  i:dim(overlap.proteins.L)[1]]/overlap.proteins.L[i,i]
                  }
                  # plot
                  heatmap.2( overlap.proteins.norm.L, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.proteins.L, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='LIGHT protein groups', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                \n@\n\\end{center}", sep="")
              SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

               #####################################
               # MEDIUM proteins
               #####################################
               if(silac.type == "triplets"){

                 SweaveFile <- paste(SweaveFile, "
                   \\begin{frame}
                   \\frametitle{Overlap between Experiments (MEDIUM protein groups)}

                   \\setkeys{Gin}{width=0.75\\textwidth}\n
                   \\begin{center}
                    \n<<overlapPlot_PGM, echo=F, fig=T, height=8, width=10>>=
                    overlap.proteins.M <- overlap[['overlap proteins M']]
                    # normalization
                    overlap.proteins.norm.M <- overlap.proteins.M
                    for(i in 1:dim(overlap.proteins.M)[1]){
                      overlap.proteins.norm.M[ i:dim(overlap.proteins.M)[1], i:dim(overlap.proteins.M)[1]  ] <- overlap.proteins.norm.M[ i:dim(overlap.proteins.M)[1], i:dim(overlap.proteins.M)[1] ] <- overlap.proteins.M[ i:dim(overlap.proteins.M)[1],  i:dim(overlap.proteins.M)[1]]/overlap.proteins.M[i,i]
                    }
                    # plot
                    heatmap.2( overlap.proteins.norm.M, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.proteins.M, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='MEDIUM protein groups', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                  \n@\n\\end{center}", sep="")
                SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

               }
               ##################################
               # HEAVY proteins
               ##################################
               SweaveFile <- paste(SweaveFile, "
                \\begin{frame}
                \\frametitle{Overlap between Experiments (HEAVY protein groups)}

                \\setkeys{Gin}{width=0.75\\textwidth}\n
                  \\begin{center}
                  \n<<overlapPlot_PGH, echo=F, fig=T, height=8, width=10>>=
                  overlap.proteins.H <- overlap[['overlap proteins H']]
                  # normalization
                  overlap.proteins.norm.H <- overlap.proteins.H
                  for(i in 1:dim(overlap.proteins.H)[1]){
                      overlap.proteins.norm.H[ i:dim(overlap.proteins.H)[1], i:dim(overlap.proteins.H)[1]  ] <- overlap.proteins.norm.H[ i:dim(overlap.proteins.H)[1], i:dim(overlap.proteins.H)[1] ] <- overlap.proteins.H[ i:dim(overlap.proteins.H)[1],  i:dim(overlap.proteins.H)[1]]/overlap.proteins.H[i,i]
                  }
                  # plot
                  heatmap.2( overlap.proteins.norm.H, Rowv=F, Colv='Rowv', dendrogram='none', trace='none', cellnote=overlap.proteins.H, margins=c(10,10), col=topo.colors(20), density.info='none', notecol='black', notecex=1.5, main='HEAVY protein groups', cexRow=1.5, cexCol=1.5, breaks=seq(0,1,0.05), symkey=F)
                \n@\n\\end{center}", sep="")
              SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

           }

           #############################################################
           #
           #          detected peptides per number of experiments
           #
           ##############################################################
           if(quant & (MQversion.num <  12000)){
               #################################################
               # SILAC: make a side by side barchart
               #################################################

              if(silac.type == "doublets"){
                  SweaveFile <- paste(SweaveFile, "
                     \n<<pepPerExp_doubleSILAC, echo=F>>=
                       pepPerExpMat <- rbind(overlap[['# peptides']], overlap[['# peptides L']], overlap[['# peptides H']] )
                       rownames(pepPerExpMat) <- c('all', 'light', 'heavy')
                     \n@\n
                  ",sep="")
              } else if(silac.type == "triplets"){
                  SweaveFile <- paste(SweaveFile, "
                     \n<<pepPerExp_tripleSILAC, echo=F>>=
                       pepPerExpMat <- rbind(overlap[['# peptides']], overlap[['# peptides L']], overlap[['# peptides M']], overlap[['# peptides H']] )
                       rownames(pepPerExpMat) <- c('all', 'light', 'medium', 'heavy')
                     \n@\n
                  ",sep="")
              }
              ########################################
              # make the plot
              ########################################
              SweaveFile <- paste(SweaveFile, "
                \\begin{frame}
                \\frametitle{Shared Peptides between Experiments}

                \\setkeys{Gin}{width=1\\textwidth}\n
                  \\begin{center}
                \n<<pepPerExp_plot, echo=F, fig=T, height=7, width=14>>=
                # colors
                if(!is.null(dim(pepPerExpMat)) ){
                       col.tmp <- grey(seq(0.1,0.9,length.out=dim(pepPerExpMat)[1]))
                }
                # plot
                fancyBarplot(pepPerExpMat, main='Peptide sequences detected in N experiments', xlab='Number of experiments (N)', ylab='Number', col=col.tmp)
                legend('top', legend=rownames(pepPerExpMat), fill=col.tmp, bty='n')
                \n@\n\\end{center}", sep="")
              SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")

           } else{
               #############################################
               #  no SILAC
               #############################################
               SweaveFile <- paste(SweaveFile, "
                    \\begin{frame}
                    \\frametitle{Shared Peptides between Experiments}

                    \\setkeys{Gin}{width=0.75\\textwidth}\n
                      \\begin{center}
                     \n<<pepPerExp_plot, echo=F, fig=T, height=8, width=10>>=
                     pepPerExp <- overlap[['# peptides']]

                     # plot
                     fancyBarplot(pepPerExp, main='Peptide sequences detected in N experiments', xlab='Number of experiments (N)', ylab='Number')
                   \n@\n\\end{center}", sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")
           }



           ##################################################################
           #
           #          detected protein groups per experiment
           #
           ##################################################################
           if(quant & (MQversion.num <  12000) ){
               #################################################
               # SILAC: make a side by side barchart
               #################################################
              if(silac.type == "doublets"){

                  SweaveFile <- paste(SweaveFile, "
                     \n<<protPerExp_doubleSILAC, echo=F>>=
                       protPerExpMat <- rbind(overlap[['# proteins']], overlap[['# proteins L']], overlap[['# proteins H']] )
                       rownames(protPerExpMat) <- c('all', 'light', 'heavy')
                     \n@\n
                  ",sep="")

              } else if(silac.type == "triplets"){

                  SweaveFile <- paste(SweaveFile, "
                     \n<<protpPerExp_tripleSILAC, echo=F>>=
                       protPerExpMat <- rbind(overlap[['# proteins']], overlap[['# proteins L']], overlap[['# proteins M']], overlap[['# proteins H']] )
                       rownames(protPerExpMat) <- c('all', 'light', 'medium', 'heavy')
                     \n@\n
                  ",sep="")
              }
              ########################################
              # make the plot
              ########################################
              SweaveFile <- paste(SweaveFile, "
                \\begin{frame}
                \\frametitle{Shared Protein Groups between Experiments}

                \\setkeys{Gin}{width=1\\textwidth}\n
                  \\begin{center}
                \n<<protPerExp_plot, echo=F, fig=T, height=7, width=12>>=
                # plot
                fancyBarplot(protPerExpMat, main='Protein groups detected in N experiments', xlab='Number of experiments (N)', ylab='Number', col=col.tmp)
                legend('top', legend=rownames(protPerExpMat), fill=col.tmp, bty='n' )
                \n@\n\\end{center}", sep="")

              SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")


           } else {
               ##################################
               # no SILAC
               ##################################
               SweaveFile <- paste(SweaveFile, "
                    \\begin{frame}
                    \\frametitle{Shared Protein Groups between Experiments}

                    \\setkeys{Gin}{width=0.75\\textwidth}\n
                      \\begin{center}
                      \n<<protPerExp_plot, echo=F, fig=T, height=8, width=10>>=
                      protPerExp <- overlap[['# proteins']]

                      # plot
                      fancyBarplot(protPerExp, main='Protein groups detected in N experiments', xlab='N', ylab='Number')
                   \n@\n\\end{center}", sep="")
               SweaveFile <- paste( SweaveFile ,"\\end{frame}\n", sep="")
           }
        }
    }


    #########################################
    # separation
    #########################################
    if(separation){
      if(!is.null(experiments))
      {
          for(ex in experiments)
          {
            SweaveFile = paste(SweaveFile, "\n\\begin{frame}
               \\frametitle{Separation: ",experiments.tex[ex],"}\n
               \\Sexpr{sum(  ed[, 'Experiment'] %in% '",experiments.dot[ex],"' )} Raw Files
               \\begin{center}
               \\setkeys{Gin}{width=0.6\\textwidth}
               \n<<separation-",experiments.tex[ex],", echo=F, fig=T, height=8, width=10>>=\n
               checkSeparation( evidence[which(evidence[, grep('^Experiment$', colnames(evidence), ignore.case=T)] == '",experiments.dot[ex],"') , ] )\n@
                \\end{center}
                \\end{frame}\n", sep="")
          }
      }
      else{

         SweaveFile = paste(SweaveFile, "\n\\begin{frame}
               \\frametitle{Separation}\n
               \\begin{center}
               \\setkeys{Gin}{width=0.6\\textwidth}
               \n<<separation, echo=F, fig=T, height=8, width=10>>=\n
               checkSeparation( evidence )\n@
                \\end{center}
                \\end{frame}\n", sep="")
      }
    }
    #############################################
    ## Histogram of m/z values
    ## 20150112 - switched to from MS/MS to evidences
    #############################################
    SweaveFile = paste(SweaveFile, "\n\\begin{frame}
               \\frametitle{Histogram of m/z values}\n
               Based on identified evidences (without con/rev).
               \\begin{center}
               \\setkeys{Gin}{width=0.6\\textwidth}
               \n<<mzHist, echo=F, fig=T, height=8, width=10>>=\n
               msms.mz <- evidence$m.z
               hist(msms.mz, xlim=c(300,2000), breaks=17, axes=F, col='grey', main='', xlab='m/z')
               axis(2)
               axis(1, at=seq(300, 2000, 100))
               q.tmp <- quantile(msms.mz, na.rm=T)
               abline(v=q.tmp, lwd=3, lty='dashed', col='grey')
               legend('topright', legend=paste(names(q.tmp), ': ',round(q.tmp,2)), bty='n', title='Quartiles', cex=2 )
               rm(q.tmp)
               \n@\n
                \\end{center}
                \\end{frame}\n", sep="")


    #############################################
    # Charge state & peptide length distribution
    #############################################
    SweaveFile = paste(SweaveFile, "\\section{Evidence Table}\n\\begin{frame}
               \\frametitle{Evidence table}\n
               Based on evidences without contaminants/reverse hits (\\Sexpr{dim(evidence)[1]})
               \\begin{center}
               \\setkeys{Gin}{width=0.9\\textwidth}
               \n<<charge_state, echo=F, fig=T, height=7, width=12>>=\n
               par(mfrow=c(1,2))
               ###################
               # peptide length
               ###################
               fancyDensPlot(evidence$Length, main='Peptide length', xlab='Length in amino acids')
               legend('topright', legend=c(paste('Mean:', round(mean(evidence$Length),2)), paste('Median:', round(median(evidence$Length),2)) ), bty='n', cex=2)
               ###################
               # charge state
               ###################
               cs <- table(evidence[ ,grep('^Charge$', colnames(evidence), ignore.case=T)])
               fancyBarplot(as.vector(cs), names=names(cs), xlab='Charge', ylab='Counts', main='Charge state distribution')
               legend('topright', legend=paste('charge ', names(cs), ': ', round(100*(cs/sum(cs)),1), ' %', sep=''), bty='n', cex=1.5 )

               rm(cs)

               \n@\n
                \\end{center}
                \\end{frame}\n", sep="")

    #############################################
    # Evidence type & MS/MS count
    #############################################
    SweaveFile = paste(SweaveFile, "\n\\begin{frame}
               \\frametitle{Evidence table}\n
               Based on evidences without contaminants/reverse hits (\\Sexpr{dim(evidence)[1]})
               \\begin{center}
               \\setkeys{Gin}{width=0.9\\textwidth}
               \n<<evidence_type, echo=F, fig=T, height=7, width=12>>=\n
               par(mfrow=c(1,2))
               ###################
               # evidence type
               ###################
               et <- table(evidence$Type)
               fancyBarplot(as.vector(et), names=names(et), main='Evidence type')

               ###################
               # MSMS count
               ###################
               emsms <- table(evidence$MS.MS.Count)

               max.msms <- 10
               emsms.tmp <- emsms[  which(as.numeric(names(emsms)) < max.msms)  ]
               emsms.tmp <- c(emsms.tmp, sum( emsms[  which(as.numeric(names(emsms)) >= max.msms)  ]  ) )
               names(emsms.tmp)[length(emsms.tmp)] <- paste('>', max.msms)

               emsms <- emsms.tmp

               fancyBarplot(as.vector(emsms), names=names(emsms), xlab='MS/MS count', ylab='Counts', main='MS/MS count')
               legend('topright', legend=paste( names(emsms), ': ', round(100*(emsms/sum(emsms)),1), ' %', sep=''), bty='n', cex=1.5 )

               rm(et, emsms)
               \n@\n
                \\end{center}
                \\end{frame}\n", sep="")



    #########################################
    #   Mass Accuracy
    #########################################
    if(ma){

      #######################
      #  calibrated
      #######################
      SweaveFile = paste(SweaveFile, "\\section{Mass Deviation}\n\\begin{frame}
               \\frametitle{Precursor Ion Mass Deviation}\n
               \\alert{Calibrated} mass error of evidences without contaminants/reverse hits (\\Sexpr{dim(evidence)[1]}).
               \\begin{center}
               \\setkeys{Gin}{width=0.95\\textwidth}
               \n<<masserror, echo=F, fig=T, height=8, width=16>>=\n
               par(mfrow=c(1,2))\n

               me <- evidence[, 'Mass.Error..ppm.']
               # check whether there are 'n. def.'s ... (version 1.1.1.26)
               nd.idx <- grep('n. def.', me)
               if(length(nd.idx)>0)
                      me[nd.idx] <- NA
               me <- as.numeric(me)

               fancyDensPlot( me , main='Mass error (calibrated)', xlab='ppm', cex.axis=1.5, cex.lab=1.5)
               legend('topright', legend=c(paste('mean: ', round(mean(me, na.rm=T),3), sep=''),  paste('median: ', round(median(me, na.rm=T),3), sep=''), paste('sd: ', round(sd(me, na.rm=T),3),sep='') ), cex=1.5  ) \n
               if( sum(is.na(me)) > 0 ) legend('topleft', legend=paste('NA/NAN values:', sum(is.nan(me))), cex=1.2 )
               fancyDensPlot(abs( me ), main='Absolute mass error (calibrated)', xlab='ppm', cex.axis=1.5, cex.lab=1.5)
               legend('topright', legend=c(paste('mean: ', round(mean(abs(me), na.rm=T),3), sep=''),  paste('median: ', round(median(abs(me), na.rm=T),3), sep=''), paste('99th %', round(quantile(abs(me), c(0.99), na.rm=T),3),sep='') ), cex=1.5  ) \n

               par(mfrow=c(1,1))\n
               rm(me)
               \n@
                \\end{center}
                \\end{frame}\n", sep="")


      ########################
      #  uncalibrated
      ########################
      SweaveFile = paste(SweaveFile, "\n\\begin{frame}
               \\frametitle{Precursor Ion Mass Deviation}\n
               \\alert{Uncalibrated} mass error of evidences without contaminants/reverse hits (\\Sexpr{dim(evidence)[1]}).
               \\begin{center}
               \\setkeys{Gin}{width=0.95\\textwidth}
               \n<<masserror_uncal, echo=F, fig=T, height=8, width=16>>=\n
               par(mfrow=c(1,2))\n

               me <- evidence[, 'Uncalibrated.Mass.Error..ppm.']
               # check whether there are 'n. def.'s ... (version 1.1.1.26)
               nd.idx <- grep('n. def.', me)
               if(length(nd.idx)>0)
                      me[nd.idx] <- NA
               me <- as.numeric(me)

               fancyDensPlot( me , main='Mass error (uncalibrated)', xlab='ppm', cex.axis=1.5, cex.lab=1.5)
               legend('topright', legend=c(paste('mean: ', round(mean(me, na.rm=T),3), sep=''),  paste('median: ', round(median(me, na.rm=T),3), sep=''), paste('sd: ', round(sd(me, na.rm=T),3),sep='') ), cex=1.5  ) \n
               if( sum(is.na(me)) > 0 ) legend('topleft', legend=paste('NA/NAN values:', sum(is.na(me))), cex=1.2 )
               fancyDensPlot(abs( me ), main='Absolute mass error (uncalibrated)', xlab='ppm', cex.axis=1.5, cex.lab=1.5)
               legend('topright', legend=c(paste('mean: ', round(mean(abs(me), na.rm=T),3), sep=''),  paste('median: ', round(median(abs(me), na.rm=T),3), sep=''), paste('99th %', round(quantile(abs(me), c(0.99), na.rm=T),3),sep='') ), cex=1.5  ) \n

               par(mfrow=c(1,1))\n
               rm(me)

               \n@
                \\end{center}
                \\end{frame}\n", sep="")


      if(msms){

          ################################
          ##  MS/MS fragment mass error
          ################################
          for(ft in names(frag.types)){

              SweaveFile = paste(SweaveFile, "\n\\begin{frame}
                   \\frametitle{Fragment Ion Mass Deviation \\alert{",ft,"}}\n
                   \\begin{figure}
                   \\begin{center}
                   \\setkeys{Gin}{width=0.7\\textwidth}
                   \n<<masserror_msms_ppm_", ft,", echo=F, fig=F>>=

                      msms.me.Da.tmp <- as.numeric(unlist(strsplit(msms.me.Da[which(msms.frag == '", ft,"')] , ';')))

                      if( 'msms.me.ppm' %in% ls() ){
                          msms.me.ppm.tmp <- as.numeric(unlist(strsplit(msms.me.ppm[which(msms.frag == '", ft,"')] , ';')))
                      } else {
                          # number of fragment ions per spectrum
                          msms.frag.N.tmp <- unlist(lapply( strsplit(msms.me.Da[ which(msms.frag == '", ft,"') ], ';'), length) )

                          #msms.masses.tmp <- as.numeric( rep(msms.masses[ which(msms.frag == '", ft,"') ], msms.frag.N.tmp) )
                          msms.masses.tmp <- as.numeric( unlist(strsplit(msms.masses[which(msms.frag == '", ft,"')] , ';')) )

                          # ppm
                          msms.me.ppm.tmp <- msms.me.Da.tmp / (msms.masses.tmp * 1e-6)
                      }



                      dens.col.tmp <- densCols(msms.me.ppm.tmp, msms.me.Da.tmp, colramp=colorRampPalette(c('lightblue','green','lightgreen', 'yellow', 'orange', 'red', 'darkred')))

                      png('pic/plot-masserror_msms_ppm_", ft,".png', width=480, height=480, res=100 )

                      plot(msms.me.ppm.tmp, msms.me.Da.tmp, main='ppm vs. Da', xlab='Mass error (ppm)', ylab='Mass error (Da)', pch=20, col=dens.col.tmp)
                      abline(v=0, h=0, lwd=2, lty='dashed', col='grey')
                      legend('topleft', legend=c(paste('abs. average=', round(mean( abs(msms.me.Da.tmp)),2),' Da', sep='')), bty='n', title='Da', cex=0.8)
                      legend('bottomright', legend=c(paste('abs. average=', round(mean( abs(msms.me.ppm.tmp)),2),' ppm', sep='')), bty='n', title='ppm', cex=0.8)

                      invisible( dev.off() )

                      #rm(msms.me.ppm.tmp, msms.me.Da.tmp, dens.col.tmp)

                   \n@\n
                   \\includegraphics{pic/plot-masserror_msms_ppm_", ft,"}
                   \\end{center}
                   \\end{figure}

                     ",frag.types[ft]," spectra/\\Sexpr{ length(msms.me.Da.tmp)} fragments\n
                     \\Sexpr{ round( length(msms.me.Da.tmp) /",frag.types[ft],", 1)} fragments per spectrum on average

                   \\end{frame}\n", sep="")
          }

      } ## end if msms

  }
    #########################################
    #  PEP and Mascot Score distributions
    #########################################
    if(pep){
      if(msms){
      SweaveFile = paste(SweaveFile, "\\section{Score Distributions}\n\\begin{frame}
               \\frametitle{PEP Distribution}\n
               Distribution of posterior error probabilities based on all identifications without contaminants/reverse hits.
               \\begin{center}
               \\setkeys{Gin}{width=0.9\\textwidth}
               \n<<PEPdist, echo=F, fig=F, height=8, width=16>>=\n
               png('pic/plot-PEPdist.png', width=960, height=480, res=100 )
               par(mar=c(8,5,2,2))
               boxplot(as.numeric(proteinGroups[, grep('^PEP$', colnames(proteinGroups), ignore.case=T)]), as.numeric(peptides[ , grep('^PEP$', colnames(peptides), ignore.case=T) ]), as.numeric(evidence[, grep('^PEP$', colnames(evidence), ignore.case=T)]), as.numeric(msmsPEP), names=c('protein groups', 'peptides', 'evidences', 'msms'), col='lightblue', ylab='PEP', cex.axis=1.2, las=2)
               invisible( dev.off() )
               \n@
                \\includegraphics{pic/plot-PEPdist}
                \\end{center}
                \\end{frame}\n", sep="")


       } else {
      SweaveFile = paste(SweaveFile, "\\section{Score Distributions}\n\\begin{frame}
               \\frametitle{PEP Distribution}\n
               Distribution of posterior error probabilities.
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<PEPdist, echo=F, fig=T, height=8, width=16>>=\n
               boxplot( as.numeric(proteinGroups[, grep('^PEP$', colnames(proteinGroups), ignore.case=T)]), as.numeric(peptides[ , grep('^PEP$', colnames(peptides), ignore.case=T) ]), as.numeric(evidence[, grep('^PEP$', colnames(evidence), ignore.case=T)]), names=c('protein groups', 'peptides', 'evidences'), col='lightblue', ylab='PEP', cex.axis=1.5)
               \n@
                \\end{center}
                \\end{frame}\n", sep="")

       }

      #############################
      # Mascot score (evidence.txt)
      #############################
       if( length( grep("^1.0", MQversion) ) > 0){
           SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{Mascot Distribution}\n
               Mascot score distribution of evidences without contaminants/reverse hits.
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<Mascotdist, echo=F, fig=T, height=6, width=12>>=\n
               mascot.score <- as.numeric(evidence[, 'Mascot.Score'])
               mascot.score.mean <- mean(mascot.score, na.rm=T)
               mascot.score.median <- median(mascot.score, na.rm=T)
               mascot.score.99 <- quantile( mascot.score, 0.99, na.rm=T)
               par(mfrow=c(1,2))
               plot( density( mascot.score, na.rm=T  ), xlab='Mascot score', lwd=3, main='density plot' )
               legend('topright', legend=c(paste('mean:', round(mascot.score.mean,2) ), paste('median:', round(mascot.score.median,2) ), paste('99th percentile:', round(mascot.score.99,2) ) ), cex=1.5   )
               boxplot( mascot.score, main='boxplot' )
               \n@
                \\end{center}
                \\end{frame}\n", sep="")
       } else{
           SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{Andromeda Distribution}\n
               Andromeda score distribution of evidences without contaminants / reverse hits.
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<Andromedadist, echo=F, fig=T, height=6, width=12>>=\n
               score <- as.numeric( evidence[, 'Score'])
               score.mean <- mean(score, na.rm=T)
               score.median <- median(score, na.rm=T)
               score.99 <- quantile( score, 0.99, na.rm=T)
               par(mfrow=c(1,2))
               plot( density( score, na.rm=T  ), xlab='Andromeda score', lwd=3, main='density plot' )
               legend('topright', legend=c(paste('mean:', round(score.mean,2) ), paste('median:', round(score.median,2) ), paste('99th percentile:', round(score.99,2) ) )   )
               boxplot( score, main='boxplot' )
               \n@
                \\end{center}
                \\end{frame}\n", sep="")

           ###############################################################
           ## separate for each experiment
           ## SweaveFile = paste(SweaveFile, "\\begin{frame}
           ##     \\frametitle{Andromeda Distribution II}\n
           ##     Andromeda score distribution of evidences without contaminants / reverse hits.
           ##     \\begin{center}
           ##     \\setkeys{Gin}{width=0.8\\textwidth}
           ##     \n<<Andromedadist_sep, echo=F, fig=T, height=6, width=12>>=\n

           ##     scoreL <- tapply( evidence[, 'Score'], evidence[, 'Experiment'], function(x) as.numeric(x))

           ##     ## score <- as.numeric( evidence[, 'Score'])
           ##     scoreL.mean <- lapply(scoreL, mean, na.rm=T)
           ##     score.median <- lapply(scoreL, median, na.rm=T)
           ##     ##score.99 <- quantile( score, 0.99, na.rm=T)
           ##     ##par(mfrow=c(1,2))

           ##     boxplot(scoreL)
           ##     ##plot( density( score, na.rm=T  ), xlab='Andromeda score', lwd=3, main='density plot' )
           ##     ##legend('topright', legend=c(paste('mean:', round(score.mean,2) ), paste('median:', round(score.median,2) ), paste('99th percentile:', round(score.99,2) ) )   )
           ##     ##boxplot( score, main='boxplot' )
           ##     \n@
           ##      \\end{center}
           ##      \\end{frame}\n", sep="")
       }


    }
    #####################################################
    #  label free quantitaion: only for two experiments
    #####################################################
    if(lfq){

        ######################################
        # determine a ratio
        ######################################
        SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{Label Free Quantitation}\n

               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}

               \n<<LFQ, echo=F, fig=T, height=7, width=7>>=\n
                   if(length(grep('^LFQ', colnames(proteinGroups))) == 2){
                           LFQ <- proteinGroups[, grep('^LFQ', colnames(proteinGroups))]
                           LFQ.ratio <- LFQ[, 1]/LFQ[,2]
                   }
                   if(length(grep('^Norm.Intensity', colnames(proteinGroups))) == 2){
                           LFQ <- proteinGroups[, grep('^Norm.Intensity', colnames(proteinGroups))]
                           LFQ.ratio <- LFQ[, 1]/LFQ[,2]
                   }
                   ####################################
                   # append to protein groups table
                   ####################################
                   LFQ.ratio[is.infinite(LFQ.ratio)] <- NA
                   LFQ.ratio[LFQ.ratio == 0] <- NA

                   proteinGroups.lfq <- cbind(proteinGroups, LFQ.ratio)
                   colnames(proteinGroups.lfq)[dim(proteinGroups.lfq)[2]] <- paste(colnames(LFQ)[1], '/', colnames(LFQ)[2], sep='')

                   ####################################
                   # export that table
                   ####################################
                   write.table(proteinGroups.lfq, file='proteinGroups_LFQ.txt', sep='\t', quote=F, col.names=NA, na='')

                   ####################################
                   # ratio vs. intensity plot
                   ####################################
                   tmp <- sigRatioPlot.ph(rat=LFQ.ratio, int=proteinGroups.lfq[, 'Intensity'],  p = rep(1, dim(proteinGroups.lfq)[1]),main=paste( colnames(proteinGroups.lfq)[dim(proteinGroups.lfq)[2]] ), legend.left=F, legend.right=F )

               \n@\n\\end{center}\\end{frame}\n", sep="")
        }


    #########################################
    #  SILAC
    #########################################
    if(quant){

        ############################################
        #           overview
        ############################################
        # mixing error
        # ...
        SweaveFile <- paste( SweaveFile, "\\section{SILAC}\n")

        if(!is.null(experiments)){

           if(silac.type=="doublets"){

             SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{SILAC Mixing Error}\n
                \\alert{Median} of unnormalized \\alert{evidence} ratios on \\alert{raw scale}.
                \n<<mixingerror, echo=F, results=tex>>=
                tab.exp.tmp <- table(evidence$Experiment)
                mix.error.tab <- matrix( 0, nrow=3, ncol=length(tab.exp.tmp), dimnames=list(c('Median(H/L)', 'N detected', 'N quantified'), names(tab.exp.tmp))  )
                mix.error.tab['N detected',] <- tab.exp.tmp

                # calculate mixing error
                for(ee in names(tab.exp.tmp)){
                       rat.tmp <- evidence[which(evidence[ ,grep('^Experiment$', colnames(evidence), ignore.case=T)] == ee), grep('^Ratio.H.L$', colnames(evidence), ignore.case=T)]
                       rat.tmp[nchar(rat.tmp) == 0] <- NA
                       #rat.tmp <- log(as.numeric(rat.tmp),2)
                       rat.tmp <- as.numeric(rat.tmp)

                       mix.error.tab[1, ee] <- round( median(rat.tmp, na.rm=T),2)
                       mix.error.tab['N quantified', ee] <- sum(!is.na(rat.tmp) )
                       rm(rat.tmp)
                }
                # latex table
                print(xtable(t(mix.error.tab), digits=c(0,2,0,0) ), size='tiny' )

                rm(tab.exp.tmp)
                \n@\n
               \\end{frame}\n", sep="")
           } # end if doublets

           if(silac.type=="triplets"){

             SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{SILAC Mixing Error}\n
               \\alert{Median} of unnormalized \\alert{evidence} ratios on \\alert{raw scale}.
                \n<<mixingerror, echo=F, results=tex>>=
                tab.exp.tmp <- table(evidence$Experiment)
                mix.error.tab <- matrix( 0, nrow=5, ncol=length(tab.exp.tmp), dimnames=list(c('M/L', 'H/L', 'H/M', 'N detected', 'N quantified'), names(tab.exp.tmp))  )
                mix.error.tab['N detected',] <- tab.exp.tmp

                # calculate mixing error
                for(ee in names(tab.exp.tmp)){
                       # /H/L
                       rat.tmp.hl <- evidence[which(evidence[,grep('^Experiment$', colnames(evidence), ignore.case=T)] == ee),  grep('^Ratio.H.L$', colnames(evidence), ignore.case=T)]
                       rat.tmp.hl[nchar(rat.tmp.hl) == 0] <- NA
                       #rat.tmp.hl <- log(as.numeric(rat.tmp.hl),2)
                       rat.tmp.hl <- as.numeric(rat.tmp.hl)

                       # M/L
                       rat.tmp.ml <- evidence[which(evidence[,grep('^Experiment$', colnames(evidence), ignore.case=T)] == ee), grep('^Ratio.M.L$', colnames(evidence), ignore.case=T)]
                       rat.tmp.ml[nchar(rat.tmp.ml) == 0] <- NA
                       #rat.tmp.ml <- log(as.numeric(rat.tmp.ml),2)
                       rat.tmp.ml <- as.numeric(rat.tmp.ml)

                       # H/M
                       rat.tmp.hm <- evidence[which(evidence[,grep('^Experiment$', colnames(evidence), ignore.case=T)] == ee),  grep('^Ratio.H.M$', colnames(evidence), ignore.case=T)]
                       rat.tmp.hm[nchar(rat.tmp.hm) == 0] <- NA
                       #rat.tmp.hm <- log(as.numeric(rat.tmp.hm),2)
                       rat.tmp.hm <- as.numeric(rat.tmp.hm)


                       mix.error.tab['M/L', ee] <- round( median(rat.tmp.ml, na.rm=T),2)
                       mix.error.tab['H/L', ee] <- round( median(rat.tmp.hl, na.rm=T),2)
                       mix.error.tab['H/M', ee] <- round( median(rat.tmp.hm, na.rm=T),2)
                       mix.error.tab['N quantified', ee] <- sum(!is.na(rat.tmp.hl) )
                       rm(rat.tmp.hl,rat.tmp.ml,rat.tmp.hm)
                }
                # latex table
                print(xtable(t(mix.error.tab), digits=c(0,2,2,2,0,0)), size='tiny' )

                rm(tab.exp.tmp )
                \n@\n
               \\end{frame}\n", sep="")
           } # end if doublets


        } else {   # end if experiments
            #################################
            # no experiments: doublets
            #################################
            if(silac.type=="doublets"){

             SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{SILAC Mixing Error}\n
               \\alert{Median} of unnormalized \\alert{evidence} ratios on \\alert{raw scale}.

                \n<<mixingerror, echo=F, results=tex>>=
                mix.error.tab <- matrix( 0, nrow=3, ncol=1, dimnames=list(c('H/L', 'N detected', 'N quantified'), 'Ratio.H.L')  )
                mix.error.tab['N detected',1] <- dim(evidence)[1]

                # calculate mixing error
                rat.tmp <- evidence[,  grep('^Ratio.H.L$', colnames(evidence), ignore.case=T)]
                rat.tmp[nchar(rat.tmp) == 0] <- NA
                #rat.tmp <- log(as.numeric(rat.tmp),2)
                rat.tmp <- as.numeric(rat.tmp)

                mix.error.tab[1, 1] <- round( median(rat.tmp, na.rm=T),2)
                mix.error.tab['N quantified', 1] <- sum(!is.na( rat.tmp) )

                # latex table
                print(xtable(t(mix.error.tab)), size='tiny')

                rm( rat.tmp)
                \n@\n
               \\end{frame}\n", sep="")
           } # end if doublets
            #################################
            # no experiments: triplets
            #################################
             if(silac.type=="triplets"){
               SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{SILAC Mixing Error}\n
                \\alert{Median} of unnormalized \\alert{evidence} ratios on \\alert{raw scale}.

                \n<<mixingerror, echo=F, results=tex>>=
                mix.error.tab <- matrix( 0, nrow=5, ncol=1, dimnames=list(c('M/L', 'H/L', 'H/M', 'N detected', 'N quantified'), 'Ratio.H.L')  )
                mix.error.tab['N detected',1] <- dim(evidence)[1]

                # calculate mixing error
                # H/L
                rat.tmp.hl <- evidence[,  grep('^Ratio.H.L$', colnames(evidence), ignore.case=T)]
                rat.tmp.hl[nchar(rat.tmp.hl) == 0] <- NA
                #rat.tmp.hl <- log(as.numeric(rat.tmp.hl),2)
                rat.tmp.hl <- as.numeric(rat.tmp.hl)

                # M/L
                rat.tmp.ml <- evidence[,  grep('^Ratio.M.L$', colnames(evidence), ignore.case=T)]
                rat.tmp.ml[nchar(rat.tmp.ml) == 0] <- NA
                #rat.tmp.ml <- log(as.numeric(rat.tmp.ml),2)
                rat.tmp.ml <- as.numeric(rat.tmp.ml)

                # H/M
                rat.tmp.hm <- evidence[,  grep('^Ratio.H.M$', colnames(evidence), ignore.case=T)]
                rat.tmp.hm[nchar(rat.tmp.hm) == 0] <- NA
               #rat.tmp.hm <- log(as.numeric(rat.tmp.hm),2)
                rat.tmp.hm <- as.numeric(rat.tmp.hm)

                mix.error.tab['H/L', 1] <- round( median(rat.tmp.hl, na.rm=T),2)
                mix.error.tab['H/M', 1] <- round( median(rat.tmp.hm, na.rm=T),2)
                mix.error.tab['M/L', 1] <- round( median(rat.tmp.ml, na.rm=T),2)
                mix.error.tab['N quantified', 1] <- sum(!is.na( rat.tmp.hl) )

                # latex table
                print(xtable(t(mix.error.tab)), size='tiny')

                rm( rat.tmp.hl, rat.tmp.hm, rat.tmp.ml)
                \n@\n
               \\end{frame}\n", sep="")
           } # end if triplets



        }

        ############################################
        ##   detected and quantified protein groups
        ##
        ## 20141010 disabled since 'Experiment' column
        ## is not contained in the protein groups table
        ## anymore which is used by function 'quantPG'
        ############################################
        #if(!is.null(experiments)){
        if( FALSE ){
           SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{Quantified protein groups}\n
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<quantPG, echo=F, fig=T, height=8, width=16>>=\n
               silac.quant <- getProteinState(proteinGroups, minPep=",min.peptides,", pep.col='",ifelse(MQversion.num >= 11126, "Peptides.", "Peptides..seq.."),"',cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, cex.main=2)

               # determine the number of protein groups quantified in all experiments
               pg.ratmat <- proteinGroups[ , paste('Ratio.H.L', experiments.dot, sep='.')]
               if(!is.null(dim(pg.ratmat)[2])){
                   if(dim(pg.ratmat)[2] > 1){

                          pg.ratmat.ol <- intersect( rownames(pg.ratmat[which(!is.na( pg.ratmat[, 1]) ), ]) ,  rownames(pg.ratmat[which(!is.na( pg.ratmat[, 2]) ), ]) )
                          if(dim(pg.ratmat)[2] > 2){
                              for(i.ol in 3:dim(pg.ratmat)[2])
                                 pg.ratmat.ol <- intersect(pg.ratmat.ol,  rownames(pg.ratmat[which(!is.na( pg.ratmat[, i.ol]) ), ]) )
                          }
                   }
               }
               \n@
                \\end{center}", sep="")

           if(length(experiments) > 1){
              SweaveFile <- paste(SweaveFile,"
                     Protein groups quantified in all \\Sexpr{length(experiments)} experiments: \\Sexpr{length(pg.ratmat.ol)}
                     \\end{frame}\n", sep="")
           } else {
               SweaveFile <- paste(SweaveFile,"\n\\end{frame}\n", sep="")
           }
        } # end if is.null(experiments)

        ############################################
        # SILAC state of evidences
        ############################################
        SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{SILAC state of evidences}\n
               \\begin{center}
               \\setkeys{Gin}{width=0.7\\textwidth}
               \n<<silac_state, echo=F, fig=T, height=8, width=8>>=\n
               if('SILAC.State' %in% colnames(evidence)){
                   sstate <- table(evidence$SILAC.State)
               } else{
                   sstate <- table(evidence$Labeling.State)
               }
               fancyBarplot(as.vector(sstate), names=names(sstate), xlab='SILAC state', ylab='Number of evidences', main='', ylim=c(0, max(sstate)+(max(sstate)*0.35) ), srt=0)
               legend('top', legend=paste( names(sstate), ': ', round(100*(sstate/sum(sstate)),1), ' %', sep=''), bty='n' )
               \n@
                \\end{center}
                \\end{frame}\n", sep="")



        SweaveFile = paste(SweaveFile, "\\section{Regulated Protein Groups}\n", sep="")
        fig.ratio <- ifelse(norm.ratio, "Normalized ratio", "Ratio")

        ############################################
        #   significantly regulated protein groups
        ############################################
        if(!is.null(experiments) ){

            ##########################################
            # determine a common xlim for all plots
            ##########################################
            SweaveFile = paste(SweaveFile, "
                       \n<<xlim_pg, echo=F>>=\n
                         if(norm.ratio){
                             rat.mat <- proteinGroups[, colnames(proteinGroups)[grep(paste('^Ratio.[H|M].[L|M].Normalized.(',paste( experiments.dot, collapse='|', sep=''), ')$',sep=''), colnames(proteinGroups), ignore.case=T)]]
            } else {
                 rat.mat <- proteinGroups[, colnames(proteinGroups)[grep(paste('^Ratio.[H|M].[L|M].(',paste( experiments.dot, collapse='|', sep=''), ')$',sep=''), colnames(proteinGroups), ignore.case=T)]]
            }
            xlim.rat.int =   range(log(unlist(rat.mat),2), na.rm=T)
            xlim.rat.int = c(-max(abs(xlim.rat.int)),  max(abs(xlim.rat.int)))
                       \n@\n", sep="")


            ##########################################
            # loop over the experiments
            ##########################################

            for(ex in experiments){

               SweaveFile = paste(SweaveFile, "\\begin{frame}
                    \\frametitle{Protein groups: experiment ",experiments.tex[ex],"}\n

                     \\begin{center}

                      \\setkeys{Gin}{width=0.6\\textwidth}
                       \n<<ratioVsint-",experiments.tex[ex],", echo=F, fig=T, height=7, width=7>>=\n
                       tmp <- sigRatioPlot.ph(table=proteinGroups,  experiment='",experiments.dot[ex],"',  norm=",norm.ratio,", pcut=", psig,", main='Heavy to light', silac.state='H.L', adjust='", padjust,"', xlim= xlim.rat.int, reset.par=T )
                       \n@\n
                     \\end{center}

                     \\end{frame}\n", sep="")

               if(silac.type=="triplets"){

                  # heavy to medium
                  SweaveFile = paste(SweaveFile, "\\begin{frame}
                      \\frametitle{Protein groups: experiment ",experiments.tex[ex],"}\n

                         \\begin{center}
                           \\setkeys{Gin}{width=0.6\\textwidth}
                            \n<<ratioVsint-",experiments.tex[ex],"hm, echo=F, fig=T, height=7, width=7>>=\n
                            tmp <- sigRatioPlot.ph(table=proteinGroups, experiment='",experiments.dot[ex],"',  norm=",norm.ratio,", pcut=", psig,", main='Heavy to medium', silac.state='H.M', adjust='",padjust,"', xlim= xlim.rat.int, reset.par=T)
                            \n@
                             \n
                           \\end{center}
                            \\end{frame}\n", sep="")
                  # medium to light
                  SweaveFile = paste(SweaveFile, "\\begin{frame}
                      \\frametitle{Protein groups: experiment ",experiments.tex[ex],"}\n

                         \\begin{center}
                           \\setkeys{Gin}{width=0.6\\textwidth}
                            \n<<ratioVsint-",experiments.tex[ex],"ml, echo=F, fig=T, height=7, width=7>>=\n
                            tmp <- sigRatioPlot.ph(table=proteinGroups,  experiment='",experiments.dot[ex],"',  norm=",norm.ratio,", pcut=",psig,", main='Medium to light', silac.state='M.L', adjust='",padjust,"', xlim= xlim.rat.int, , reset.par=T)
                            \n@
                            \n
                            \\end{center}
                            \\end{frame}\n", sep="")

                  } # end if, silac.state

            } # end for, experiments

        ##############################################################
        #    if there are no experiments defined ...
        ##############################################################
        } else {

            ########################################
            # determine common x-lim
            ########################################
             SweaveFile = paste(SweaveFile, "
                       \n<<xlim_pg, echo=F>>=\n
                         if(norm.ratio){
                               rat.mat <- proteinGroups[, colnames(proteinGroups)[grep(paste('^Ratio.[H|M].[L|M].Normalized',sep=''), colnames(proteinGroups), ignore.case=T)]]
                          } else {
                            rat.mat <- proteinGroups[, colnames(proteinGroups)[grep(paste('^Ratio.[H|M].[L|M]$',sep=''), colnames(proteinGroups), ignore.case=T)]]
                          }
                       xlim.rat.int =   range(log(unlist(rat.mat),2), na.rm=T)
                       xlim.rat.int = c(-max(abs(xlim.rat.int)),  max(abs(xlim.rat.int)))
                       \n@\n", sep="")


          SweaveFile = paste(SweaveFile, "\\begin{frame}
                    \\frametitle{Protein groups}\n

                     \\begin{center}
                      \\setkeys{Gin}{width=0.6\\textwidth}
                       \n<<ratioVsint, echo=F, fig=T, height=7, width=7>>=\n
                       tmp <- sigRatioPlot.ph(table=proteinGroups, experiment='', norm=",norm.ratio,", pcut=",psig,", main='Heavy to light', silac.state='H.L', adjust='",padjust,"', xlim= xlim.rat.int, reset.par=T)
                       \n@
                       \n
                     \\end{center}
                     \\end{frame}\n", sep="")

               if(silac.type=="triplets"){
                  # heavy to medium
                  SweaveFile = paste(SweaveFile, "\\begin{frame}
                      \\frametitle{Protein groups}\n

                         \\begin{center}
                           \\setkeys{Gin}{width=0.6\\textwidth}
                            \n<<ratioVsint-hm, echo=F, fig=T, height=7, width=7>>=\n
                            tmp <- sigRatioPlot.ph(table=proteinGroups,  experiment='',  norm=",norm.ratio,", pcut=", psig,", main='Heavy to medium', silac.state='H.M', adjust='",padjust,"', xlim= xlim.rat.int, reset.par=T)
                            \n@
                            \\end{center}
                            \\end{frame}\n", sep="")
                  # medium to light
                  SweaveFile = paste(SweaveFile, "\\begin{frame}
                      \\frametitle{Protein groups}\n

                         \\begin{center}
                           \\setkeys{Gin}{width=0.6\\textwidth}
                            \n<<ratioVsint-ml, echo=F, fig=T, height=7, width=7>>=\n
                            tmp <- sigRatioPlot.ph(table=proteinGroups,  experiment='', norm=",norm.ratio,", pcut=",psig,", main='Medium to light', silac.state='M.L', adjust='", padjust,"', xlim= xlim.rat.int, reset.par=T)
                            \n@

                            \\end{center}
                            \\end{frame}\n", sep="")

                  } # end if, silac.state
        }
    }
    ###############################################################################
    #
    #               regulated sites: ratio vs intensity
    #
    ###############################################################################
    if(quant & (length(siteTabs) > 0)){
        for(s in names(siteTabs))
            SweaveFile <- paste(SweaveFile, summaryPDF.Sites( s, what="ratios",  mod.header=mod.header.names[s], phLscore=phLscore, phPEP=phPEP, quant=quant, experiments=experiments, experiments.tex=experiments.tex, experiments.dot = experiments.dot, silac.type=silac.type, psig=psig, padjust=padjust ), sep="")
    }


    #####################################################
    #
    #         incorporation check
    #
    #####################################################
    if(ic ){

        SweaveFile <- paste(SweaveFile,
                               "\\section{Incorporation check}\n", sep="")

        # no experimental design
        if(is.null(experiments)){

            SweaveFile <- paste(SweaveFile,
                               "\\begin{frame}
                                \\frametitle{ Inorporation Check}\n
                                Using \\texttt{evidence.txt} table.
                                  \\begin{center}
                                    \\setkeys{Gin}{width=0.95\\textwidth}
                                    \n<<IC, echo=F, fig=T, height=7, width=14>>=\n

                                     tmp <- checkIncorporation( evidence )

                                  \n@\n
                                  \\end{center}
                                \n\\end{frame}\n",sep="")
        } else {


            # loop over the experiments
            for(ex in experiments){

                  SweaveFile <- paste(SweaveFile,
                               "\\begin{frame}
                                \\frametitle{ Inorporation Check (",experiments.tex[ex],")}\n
                                Using \\texttt{ evidence.txt} table.
                                  \\begin{center}
                                    \\setkeys{Gin}{width=0.95\\textwidth}
                                    \n<<IC_",ex,", echo=F, fig=T, height=7, width=14>>=\n

                                      # get the peptides detected in the experiment (use the 'Experiment.'-columns)
                                      #ex.idx <- rownames(peptides)[ which(!is.na( peptides[, paste('Experiment.','", experiments.dot[ex],"', sep='' ) ]  ) ) ]
                                      ex.idx <- grep( '^", experiments.dot[ex],"$', evidence[, grep('^Experiment$',colnames(evidence), ignore.case=T )] )
                                      tmp <- checkIncorporation( evidence[ex.idx, ] )

                                  \n@\n
                                  \\end{center}
                                \n\\end{frame}\n",sep="")

            } # end for ex

        } # end else

    } # end if ic



    ##########################################
    # low-level qc
    ##########################################
    if(low.level.qc){


       #################################################################
       #
       #                 number of MS and MSMS scans
       #
       #################################################################
       SweaveFile <- paste(SweaveFile, "\\section{MS and MS/MS Scans }
                   \n\n",sep="")



        # if there are no experiments defined...
        if(is.null(experiments)){

            SweaveFile <- paste(SweaveFile, "
              \\begin{frame}
               \\frametitle{MS Scans}
                   \\begin{center}
                   \\setkeys{Gin}{width=0.7\\textwidth}\n
                   \n<<msScansBarplot, echo=F, fig=TRUE, height=10, width=15>>=

                      # MS scans per raw files
                      msScansAll <- table( msScans[ grep('^Raw.File$', colnames(msScans), ignore.case=T )])

                      # MSMS
                      msScansAll.MSMS <- vector('numeric', length(msScansAll))
                      names(msScansAll.MSMS) <- names(msScansAll)
                      for(r in names(msScansAll))
                            msScansAll.MSMS[r] <- sum(msScans[ which(msScans[, grep('^Raw.File$', colnames(msScans), ignore.case=T )] == r), grep('^MS.MS.Count$', colnames(msScans), ignore.case=T )], na.rm=T )


                      par(mar=c( 20,4,4,2) )
                      fancyBarplot( rbind(msScansAll, msScansAll.MSMS), main='Number of MS and MS/MS scans', las=2, cex.main=1.5)

                      par(mar=c(5,4,4,2))
                   \n@
                   \\end{center}

               \\end{frame}\n\n",sep="")

        } else { # end if experiments

           for(ex in experiments){


                SweaveFile <- paste(SweaveFile, "\n
                  \\begin{frame}
                   \\frametitle{MS Scans: ", experiments.tex[ex],"}
                       \\begin{center}
                       \\setkeys{Gin}{width=0.7\\textwidth}\n
                       \n<<msScansBarplot",experiments.tex[ex] ,", echo=F, fig=TRUE, height=10, width=15>>=
                          #ex.idx <- which( msScans[, grep('^Experiment$', colnames(msScans), ignore.case=T )] == '",experiments.dot[ex],"' )
                          ex.idx <- grep( '^",experiments.dot[ex],"$', msScans[, grep('^Experiment$', colnames(msScans), ignore.case=T )])

                          msScansEx <- msScans[ex.idx, ]
                          # MS scans per raw file
                          msScansEx.MS <- table( msScansEx[ ,grep('^Raw.File$', colnames(msScansEx), ignore.case=T)] )

                          # MSMS scans per raw file
                          msScansEx.MSMS <- vector('numeric', length(msScansEx.MS))
                          names(msScansEx.MSMS) <- names(msScansEx.MS)
                          for(r in names(msScansEx.MS))
                               msScansEx.MSMS[r] <- sum(msScansEx[ which( msScansEx[ , grep('^Raw.File$', colnames(msScansEx), ignore.case=T ) ] == r), grep('^MS.MS.Count$', colnames(msScansEx), ignore.case=T )], na.rm=T )

                          par(mar=c( 20,4,4,2))

                          fancyBarplot(rbind(msScansEx.MS, msScansEx.MSMS), main='Number of MS and MS/MS scans', las=2, cex.main=1.5)

                          par(mar=c(5,4,4,2))
                          rm(msScansEx, ex.idx)
                       \n@
                       \\end{center}
                   \\end{frame}\n\n",sep="")



            } # end for experiments
        }

       #################################################################
       #
       #                 ion injection time
       #
       #################################################################
       SweaveFile <- paste(SweaveFile, "\\section{Ion Injection times }
                   \n\n",sep="")



        # if there are no experiments defined...
        if(is.null(experiments)){

            SweaveFile <- paste(SweaveFile, "
              \\begin{frame}
               \\frametitle{MS Ion Injection Times}
                   \\begin{center}
                   \\setkeys{Gin}{width=0.7\\textwidth}\n
                   \n<<msScansIonInjection, echo=F, fig=TRUE, height=10, width=15>>=

                      # ion injection times per raw file
                      msScans.ioninjection <- vector('list', length(msScansAll))
                      names(msScans.ioninjection) <- names(msScansAll)
                      for(r in names(msScans.ioninjection))
                            msScans.ioninjection[[r]] <- msScans[which(msScans[, grep('^Raw.File$', colnames(msScans), ignore.case=T )] == r), grep('^Ion.Injection.Time$', colnames(msScans), ignore.case=T ) ]


                      par(mar=c( 20,4,4,2))
                      #X11()
                      boxplot( msScans.ioninjection, main='Ion injection time', las=2, cex.main=1.5)
                      par(mar=c(5,4,4,2))
                   \n@\n\n
                   \\end{center}

               \\end{frame}\n\n",sep="")

        } else { # end if experiments

           for(ex in experiments){


                SweaveFile <- paste(SweaveFile, "\n
                  \\begin{frame}
                   \\frametitle{MS Ion Injection Times: ", experiments.tex[ex],"}
                       \\begin{center}
                       \\setkeys{Gin}{width=0.7\\textwidth}\n
                       \n<<msScansIonInjection",experiments.tex[ex] ,", echo=F, fig=TRUE, height=10, width=15>>=
#                         ## ex.idx <- which(msScans[, grep('^Experiment$', colnames(msScans), ignore.case=T )] == '",experiments.dot[ex],"' )
                          ex.idx <- grep( '^",experiments.dot[ex],"$', msScans[, grep('^Experiment$', colnames(msScans), ignore.case=T )])
                          msScansEx <- msScans[ex.idx, ]

                          # ion injection times per raw file
                          msScans.ioninjection <- vector('list', length(unique(msScansEx[ ,grep('^Raw.File$', colnames(msScansEx), ignore.case=T) ])))
                          names(msScans.ioninjection) <- unique(msScansEx[, grep('^Raw.File$', colnames(msScansEx), ignore.case=T )])
                          for(r in names(msScans.ioninjection))
                                msScans.ioninjection[[r]] <- msScansEx[ which(msScansEx[ ,grep('^Raw.File$', colnames(msScansEx), ignore.case=T )] == r),  grep('^Ion.Injection.Time$', colnames(msScansEx), ignore.case=T )]


                          par(mar=c( 20,4,4,2))

                          boxplot( msScans.ioninjection, main='Ion injection time', las=2, cex.main=1.5)
                          par(mar=c(5,4,4,2))

                       \n@\n\n
                       \\end{center}
                   \\end{frame}\n\n",sep="")



            } # end for experiments
        } # end else

       #################################################################
       #
       #                 cycle time
       #
       #################################################################
       SweaveFile <- paste(SweaveFile, "\\section{Cycle Times }
                   \n\n",sep="")



        # if there are no experiments defined...
        if(is.null(experiments)){

            SweaveFile <- paste(SweaveFile, "
              \\begin{frame}
               \\frametitle{MS Cycle Times}
                   \\begin{center}
                   \\setkeys{Gin}{width=0.7\\textwidth}\n
                   \n<<msScansCycleTime, echo=F, fig=TRUE, height=10, width=15>>=

                      # cycle times per raw file
                      msScans.cycletime <- vector('list', length(msScansAll))
                      names(msScans.cycletime) <- names(msScansAll)
                      for(r in names(msScans.cycletime))
                            msScans.cycletime[[r]] <- msScans[which(msScans[, grep('^Raw.File$', colnames(msScans), ignore.case=T )] == r), grep('^Cycle.Time$', colnames(msScans), ignore.case=T )]

                      par(mar=c( 20,4,4,2))

                      boxplot( msScans.cycletime, main='Cycle time', las=2, cex.main=1.5)
                      par(mar=c(5,4,4,2))
                   \n@\n\n
                   \\end{center}

               \\end{frame}\n\n",sep="")

        } else { # end if experiments

           for(ex in experiments){


                SweaveFile <- paste(SweaveFile, "\n
                  \\begin{frame}
                   \\frametitle{MS Cycle Times: ", experiments.tex[ex],"}
                       \\begin{center}
                       \\setkeys{Gin}{width=0.7\\textwidth}\n
                       \n<<msScansCycleTime",experiments.tex[ex] ,", echo=F, fig=TRUE, height=10, width=15>>=
#                          ex.idx <- which(msScans[,grep('^Experiment$', colnames(msScans), ignore.case=T )] == '",experiments.dot[ex],"' )
                          ex.idx <- grep( '^",experiments.dot[ex],"$', msScans[, grep('^Experiment$', colnames(msScans), ignore.case=T )])
                          msScansEx <- msScans[ex.idx, ]

                          # cycle times per raw file
                          msScans.cycletime <- vector('list', length(unique(msScansEx[,grep('^Raw.File$', colnames(msScansEx), ignore.case=T )])))
                          names(msScans.cycletime) <- unique(msScansEx[, grep('^Raw.File$', colnames(msScansEx), ignore.case=T )])
                          for(r in names(msScans.cycletime))
                                msScans.cycletime[[r]] <- msScansEx[which(msScansEx[, grep('^Raw.File$', colnames(msScansEx), ignore.case=T )] == r), grep('^Cycle.Time$', colnames(msScansEx), ignore.case=T )]


                          par(mar=c( 20,4,4,2))
                          #X11()
                          boxplot( msScans.cycletime, main='Cycle time', las=2, cex.main=1.5)
                          par(mar=c(5,4,4,2))

                       \n@\n\n
                       \\end{center}
                   \\end{frame}\n\n",sep="")



            } # end for experiments
        }

    } # end if qc



    #########################################
    #  used parameters
    #########################################
    SweaveFile <- paste(SweaveFile, "\\section{Parameters}\n",  sep="")
    for(i in 1:nParamPages){

         SweaveFile <- paste(SweaveFile, "\\begin{frame}
                          \\frametitle{Parameters}\n",  sep="")

         SweaveFile <- paste(SweaveFile, "<<parameters, echo=F, results=tex>>=\n
                         print(xtable(param[(1 + (",i,"-1)*20 ):min(dim(param)[1], ",i,"*20 ), ] ), size='tiny', align=c('llc'))\n@\n\n",sep="")

         SweaveFile <- paste(SweaveFile, "\\end{frame}", sep="")
    }



   #SweaveFile <- paste(SweaveFile, summaryPDF.table.desc())
   #SweaveFile <- paste(SweaveFile, "\n\\includepdf[pages={1,2}]{h:/Projects/MaxQuantParser/tables.pdf}")


   ###########################################
   # end of Rnw file
   ###########################################
   SweaveFile <- paste(SweaveFile, "\n\n\\end{document}\n", sep="")


   ###############################################################
   #             execute & compile
   ###############################################################
   dir.create( outDir  )
   setwd(outDir)
   dir.create( "pic"  )

   # write the Sweave file
   cat( SweaveFile, file = paste("summary_",title,".Rnw", sep="") )

   # process code chunks
   require(tools)
   Sweave( paste("summary_",title,".Rnw", sep="") )

   # compile tex file
   texi2dvi(paste("summary_",title,".tex", sep=""), pdf=T)

   ################################################################
   # clean up
   ################################################################
   if(clean.up){

      unlink("pic", recursive=T)
      file.remove( dir()[-grep("(*.txt$|*.pdf$)", dir()) ])

   }
   ################################################################
   # calculate running time
   ################################################################
   dt <- difftime( Sys.time(), time.start, units="secs")
   if(dt <= 60)
       cat("\n\npdf-file produced in ", dt, " seconds\n\n\n" )
   if(dt > 60)
       cat("\n\npdf-file produced in ", dt/60, " minutes\n\n\n" )

   setwd(d)
}


#####################################################################################
#
#  generate latex code for modification site tables
#
#
# mod               - string, modification name, e.g. Phospho, Oxidation
# what              -
# mod.header        - string, name of modification as in the table header, e.g. Phospho..STY.
# phLscore          - numeric, localization score threshold
# phPEP             - numeric, PEP threshold  -> OBSOLETE
# quant             - logical, SILAC or not?
# experiments       - character vector, original experiment names
#
# exeperiments.dot  - character vector, all special characters replaced by '.'
# expermiments.tex  - character vector,
# silac.type        - character, singlets, doublets, triplets
# psig              - numeric, threshold for significance B
# padjust           - character, specifies correction method of significance B values
#
# changelog:  20110331 implementation
#             20110701 significance B values for unnormalized ratios are calculated
#             20110902 H/M ratios are plotted as well
#             20111126 - in case there are no sites to plot, a dummy plot is produced using
#                        the 'scatterplot' function
#                      - before I used 'plot' which sometimes did not work (figure margins too large)
#             20120210 - rat vs. int: if normalized SILAC ratios are used, it is written in the title of the figure
#####################################################################################
summaryPDF.Sites <- function(mod, what=c("general", "ratios"), mod.header, phLscore, phPEP, quant, experiments, experiments.dot, experiments.tex, silac.type="singlets", psig, padjust, norm.ratio){

      #######################################################
      # determine what kind of tex code should be produced
      #######################################################
      what = match.arg(what)

      ####################################################################################
      #
      #                   general numbers and ovelap
      #
      ####################################################################################
      if(what == "general"){
            ############################################
            # non-redundant unlocalized sites:
            #
            ############################################
            SweaveFile = paste( "\\section{",mod," Table}\n\\begin{frame}
               \\frametitle{", mod," Table}\n

               \\begin{description}
                   \\item[Localized sites:] localization probability $\\geq$ ",phLscore," \\\\
               \\end{description}

               Site redundancy:

               \\begin{center}
               \n<<",mod,"NRsites, echo=F, fig=F, results=tex, height=7, width=14>>=\n

               #######################
               # unlocalized sites
               #######################
               ",mod,"Sites.unloc <- ",mod,"Sites[ ",mod,"Sites[, grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )] < ", phLscore ,"  ,  ]

               # get ids of non-redundant sites
               ",mod,".nr.sites.list.unloc <- getNRSites(",mod,"Sites.unloc, mod='",mod.header,"')

               # get the corresponding part of the site table
               ",mod,"Sites.unloc.nr <- ",mod,"Sites.unloc[ unlist(",mod,".nr.sites.list.unloc), ]


               ##################################################################
               # export that table: localized + non-redundant unlocalized sites
               ##################################################################
               ",mod,"Sites.nr <- rbind(",mod,"Sites.loc, ",mod,"Sites.unloc.nr)

               if(!('", mod,"Sites_nonredundant.txt' %in% dir() ) )
                       write.table(",mod,"Sites.nr, '", mod,"Sites_nonredundant.txt', sep='\\t', quote=F, na='', col.names=NA)


               ###################################
               # put some numbers in a table
               ###################################
               ",mod,"nr.sites.mat <- matrix(0, nrow=3, ncol=2, dimnames=list(c( 'localized', 'unlocalized', 'total'), c('all sites', 'non-redundant sites')  ))

               ",mod,"nr.sites.mat[1,1] <- dim(",mod,"Sites.loc)[1]
               ",mod,"nr.sites.mat[1,2] <- dim(",mod,"Sites.loc)[1]

               ",mod,"nr.sites.mat[2,1] <- dim(",mod,"Sites.unloc)[1]
               ",mod,"nr.sites.mat[2,2] <- dim(",mod,"Sites.unloc.nr)[1]

               ",mod,"nr.sites.mat[3,1] <- ",mod,"nr.sites.mat[1,1] + ",mod,"nr.sites.mat[2,1]
               ",mod,"nr.sites.mat[3,2] <- ",mod,"nr.sites.mat[1,2] + ",mod,"nr.sites.mat[2,2]


               print( xtable(",mod,"nr.sites.mat, digits=c(0,0,0)) , size='tiny')

               \n@
               \\end{center}

                Final table contains all localized and all non-redundant unlocalized sites (\\Sexpr{dim(",mod,"Sites.nr)[1]}).

                %\\begin{center}
                %    \\Sexpr{nrow(",mod,"Sites.loc)} + \\Sexpr{nrow(",mod,"Sites.unloc.nr)} = \\Sexpr{dim(",mod,"Sites.nr)[1]}.
                %\\end{center}

                \\end{frame}\n", sep="")


            ############################################
            # site composition
            # - only if modification can occur on several residues
            ############################################
            if( mod == "Phospho"){

                SweaveFile <- paste(SweaveFile, "\\begin{frame}
                          \\frametitle{", mod, " Table}\n
                          \\begin{description}
                             \\item[Localized sites:] localization probability $\\geq$ ",phLscore," \\\\
                          \\end{description}\n
                          ",  sep="")

                SweaveFile <- paste(SweaveFile, "\n<<",mod,"CompositionTab, echo=F, results=tex>>=\n
                         ",mod,".comp.all <- table(",mod,"Sites[, grep( '^Amino.Acid$', colnames(",mod,"Sites), ignore.case=T) ])

                         ",mod,".comp.L75 <- table(",mod,"Sites.loc[, grep( '^Amino.Acid$', colnames(",mod,"Sites), ignore.case=T)])
                         ",mod,".comp.nr <- table(",mod,"Sites.nr[, grep( '^Amino.Acid$', colnames(",mod,"Sites), ignore.case=T)])

                         # table containing the STY compositions of the dataset filtered differentially
                         ",mod,".comp.tab <- rbind(",mod,".comp.all, ",mod,".comp.L75 )
                         ",mod,".comp.tab <- rbind(",mod,".comp.tab, ",mod,".comp.nr)

                         # add the total numbers
                         ",mod,".comp.tab <- cbind(",mod,".comp.tab, c(sum(",mod,".comp.all), sum(",mod,".comp.L75),  sum(",mod,".comp.nr) )  )

                         # add percentage
                         for(i in 1:dim(",mod,".comp.tab)[1])
                               for(j in 1:(dim(",mod,".comp.tab)[2]-1))
                                        ",mod,".comp.tab[i,j] <- paste(",mod,".comp.tab[i,j], '(', round((as.numeric(",mod,".comp.tab[i,j])/as.numeric(",mod,".comp.tab[i, dim(",mod,".comp.tab)[2]]))*100,2), '%)' )

                         # row- and column names
                         rownames(",mod,".comp.tab) <- c('All sites', 'Localized', 'non-redundant')
                         colnames(",mod,".comp.tab)[ dim(",mod,".comp.tab)[2]  ] <- 'Total'

                         print(xtable( ",mod,".comp.tab ), size='tiny')

                        ############################################################################
                        \n\n@\n\n",sep="")




                SweaveFile <- paste( SweaveFile, "\\end{frame}\n", sep="")
            } # end if Phospho

            ###########################################
            # distribution of localization score
            ############################################
            SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{",mod," Table}\n
               Distribution of localization probabilities.

               \\begin{center}
               \\setkeys{Gin}{width=0.9\\textwidth}
               \n<<",mod,"LOCdist, echo=F, fig=T, height=4, width=9>>=\n
               par(mfrow=c(1,3))
               # all sites
               hist( ",mod,"Sites[, grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )], main='All sites', col='lightblue', cex.axis=1.5, xlab='Localization probability', breaks=40 )
               legend('topleft', legend=c(paste('median:', round(median(",mod,"Sites[,  grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )]),3))  ), bty='n', cex=1.5 )
               # localized sites
               hist( ",mod,"Sites.loc[,  grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )], main='Localized sites', col='lightblue', cex.axis=1.5, xlab='Localization probability', breaks=40 )
               legend('topleft', legend=c(paste('median:', round(median(",mod,"Sites.loc[,  grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )]),3))  ), bty='n', cex=1.5 )

               # non-redundant
               hist( ",mod,"Sites.nr[, grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )], main='Non-redundant sites', col='lightblue', cex.axis=1.5, xlab='Localization probability', breaks=40 )
               legend('topleft', legend=c(paste('median:', round(median(",mod,"Sites.nr[, grep('^Localization.Prob$', colnames(",mod,"Sites), ignore.case=T )]),3))  ), bty='n', cex=1.5 )
               \n@
                \\end{center}
                \\end{frame}\n", sep="")


            ############################################
            # site PEP distribution
            ############################################
            SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{",mod," Table}\n
               Distribution of posterior error probabilities.
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<",mod,"PEPdist, echo=F, fig=T, height=8, width=16>>=\n
               boxplot(",mod,"Sites[, grep('^PEP$', colnames(",mod,"Sites), ignore.case=T)], ",mod,"Sites.loc[ , grep('^PEP$', colnames(",mod,"Sites.loc), ignore.case=T) ], ",mod,"Sites.nr[, grep('^PEP$', colnames(",mod,"Sites.nr), ignore.case=T)], names=c('all sites', 'localized', 'non-redundant'), col='lightblue', ylab='PEP', cex.axis=1.5)
               \n@
                \\end{center}
                \\end{frame}\n", sep="")


            ###############################################
            # PEP distribution: modified vs unmodified peptides
            ###############################################
            SweaveFile = paste(SweaveFile, "\\begin{frame}
               \\frametitle{PEP Distribution: ",mod," vs. Unmodified}\n
               Based on the \\texttt{modifiedPeptides.txt} file.
               \\begin{center}
               \\setkeys{Gin}{width=0.8\\textwidth}
               \n<<",mod,"PEPdistModvsUnmod, echo=F, fig=T, height=8, width=16>>=\n
               PEPunmod <- as.numeric( modifiedPep[grep('^Unmodified', modifiedPep[, grep( '^Modifications$', colnames(modifiedPep), ignore.case=T  )]), grep('^PEP$', colnames(modifiedPep), ignore.case=T)])
               PEP",mod," <- as.numeric(modifiedPep[grep('",mod,"', modifiedPep[, grep('^Modifications$', colnames(modifiedPep), ignore.case=T)]), grep('^PEP$', colnames(modifiedPep), ignore.case=T)])

               boxplot(PEPunmod, PEP",mod,", names=c(paste('unmodifed peptides (', length(PEPunmod), ')'), paste('",mod,"-peptides (', length(PEP",mod,"), ')')), col='lightblue', ylab='PEP', cex.axis=1.5)
               \n@
                \\end{center}
                \\end{frame}\n", sep="")

            ###############################################
            # overlap of sites between experiments
            ###############################################
            if(length(experiments) < 5 & length(experiments) > 1){

                SweaveFile = paste(SweaveFile, "\\begin{frame}
                   \\frametitle{",mod," Table}\n
                   Overlap between experiments: detected sites
                   \\begin{center}
                   \\setkeys{Gin}{width=1.0\\textwidth}
                   \n<<",mod,"Overlap, echo=F, fig=T, height=5, width=15>>=\n
                   # calculate the overlap
                   ol.",mod,".Sites <- siteOverlap(",mod,"Sites, experiments=experiments.dot, phPEP=",phPEP,", phLoc=",phLscore,", quant=", quant,", quantCol.prefix=colnames(",mod,"Sites)[grep('^Ratio.H.L.Normalized$', colnames(",mod,"Sites), ignore.case=T)])
                   ol.",mod,".Sites.nr <- siteOverlap(",mod,"Sites.nr, experiments=experiments.dot, phPEP=",phPEP,", phLoc=",phLscore,", quant=", quant,", quantCol.prefix=colnames(",mod,"Sites)[grep('^Ratio.H.L.Normalized$', colnames(",mod,"Sites), ignore.case=T)])


                   #######################
                   # venn diagram
                   #######################
                   ol.",mod,".all <- ol.",mod,".Sites[[1]]
                   ol.",mod,".loc <- ol.",mod,".Sites[['localized sites']]
                   ol.",mod,".nr <- ol.",mod,".Sites.nr[[1]]


                   if(length(experiments) == 3){
                      par(mfrow=c(1,3))
                      ",mod,".ol.all <- venndiagram( ol.",mod,".all[[experiments.dot[1]]], ol.",mod,".all[[experiments.dot[2]]], ol.",mod,".all[[experiments.dot[3]]], labels=names(ol.",mod,".all), title='all sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )
                      ",mod,".ol.loc <- venndiagram( ol.",mod,".loc[[experiments.dot[1]]], ol.",mod,".loc[[experiments.dot[2]]], ol.",mod,".loc[[experiments.dot[3]]], labels=names(ol.",mod,".loc), title='localized sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )
                      ",mod,".ol.nr <- venndiagram( ol.",mod,".nr[[experiments.dot[1]]], ol.",mod,".nr[[experiments.dot[2]]], ol.",mod,".nr[[experiments.dot[3]]], labels=names(ol.",mod,".nr), title='non-redundant sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )

                   }
                   if(length(experiments) == 2){
                      par(mfrow=c(1,3))
                      ",mod,".ol.all <- venndiagram( ol.",mod,".all[[experiments.dot[1]]], ol.",mod,".all[[experiments.dot[2]]], labels=names(ol.",mod,".all), title='all sites', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )
                      ",mod,".ol.loc <- venndiagram( ol.",mod,".loc[[experiments.dot[1]]], ol.",mod,".loc[[experiments.dot[2]]], labels=names(ol.",mod,".loc), title='localized sites', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )
                      ",mod,".ol.nr <- venndiagram( ol.",mod,".nr[[experiments.dot[1]]], ol.",mod,".nr[[experiments.dot[2]]], labels=names(ol.",mod,".nr), title='non-redundant sites', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )

                   }
                   if(length(experiments) == 4){
                      par(mfrow=c(1,3))
                      ",mod,".ol.all <- venndiagram( ol.",mod,".all[[experiments.dot[1]]], ol.",mod,".all[[experiments.dot[2]]], ol.",mod,".all[[experiments.dot[3]]], ol.",mod,".all[[experiments.dot[4]]], labels=names(ol.",mod,".all), title='all sites', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )
                      ",mod,".ol.loc <- venndiagram( ol.",mod,".loc[[experiments.dot[1]]], ol.",mod,".loc[[experiments.dot[2]]], ol.",mod,".loc[[experiments.dot[3]]], ol.",mod,".loc[[experiments.dot[4]]],labels=names(ol.",mod,".loc), title='localized sites', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )
                      ",mod,".ol.nr <- venndiagram( ol.",mod,".nr[[experiments.dot[1]]], ol.",mod,".nr[[experiments.dot[2]]], ol.",mod,".nr[[experiments.dot[3]]], ol.",mod,".nr[[experiments.dot[4]]],labels=names(ol.",mod,".nr), title='non-redundant', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )

                   }

                   \n@
                   \\end{center}
                   \\end{frame}\n", sep="")

            if(quant){
               ##################################
               # quantified sites
               ##################################
               SweaveFile = paste(SweaveFile, "\\begin{frame}
                   \\frametitle{",mod," Table}\n
                   Overlap between experiments: quantified sites
                   \\begin{center}
                   \\setkeys{Gin}{width=1.0\\textwidth}
                   \n<<",mod,"OverlapQuant, echo=F, fig=T, height=5, width=15>>=\n
                   #######################
                   # venn diagram
                   #######################
                   ol.",mod,".quant <- ol.",mod,".Sites[['all quantified sites']]
                   ol.",mod,".quant.loc <- ol.",mod,".Sites[['quantified localized sites']]
                   ol.",mod,".quant.nr <- ol.",mod,".Sites.nr[['all quantified sites']]

                   if(length(experiments) == 3){
                      par(mfrow=c(1,3))
                      ",mod,".ol.quant <- venndiagram( ol.",mod,".quant[[experiments.dot[1]]], ol.",mod,".quant[[experiments.dot[2]]], ol.",mod,".quant[[experiments.dot[3]]], labels=names(ol.",mod,".quant), title='all quantified sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )
                      ",mod,".ol.quant.loc <- venndiagram( ol.",mod,".quant.loc[[experiments.dot[1]]], ol.",mod,".quant.loc[[experiments.dot[2]]], ol.",mod,".quant.loc[[experiments.dot[3]]], labels=names(ol.",mod,".quant.loc), title='quantified and localized sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )
                      ",mod,".ol.quant.nr <- venndiagram( ol.",mod,".quant.nr[[experiments.dot[1]]], ol.",mod,".quant.nr[[experiments.dot[2]]], ol.",mod,".quant.nr[[experiments.dot[3]]], labels=names(ol.",mod,".quant.nr), title='quantified non-redundant sites', lwd=2, lines=1:3, lcol=1:3, cex=1.3  )
              }
                   if(length(experiments) == 2){
                      par(mfrow=c(1,3))
                      ",mod,".ol.quant <- venndiagram( ol.",mod,".quant[[experiments.dot[1]]], ol.",mod,".quant[[experiments.dot[2]]], labels=names(ol.",mod,".quant), title='all quantified sites', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )
                      ",mod,".ol.quant.loc <- venndiagram( ol.",mod,".quant.loc[[experiments.dot[1]]], ol.",mod,".quant.loc[[experiments.dot[2]]], labels=names(ol.",mod,".quant.loc), title='quantified and localized sites', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )
                      ",mod,".ol.quant.nr <- venndiagram( ol.",mod,".quant.nr[[experiments.dot[1]]], ol.",mod,".quant.nr[[experiments.dot[2]]], labels=names(ol.",mod,".quant.nr), title='quantified non-redundant sites ', lwd=2, lines=1:2, lcol=1:2, cex=1.3, type='2'  )
                   }

                   if(length(experiments) == 4){
                      par(mfrow=c(1,3))
                      ",mod,".ol.quant <- venndiagram( ol.",mod,".quant[[experiments.dot[1]]], ol.",mod,".quant[[experiments.dot[2]]], ol.",mod,".quant[[experiments.dot[3]]], ol.",mod,".quant[[experiments.dot[4]]], labels=names(ol.", mod,".quant), title='all quantified sites', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )
                      ",mod,".ol.quant.loc <- venndiagram( ol.",mod,".quant.loc[[experiments.dot[1]]], ol.",mod,".quant.loc[[experiments.dot[2]]], ol.",mod,".quant.loc[[experiments.dot[3]]],ol.",mod,".quant.loc[[experiments.dot[4]]],labels=names(ol.",mod,".quant.loc), title='quantified and localized sites', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )
                      ",mod,".ol.quant.nr <- venndiagram( ol.",mod,".quant.nr[[experiments.dot[1]]], ol.",mod,".quant.nr[[experiments.dot[2]]], ol.",mod,".quant.nr[[experiments.dot[3]]], ol.",mod,".quant.nr[[experiments.dot[4]]], labels=names(ol.",mod,".quant.nr), title='quantified non-redundant sites', lwd=2, lines=1:4, lcol=1:4, cex=1.3, type='4'  )
                   }

                   \n@
                   \\end{center}
                   \\end{frame}\n", sep="")

               }
            }
        } # end if what

      ######################################################################################
      #
      #                        ratio vs. intensity plots
      #
      #######################################################################################
      if(what == "ratios"){


          SweaveFile = paste("\\section{Regulated ",mod," Sites}\n", sep="")

             if(!is.null(experiments)){
                ##########################################
                # determine a common xlim for all plots
                ##########################################
                SweaveFile = paste(SweaveFile, "
                       \n<<xlim_",mod,", echo=F>>=\n
                         if(norm.ratio){
                             rat.mat.",mod," <- ",mod,"Sites.nr[, colnames(",mod,"Sites.nr)[grep(paste('^Ratio.[H|M].[L|M].Normalized.(',paste( experiments.dot,collapse='|', sep=''), ')$',sep=''), colnames(",mod,"Sites.nr), ignore.case=T)]]
                         } else {
                             rat.mat.",mod," <- ",mod,"Sites.nr[, colnames(",mod,"Sites.nr)[grep(paste('^Ratio.[H|M].[L|M].(',paste( experiments.dot, collapse='|', sep=''), ')$',sep=''), colnames(",mod,"Sites.nr), ignore.case=T)]]
            }
                          xlim.rat.int.",mod," =   range(log(unlist(rat.mat.",mod,"),2), na.rm=T)
                          xlim.rat.int.",mod," = c(-max(abs(xlim.rat.int.",mod,")),  max(abs(xlim.rat.int.",mod,")))

                         #############################################
                         # check whether ther are any quantified sites
                         #############################################


                       \n@\n", sep="")


                ##########################################
                # loop over the experiments
                ##########################################
                for(ex in experiments){

                    ####################################
                    # medium to light
                    ####################################
                    if(silac.type=="triplets"){

                        SweaveFile = paste(SweaveFile, "\\begin{frame}
                        \\frametitle{",mod," sites (nr): experiment ",experiments.tex[ex],"}\n

                        \\begin{center}

                          \\setkeys{Gin}{width=0.6\\textwidth}
                           \n<<",mod,"ratioVsint-ml-",experiments.tex[ex],", echo=F, fig=T, height=7, width=7>>=\n
                           # get ratios, intensities and significance values
                           if(norm.ratio){
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.M.L.Normalized.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T ) ]
                           } else {
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.M.L.",experiments.dot[ex],"$',colnames(",mod,"Sites.nr), ignore.case=T  )]
                           }
                           int <- ",mod,"Sites.nr[, grep( '^Intensity.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T )]
                           if(norm.ratio & ('Ratio.M.L.Significance.B..",experiments.dot[ex],"' %in% colnames(",mod,"Sites.nr)) ){
                                 p <- ",mod,"Sites.nr[, 'Ratio.M.L.Significance.B..",experiments.dot[ex],"']
                           } else {
                                 p <- significanceB( rat, int  )[[1]]

                                 # if there are only 2 point the calculation of significance values fails...
                                 if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                           }
                           # plot, if at least one site is quantified
                           if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                               tmp <- sigRatioPlot.ph(rat,  int,  p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Medium to light (normalized)',  'Medium to light'),  xlim=xlim.rat.int.",mod,", reset.par=T)
                           } else {
                                 scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                           }
                           rm(rat, int, p)
                           \n@
                         \\end{center}

                         \\end{frame}\n", sep="")

                   ############################################
                   #  heavy to medium
                   ############################################
                   SweaveFile = paste(SweaveFile, "\\begin{frame}
                        \\frametitle{",mod," sites (nr): experiment ",experiments.tex[ex],"}\n

                         \\begin{center}

                          \\setkeys{Gin}{width=0.6\\textwidth}
                           \n<<",mod,"ratioVsint-hm-", experiments.tex[ex],", echo=F, fig=T, height=7, width=7>>=\n
                           # get ratios, intensities and significance values
                           if(norm.ratio){
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.M.Normalized.",experiments.dot[ex],"$',colnames(",mod,"Sites.nr), ignore.case=T )]
                           } else {
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.M.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T)]
                           }
                           int <- ",mod,"Sites.nr[, grep( '^Intensity.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T)]
                           if(norm.ratio & ('Ratio.H.M.Significance.B..",experiments.dot[ex],"' %in% colnames(",mod,"Sites.nr)) ){
                                 p <- ",mod,"Sites.nr[, 'Ratio.H.M.Significance.B..",experiments.dot[ex],"']
                           } else {
                                 p <- significanceB( rat, int  )[[1]]

                                 # if there are only 2 point the calculation of significance values fails...
                                 if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                           }

                           # plot, if at least one site is quantified
                           if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                                           tmp <- sigRatioPlot.ph(rat, int, p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Heavy to Medium (normalized)',  'Heavy to Medium'), xlim=xlim.rat.int.",mod,", reset.par=T)
                            } else {
                                scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                           }
                           rm(rat, int, p)
                           \n@
                         \\end{center}

                         \\end{frame}\n", sep="")


                    }
                   ############################################
                   #  heavy to light
                   ############################################
                   SweaveFile = paste(SweaveFile, "\\begin{frame}
                        \\frametitle{",mod," sites (nr): experiment ",experiments.tex[ex],"}\n

                         \\begin{center}

                          \\setkeys{Gin}{width=0.6\\textwidth}
                           \n<<",mod,"ratioVsint-hl-",experiments.tex[ex],", echo=F, fig=T, height=7, width=7>>=\n
                           # get ratios, intensities and significance values
                           if(norm.ratio){
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.L.Normalized.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T) ]
                           } else {
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.L.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T)]
                           }
                           int <- ",mod,"Sites.nr[, grep( '^Intensity.",experiments.dot[ex],"$', colnames(",mod,"Sites.nr), ignore.case=T) ]
                           if(norm.ratio & ('Ratio.H.L.Significance.B..",experiments.dot[ex],"' %in% colnames(",mod,"Sites.nr)) ){
                                 p <- ",mod,"Sites.nr[, 'Ratio.H.L.Significance.B..",experiments.dot[ex],"']
                           } else {
                                 p <- significanceB( rat, int  )[[1]]

                                 # if there are only 2 point the calculation of significance values fails...
                                 if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                           }

                           # plot, if at least one site is quantified
                           if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                                     tmp <- sigRatioPlot.ph(rat, int, p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Heavy to light (normalized)',  'Heavy to light'), xlim=xlim.rat.int.",mod,", reset.par=T)
                            } else {
                                scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                           }
                           rm(rat, int, p)
                           \n@
                         \\end{center}

                         \\end{frame}\n", sep="")


                } # end for, experiments

            } else{
                #############################################################################
                # if there are no experiments defined...
                #############################################################################

                ########################################
                # determine common x-lim
                ########################################
                 SweaveFile = paste(SweaveFile, "
                        \n<<xlim_",mod,", echo=F>>=\n
                         if(norm.ratio){
                               rat.mat.",mod," <- ",mod,"Sites.nr[, colnames(",mod,"Sites.nr)[grep(paste('^Ratio.[H|M].[L|M].Normalized$',sep=''), colnames(",mod,"Sites.nr), ignore.case=T) ]]
                          } else {
                            rat.mat.",mod," <- ",mod,"Sites.nr[, colnames(",mod,"Sites.nr)[grep(paste('^Ratio.[H|M].[L|M]$',sep=''), colnames(",mod,"Sites.nr), ignore.case=T)]]
                          }
                       xlim.rat.int.",mod," = range(log(rat.mat.",mod,",2), na.rm=T)
                       xlim.rat.int.",mod," = c(-max(abs(xlim.rat.int.",mod,")),  max(abs(xlim.rat.int.",mod,")))

                       \n@\n", sep="")

                ##############################
                # medium to light
                ##############################
                if(silac.type=="triplets"){

                     SweaveFile = paste(SweaveFile, "\\begin{frame}
                        \\frametitle{",mod," sites (nr)}\n

                         \\begin{center}

                          \\setkeys{Gin}{width=0.6\\textwidth}
                           \n<<",mod,"ratioVsint-ml, echo=F, fig=T, height=7, width=7>>=\n
                           if(norm.ratio){
                                   rat <- ",mod,"Sites.nr[, grep('^Ratio.M.L.Normalized$', colnames(",mod,"Sites.nr), ignore.case=T)]
                           } else{
                                   rat <- ",mod,"Sites.nr[, grep('^Ratio.M.L$', colnames(",mod,"Sites.nr), ignore.case=T)]
                           }
                           int <- ",mod,"Sites.nr[, grep('^Intensity$', colnames(",mod,"Sites.nr), ignore.case=T) ]
                           if( norm.ratio & ( 'Ratio.M.L.Significance.B..' %in% colnames(",mod,"Sites.nr)) ){
                                     p <- ",mod,"Sites.nr[, 'Ratio.M.L.Significance.B..']
                           } else {
                                     p <- significanceB( rat, int  )[[1]]

                                     # if there are only 2 point the calculation of significance values fails...
                                     if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                           }
                           # plot, if at least one site is quantified
                           if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                                 tmp <- sigRatioPlot.ph(rat,  int,  p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Medium to light (normalized)',  'Medium to light'), xlim=xlim.rat.int.",mod,", reset.par=T)
                           } else {
                                 scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                           }
                           rm(rat, int, p)
                           \n@
                         \\end{center}

                         \\end{frame}\n", sep="")


                ###############################
                # heavy to medium
                ###############################
                SweaveFile = paste(SweaveFile, "\\begin{frame}
                     \\frametitle{",mod," sites (nr)}\n

                     \\begin{center}

                      \\setkeys{Gin}{width=0.6\\textwidth}
                       \n<<",mod,"ratioVsint-hm, echo=F, fig=T, height=7, width=7>>=\n
                       if(norm.ratio){
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.M.Normalized$', colnames(",mod,"Sites.nr), ignore.case=T)]
                       } else {
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.M$',colnames(",mod,"Sites.nr), ignore.case=T)]
                       }

                       int <- ",mod,"Sites.nr[, grep('^Intensity$',colnames(",mod,"Sites.nr), ignore.case=T)]
                       if('Ratio.H.M.Significance.B..' %in% colnames(",mod,"Sites.nr)){
                                 p <- ",mod,"Sites.nr[, 'Ratio.H.M.Significance.B..']
                       } else {
                                 #p <- rep(1, length(rat))
                                 p <- significanceB( rat, int  )[[1]]

                                 # if there are only 2 point the calculation of significance values fails...
                                 if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                       }

                       # plot, if at least one site is quantified
                       if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                               tmp <- sigRatioPlot.ph( rat, int,  p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Heavy to medium (normalized)',  'Heavy to medium'), xlim=xlim.rat.int.",mod,", reset.par=T)
                       } else {
                             scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                             #plot(1:1e3,1:1e3, col='white', axes=F, xlab='', ylab='')
                       }
                       rm(rat, int, p)
                      # dev.off()
                       \n@
                     \\end{center}

                     \\end{frame}\n", sep="")

                }
                ###############################
                # heavy to light
                ###############################
                SweaveFile = paste(SweaveFile, "\\begin{frame}
                     \\frametitle{",mod," sites (nr)}\n

                     \\begin{center}

                      \\setkeys{Gin}{width=0.6\\textwidth}
                       \n<<",mod,"ratioVsint-hl, echo=F, fig=T, height=7, width=7>>=\n
                       if(norm.ratio){
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.L.Normalized$', colnames(",mod,"Sites.nr), ignore.case=T)]
                       } else {
                                  rat <- ",mod,"Sites.nr[, grep( '^Ratio.H.L$', colnames(",mod,"Sites.nr), ignore.case=T)]
                       }

                       int <- ",mod,"Sites.nr[, grep( '^Intensity$', colnames(",mod,"Sites.nr), ignore.case=T)]
                       if('Ratio.H.L.Significance.B..' %in% colnames(",mod,"Sites.nr)){
                                 p <- ",mod,"Sites.nr[, 'Ratio.H.L.Significance.B..']
                       } else {
                                 #p <- rep(1, length(rat))
                                 p <- significanceB( rat, int  )[[1]]

                                 # if there are only 2 point the calculation of significance values fails...
                                 if(sum(is.na(p)) == length(p))
                                               p <- rep(1, length(p))

                       }

                       # plot, if at least one site is quantified
                       if( length( intersect(which(!is.na(rat)), which(!is.na(int)) )) > 0 ){
                               tmp <- sigRatioPlot.ph( rat, int,  p, pcut=", psig,", adjust='", padjust,"', main=ifelse(norm.ratio, 'Heavy to light (normalized)',  'Heavy to light'), xlim=xlim.rat.int.",mod,", reset.par=T)
                       } else {
                                scatterplot(1,1, grid=F, col='white', axes=F, boxplots='n')
                                #plot(1:1e3,1:1e3, col='white', axes=F, xlab='', ylab='')
                       }
                       rm(rat, int, p)
                       #dev.off()
                       \n@
                     \\end{center}

                     \\end{frame}\n", sep="")

        } # end else

      } # end what

    return(SweaveFile)
}



#####################################################################################
#
#                          add q-values and remove CON/REV
#
# table
# q
# con
# rev
# base.on
# rev.col
# con.col
#
# changelog    20100927 implementation
#              20120427 if column  "Only.identified.by.site" is presented
#                       entries marked with '+' are filterd out as well
#####################################################################################
rmConRev <- function(table, con=T, rev=T, q=T, base.on="PEP", rev.col="Reverse", con.col="Contaminant", oibs= "Only.identified.by.site"){

    # add q-values
    if(q)
        table <- q.value(table, base.on=base.on, rev.col=rev.col)

    # remove CON
    if(con){
        con.idx <- grep("\\+", table[, grep(paste("^", con.col,"$", sep=""), colnames(table), ignore.case=T)] )
        if(length(con.idx) > 0)
            table <- table[-con.idx, ]
    }
    # remove REV
    if(rev){
        rev.idx <- grep("\\+", table[, grep(paste("^", rev.col,"$", sep=""), colnames(table), ignore.case=T)] )
        if(length(rev.idx) > 0)
            table <- table[-rev.idx, ]

    }
    # remove  "Only.identified.by.site"
    if( oibs %in% colnames(table)){
        oibs.idx <- grep("\\+", table[, oibs ])
        if(length(oibs.idx) > 0)
            table <- table[-oibs.idx,]
    }

    return(table)
}


#####################################################
#
#                   fancyBarplot
#
# 20110217 added support for matrices
# 20140201 'cex.numb'
##
#####################################################
fancyBarplot <- function(x, space=0.2, ylim=c(0, max(x, na.rm=T)+ 0.15*max(x, na.rm=T)), ndigits=3, add.numb=T, cex.numb=.9, srt=90, ...)
{
    #########################################
    #         vector
    #########################################
    if(length(dim(x)) <= 1){
        # make the plot
        barplot(x, space=space, ylim=ylim, ...)

        # add numbers
        if(add.numb) text(  seq( 0.5+space, (0.5+space)+((length(x)-1)*(1+space)), 1+space ),
                 x+(0.05*max(x)), round(x, ndigits), srt=srt, pos=4, offset=-0.05, cex=cex.numb )
    } else{

    #########################################
    #        matrix
    #########################################
        offset = -0.1
        barplot(x, ylim=ylim, beside=T, ...)

        if(add.numb){
            p.tmp=1
            for(p in 1:dim(x)[2]){
                if( sum(is.na(x[, p])) < nrow(x) ){
                    text( (p.tmp : (p.tmp + dim(x)[1]-1) ) + offset, max(x[, p], na.rm=T)+0.05*( max(x[, p], na.rm=T)), round(x[,p], ndigits), adj=c(1,1), pos=4, cex=cex.numb, srt=srt  )
                }
                p.tmp <- p.tmp+dim(x)[1]+1
            }
        }

    }
}

######################################################
#
#              fancyDensPlot
#
######################################################
fancyDensPlot <- function(x, breaks=60, hist.col="cyan",  grid=F, ...)
{

    hist(x, breaks=breaks, add=F, freq=F, col=hist.col, ... )
    lines(density(x, na.rm=T),  lwd=2)
    if(grid) grid()
}


################################################################################################
#
# determine in how many raw files a certain site has been detected
#
#
# changelog: 20120207 implementation
################################################################################################
rawfiles.per.site <- function( sites, evidence, mod="Phospho..STY.."  ){

    # get all evidences having at least one site
    #evidence.sites <- evidence[ which(nchar( evidence[, paste( mod, "Site.IDs", sep="")]) > 0 ) ,  ]

    # site ids
    site.ids <- rownames(sites)

    #  get all sites in 'evidence'-table
    #evidence.site.ids.list <- strsplit( evidence.sites[, paste( mod, "Site.IDs", sep="") ], ";")
    #names(evidence.site.ids.list) <- rownames(evidence.sites)

    # raw files
    #rf <- unique(evidence$Raw.File)

    # get only sites present in 'sites'-table
    #evidence.site.ids <- lapply(evidence.site.ids.list, function(x) x[which(x %in% site.ids)])

    # list to store the raw files per detected site
    site.rf <-vector("list", dim(sites)[1])
    names(site.rf) <- site.ids

    #
    for(s in site.ids)
        site.rf[[s]] <- unique(evidence[ grep(paste("(^|;)", s, "($|;)", sep=""), evidence[, paste( mod, "Site.IDs", sep="") ])  , "Raw.File"])

    # get all sites that have been detected in only one LC-MS
    site.rf.numb <- unlist(lapply(site.rf, length))

   # site.rf.unique <- site.rf[which(site.rf.numb == 1)]

    return(site.rf.numb)
   # return(evidence.site.ids)
}

################################################################################################
#                                ratio vs intensity plot
#
# rat            - numeric, ratios on raw scale
# int            - numeric, intensities on raw scale
# p              - numeric, significance values
#
# table          - NULL or a MaxQuant table. currently only the protein groups table is supported
#                - if specified, the rat/int/p vectors are not used but the respective
#                  columns are extracted from that table
# experiment     - character, name of an experiment                                  ## only considered if 'table' != NULL
# silac.state    - character, specifies which SILAC ratio to plot                    ## only considered if 'table' != NULL
# norm           - logical, if TRUE the normalized ratios are used                   ## only considered if 'table' != NULL
# signif         - character, specifies which signifiance values are used            ## only considered if 'table' != NULL
#
# label          - character vector of same length as rat/int/p
# label.string   - character, the string is used to search 'label' via grep
# pcut           - numeric, cutoff of p
# adjust         - character, method of p-value adjustment, "none" means no correction
# col            - charater, color used for unsignificant ratios
# sig.col        - character, color used for significant ratios
# label.col      - character, color used for the label
# pch            - numeric vector, specifying the plotting symbols for unregulated and regulated features
# label.ph       - numeric, plotting symbol for the label
# alpha          - numeric, transparency of points
# boxplots       - character, see function scatterplot in the car package
# legend.left    - logical, if TRUE a legend in the top left corner is shown
# legend.right   - logical, if TRUE a legend in the top right corner is shown
# xlim           - NULL or numeric vector of length 2
# reset.par      - logical, see function scatterplot in the car package
# label.legend   - logical, if TRUE a legend at the top will be plotted depicting the
#                  marked labels
#
# value:
#      list containing the indices of up and downregulated ratios
#
# changelog: 20100219 implementation
#            20100223 some documentation
#            20100311 added transparancy to the points -> 'alpha'
#                     switched to function 'scatterplot' from the 'car' package
#            20100610 added parameter 'table', 'experiment', 'silac.state', 'norm', 'signif'
#                     the idea was to replace the old 'sigRatioPlot' function that was used by
#                     the 'summaryPDF' function to plot quantified protein groups
#            20100909 several different labels can be plotted in different colors
#            20100915 if 'norm=F'meaning the unnormalized ratios are used, no regulated
#                     protein groups are reported since the significance values are based
#                     based on normalized ratios
#            20101021 compatibility for MQ v1.1.x:
#                      - if there are no significance column in the table, all p-values are set
#                        to one
#            20110505 if no valid values are found a dummy plot will be created and the functions stops.
#            20110701 if no significance column can be found the significance B values are
#                     calculated using my implementation
#            20120705 - tickmarks at x-axis based on xlim values
#            20121121 -  'ylim' as parameter
#################################################################################################
sigRatioPlot.ph <- function(rat, int, p,  table=NULL, experiment="", silac.state=c("H.L", "M.L", "H.M"), norm=T, signif=c("B", "A"), label=NULL, pcut = 0.01, adjust=c("none", "BH"), label.string=c("Aurora", "Aurora-A"), col="black", sig.col="red", label.col="darkred", label.pch=17, label.all=F, alpha=100, pch=c(16, 16), boxplots=c("xy", "x", "y", "n"), legend.left=T, legend.right=T, xlim=NULL, ylim=NULL, reset.par=TRUE, label.legend=F, main="", ...){

    # original parameter
    opar <- par(no.readonly=TRUE)

    #############################################
    # if 'table' is specified, e.g. the protein
    # groups table
    #############################################
    if(!is.null(table)){

        ###################################
        # check whether this is a MQ table
        ###################################
        if(sum( grep("^Ratio.H.L.Normalized", colnames(table), ignore.case=T) ) == 0)
#        if( ( "Ratio.H.L.Normalized" %in% colnames(table)) == FALSE )
            stop("You did not supply a MaxQuant table!\n")

        # remove con/rev
        rm.idx <- union(grep("^\\+$", table[, grep("^Contaminant$", colnames(table), ignore.case=T )]), grep("^\\+$", table[, grep("^Reverse$", colnames(table), ignore.case=T)]) )
        if(length(rm.idx) > 0)
            table <- table[-rm.idx, ]

        # get the silac label
        silac.state <- match.arg(silac.state)

        # significance A or B
        signif <- match.arg(signif)


        # if there is an experiment defined...
        if(nchar(experiment) > 0){
                ratioColumn <- ifelse(norm, paste("Ratio.", silac.state,".Normalized.", experiment, sep=""),paste("Ratio.",silac.state,".", experiment, sep=""))
            intColumn <- paste("Intensity.", experiment, sep="")
            sigColumn <-  paste("Ratio.", silac.state, ".Significance.",signif,"..", experiment, sep="")

        } else {
            ratioColumn <- ifelse(norm, paste("Ratio.", silac.state,".Normalized", sep=""),paste("Ratio.",silac.state, sep=""))
            intColumn <- paste("Intensity", sep="")
            sigColumn <-  paste("Ratio.",silac.state,".Significance.",signif,".", sep="")
        }

        ########################################
        # get ratios and intensities
        ########################################
        rat <-  table[, grep( paste("^", ratioColumn, "$", sep=""), colnames(table), ignore.case=T)]
        int <-  table[, grep( paste("^", intColumn, "$", sep=""), colnames(table), ignore.case=T)]
        int[int == "-Inf"] <- NA

        # if 'sigColums' was not given correctly, all p values are set to one
        #if( (sigColumn %in% colnames(table)) ){
        if(sum(grepl( paste("^", sigColumn, "$", sep="") , colnames( table), ignore.case=T ) ) > 0){
            p <- table[, grep( paste("^", sigColumn, "$", sep="") , colnames( table), ignore.case=T )]
        } else {
            #p <- rep(1, length(rat))
            p <- significanceB( rat, int  )[[1]]

        }

       if(norm) main <- paste(main," (", experiment, " Normalized Ratio ",silac.state, ")", sep="")
       else main <- paste(main," (", experiment, " Ratio ",silac.state, ")", sep="")


    }


    # make it all numeric
    if(is.factor(rat))
        rat <- as.numeric(as.character(rat))
    if(is.factor(int))
        int <- as.numeric(as.character(int))
    if(is.factor(p))
        p <- as.numeric(as.character(p))


    require(car)

    ##################################
    # store the orignial index of
    # the data values before removing
    # missing values
    ##################################
    data.index <- 1:length(rat)
    names(data.index) <- 1:length(rat)

    ##################################
    # remove missing values
    ##################################
    rat.na <- which( is.na(rat) )
    #int.na <- which( is.na(int) | is.infinite(log(int, 10)) )
    int.na <- which( is.na(int))

    p.na <- which( is.na(p) )
    na.idx <- union(  union(rat.na, int.na), p.na )

   if(length(na.idx)>0){
       rat <- rat[-na.idx]
       int <- int[-na.idx]
       p <- p[-na.idx]
       data.index <- data.index[-na.idx]

       if(!is.null(label))
           label <- label[ -na.idx ]
   }

    ###################################
    # boxplots on axes
    ###################################
    boxplots <- match.arg(boxplots)

    ###################################
    # log
    ###################################
    rat.log <- log(rat, 2)
    int.log <- log(int, 10)
    int.log[is.infinite(int.log)] <- NA




    #########################################################################################
    #
    #          - if there are no valid values, stop right here.....
    #          - produce an empty plot for compatibility with PDF script
    #
    #########################################################################################
    if(length(intersect(which(!is.na(rat)),which(!is.na(int)))) < 1  ){

        plot(1,1, col='white', axes=F, xlab='', ylab='')

        return(1)
    }


    ###################################
    # number of quantified items
    ###################################
    quant.numb <- sum(!is.na(rat.log), na.rm=T)

    ###################################
    # plotting symbols
    ###################################
    pch.vec <- pch
    if(!is.null(label.string)){
        if(length(label.pch) != length(label.string) )
            label.pch <- rep(label.pch[1], length(label.string))
    }

    #########################################################
    #
    #  if unnormalized ratios are used it makes no sense
    #  to look for regulated proteins because the significance
    #  values are based on normalized ratios
    #
    #########################################################
    #if(!is.null(table) & norm){
    if(norm){
           ###################################
           #   adjust p-values
           ###################################
           adjust <- match.arg(adjust)
           p <- p.adjust(p, adjust)

           ###################################
           # determine regulated
           ###################################
           rat.log.reg <- rat.log[ p <= pcut]
           int.log.reg <- int.log[ p <= pcut]

           # index of regulated
           regulated <- ifelse(p <= pcut, 1, 0)

           # in numbers
           down.numb <- sum(rat.log.reg < 0, na.rm=T)
           up.numb <- sum(rat.log.reg > 0, na.rm=T)

           # index of regulated
           down.idx <- which((rat.log < 0) & (p <= pcut))
           up.idx <- which((rat.log > 0) & (p <= pcut))

           ###################################
           # add labels
           ###################################
           if(!is.null(label.string) & !is.null(label)){

              #################
              # multiple labels
              if(length(label.string) > 1 ){
                 label.count <- 1
                 #label.idx <- c()

                 # loop over labels
                 for(m in label.string){
                   label.idx.tmp <- grep(m, label, ignore.case=T)

                   #label.idx.tmp <- union( label.idx, grep(m, label, ignore.case=T))
                   if(length(label.idx.tmp) > 0){

                        # only significant features are labeled
                        if(!label.all)
                            label.idx.tmp <- intersect(which(regulated == 1), label.idx.tmp)
                        regulated[  label.idx.tmp ] <- label.count + 1
                   }
                   label.count = label.count + 1
                 }
                 # plotting character of labels
                 regulated <- as.factor(regulated)
                 pch.vec <- c(pch.vec, label.pch  )

             ##################
             # only a single label
             } else {
               label.idx <- grep(label.string, label, ignore.case=T)
               if(length(label.idx) > 0)
                   if(!label.all)
                       label.idx <- intersect(which(regulated == 1), label.idx)
                   regulated[ label.idx   ] <- 2

               # plotting character of labels
               regulated <- as.factor(regulated)
               pch.vec <- c(pch.vec, label.pch  )
            }

         } # end if label

    } else {
        regulated <- rep(0, length(p))
    }

#   if(length(na.idx)>0){
#       rat <- rat[-na.idx]
#       int <- int[-na.idx]
#       p <- p[-na.idx]
#       data.index <- data.index[-na.idx]#
#
#       if(!is.null(label))
#           label <- label[ -na.idx ]
#   }



    ###################################
    # colors
    ###################################
    col.org <- col
    sig.col.org <- sig.col

    col <- my.col2rgb(col, alpha)
    sig.col <- my.col2rgb(sig.col, alpha)


    if(length(label.col) < length(label.string) )
        label.col <- rep(label.col[1], length(label.string))
    for(cc in 1:length(label.col)){

        label.col[cc] <- my.col2rgb(label.col[cc], alpha)
    }

    # set the color palette: the first entry is just a dummy -> used in function 'scatterplot'
    #palette( c("black", col, sig.col, label.col[as.numeric(levels(regulated)[3:length(levels(regulated))]) - 1]))     # R version 2.10.1
    palette( c(col, sig.col, label.col[as.numeric(levels(regulated)[3:length(levels(regulated))]) - 1]))               # new R version

    ###################################
    # xlim  &  ylim
    ###################################

    if(is.null(xlim)){
        xmax <- max(abs(rat.log), na.rm=T)

        xlim <- c(-1*xmax, xmax)
    }

    # ylim
    if(is.null(ylim))
    {
        ylim <- c(min(int.log, na.rm=T),max(int.log, na.rm=T)+1)
    }

    ##########################################
    # original plotting parameters
    ##########################################
    #par.org <- par(no.readonly = TRUE)

    ##########################################
    # plot
    ##########################################
    scatterplot( rat.log ,int.log, xlim=xlim, ylim=ylim, pch=pch.vec, xlab=expression(log[2](Ratio)), ylab=expression(log[10](Intensity)), boxplots=boxplots, reg.line=F, smooth=F, groups=factor(regulated), legend.plot=F, main=main, axes=F, reset.par=F, ...)

    # axes
    #axis(1, at=seq( ceiling(-max(abs(rat.log))), floor(max(abs(rat.log))), 1 )  )
    axis(1, at=seq( ceiling(xlim[1]), floor(xlim[2]), 1 )  )
    axis(2)
    #grid(ny=seq( ceiling(-max(abs(rat.log))), floor(max(abs(rat.log))), 1 ))


    # legends
    if(legend.left){
        if(norm) legend("topleft", legend=c(paste("down:", down.numb), paste("up:", up.numb), paste("quantified:", quant.numb) ), bty="n", inset=c(0.02,0))
        else legend("topleft", legend=c(paste("quantified:", quant.numb) ), bty="n", inset=c(0.02,0))
    }
    if(legend.right & norm )
        legend("topright", legend=c( paste("p adjustment:", adjust), paste("p cutoff:", pcut)), bty="n"  )

    if(label.legend)
        legend("top", legend=label.string, col=label.col, pch=label.pch, bty="n"  )


    out <- vector("list", 2)
    names(out) <- c("down", "up")

    if(norm){
        out[[1]] <- data.index[down.idx]
        out[[2]] <- data.index[up.idx]
    }

    # default color palette
    palette("default")

    # default plotting parameters
    if(reset.par)
        par(opar)


    return(out)
}


###########################################################################################
#
#                          incorporation check
#
# file    - character specifying a file name OR a peptide/evidence table as
#           R-object
#
# changelog:  20110405 parameter 'xlim'
###########################################################################################
checkIncorporation <- function(file, xlim=c(0.2, 1))
{


    #####################################
    #         data preparation
    #####################################

    if(mode(file) == "character")
    {
     	# import the data
   	pep <- read.delim(file, row.names=1, stringsAsFactors=FALSE,  blank.lines.skip = TRUE)

    } else if( is.null(dim(file))) {

        plot(0,0,axes=F, main="")
        return(1)
	#stop("This is not a peptides/evidence file!\n")

    } else
	pep <- file

    # remove reverse hits and contaminants
    #pep <- rmConRev(pep)

    revIdx <- grep("^\\+$", pep[, "Reverse"])
    if( length(revIdx) > 0 )
        pep <- pep[-revIdx, ]
    conIdx <- grep("^\\+$", pep[, "Contaminant"])
    if( length(conIdx) > 0 )
        pep <- pep[-conIdx, ]

    # store the current table, i.e. the table containing all peptides except reverse hits and contaminants
    pepOrg <- pep

    # remove peptides that don't have a ratio
    pep <- pep[!is.na(pep[, "Ratio.H.L"]),  ]


    ################################################
    #
    ################################################


    ################################################
    #     check if there are any ratios left....
    #################################################
    if(dim(pep)[1] < 2){

	# if there are also no intensities just stop...
	if( is.na(  sum(pepOrg[, "Intensity"] > 0)  )  )
		stop("There are no ratios and no intensities written in the file you\'ve uploaded!!\nIf you are using the \'evidence.txt\' give the \'peptides.txt\' a try...\n")

	barplotInt(pepOrg)

    } else{

	####################################################
	#  ... otherwise one can just proceed as usual
	#####################################################
	# estimate density to plot
	A.dens <- density( 1 - (1/(pep[, "Ratio.H.L"]+1)) )

	# total incorporation
	incTot <- median(pep[, "Ratio.H.L"]) / ( 1 + median(pep[, "Ratio.H.L"]) )


	#######################################
	# check if one or two labels were used
	#######################################

	# two labels
	if( length(  grep("[R|K]\\.Count" ,colnames(pep))  ) == 2  ){

		# get the R/K peptides
		Rpep <- pep[ pep[, "R.Count"] > 0 & pep[, "K.Count"]==0 ,   ]
		Kpep <- pep[ pep[, "K.Count"] > 0 & pep[, "R.Count"]==0 ,   ]


		if(dim(Kpep)[1] > 1){
			K.dens <- density( 1 - (1/(Kpep[, "Ratio.H.L"]+1)) )
			# incorporation
			incK <- median(Kpep[, "Ratio.H.L"]) / ( 1 + median(Kpep[, "Ratio.H.L"]) )

			# for the legend
			KlegendText <-  paste("Lysine peptides (", dim(Kpep)[1] ,"): ", round( incK, 3 )  , sep="" )
			Kdens.max = max(K.dens$y)



		} else {
			Kdens.max = 0
 			KlegendText = paste("Lysine peptides (", dim(Kpep)[1] ,")")
		}

		if(dim(Rpep)[1] > 1){
			R.dens <- density( 1 - (1/(Rpep[, "Ratio.H.L"]+1)) )
			# incorporation
			incR <- median(Rpep[, "Ratio.H.L"]) / ( 1 + median(Rpep[, "Ratio.H.L"]) )

			# for the legend
			RlegendText = paste("Arginine peptides (", dim(Rpep)[1] ,"): ", round( incR, 3 )  , sep="" )
			Rdens.max <- max(R.dens$y)
		} else{
			Rdens.max = 0
 			RlegendText = paste("Arginine peptides (", dim(Rpep)[1] ,")")

		}


		## PLOT
		par(mfrow=c(1,2))

		# RATIOS
		plot( A.dens,  main="Incorporation Rate", xlab="1-(1/(ratio+1 ))", ylim=c(0, max(A.dens$y, Kdens.max, Rdens.max)  ), xlim=xlim )
		if(Rdens.max > 0)
			lines( R.dens, col="red" )
		if(Kdens.max > 0)
			lines( K.dens, col="green" )

		legend("topleft", legend=c( paste("all peptides (", dim(pep)[1] ,"): ", round( incTot,3 )  , sep=""), RlegendText, KlegendText  ), fill=c("black", "red", "green") )

		# INTENSITIES
		barplotInt(pepOrg)
		par(mfrow=c(1,1))



	} else if( length(  grep("[R|K].Count" ,colnames(pep))  ) == 1  ){

		par(mfrow=c(1,2))
		# RATIO
		plot( A.dens,  main="Incorporation Rate", xlab="1-(1/(ratio+1 ))" )
		legend("topleft", legend=c( paste("all peptides (", dim(pep)[1] ,"): ", round( median(pep[, "Ratio.H.L"])/(1+ median(pep[, "Ratio.H.L"])),3 )  , sep="") ), fill="black" )

		# INTENSITIES
		barplotInt(pepOrg)

		par(mfrow=c(1,1))


	} else stop("No SILAC label found!")

    }


    #######################################
    #            some output
    #######################################
   output <- list()
   output[[1]] <- pepOrg
   names(output) <- c("table")

   return(output)

}
#######################################################
# help function for 'checkIncoproration'
#    barplot: number of intensities != NA
#
#######################################################
barplotInt <- function(pepOrg)
{
	# intensities
    	intAll <- sum( pepOrg[, "Intensity"] > 0, na.rm=T)
    	intL <- sum(pepOrg[, "Intensity.L"] > 0, na.rm=T)
    	intH <- sum(pepOrg[, "Intensity.H"] > 0, na.rm=T)


	x <- c(intAll, intL, intH)
	names(x) <- c("Intensity", "Intensity L", "Intensity H")

	# check if there is an 'NA' entry
	x[which(is.na(x))] <- 0
	ymax <- max( x)
	ymin=0

	# ceck if all entries are zero
	if(sum(x)==0) ymax=1

	barplot(x, main="Number of peptides/evidences having intensities above zero", ylim=c(ymin, ymax + ymax*.1) )
	for(i in 1:length(x))
		text(i, x[i]+ymax*.05, x[i], cex=1.5)

}


############################################################################################
#                               q.value
#
# - compute q-values for each item in a MaxQuant table and add
#   a new column to the table
#
#   ways to calculate significance (Käll et al., JPR 2008)  (not complete...)
#
#   simple FDR:       - for a given score threshold, count the number of decoy hits above the threshold
#                       and the number of target hits above the threshold
#                     - estimate FDR by computing the ratio of these two values
#   q-values:         - minimal FDR
#
#
#
# changelog:    20100820 implementation
#               20110322 changed the overall number N to consider only forward hits
#               20110414 changed calculation of q-values
#                         - # reverse hits above score of current reverse hit divided by
#                           # number of forward hits above score of current reverse hit
#               20110414 bugfix: q-values for target hits coming after the last reverse are now computed...
############################################################################################
q.value <- function(tab, rev.col="Reverse", base.on="PEP", decreasing=F ){

    #base.on <- match.arg(base.on)

    # number of table entries
    N <- dim(tab)[1]                                                    # forward and reverse
    N.rev <- length(grep("^\\+$", tab[, rev.col]))                      # reverse
    N.fwd <- N - N.rev                                                  # forward



    # original ordering
    ord.org <- rownames(tab)

    # order accoring to some score
    tab <- tab[order( tab[, base.on], decreasing=decreasing),  ]

    # initialize a vector for storing the q-values
    q.val <- vector("numeric", N)
    names(q.val) <- rownames(tab)

    # get the index of reverse hits
    rev.idx <- grep("\\+", tab[, rev.col])

    # if there are no reverse hits report a warning
    if(length(rev.idx) == 0){

        warning("\nNo reverse hits found!!\n")
        q.val <- rep(0, N)

    } else{

        rev.count <- 1

        # loop over the reverse hits
        for(r in rev.idx){

            # if its the first reverse hit
            if(r == rev.idx[1]){

               q.val[1:(r-1)] = 0

            } else {

                # number of reverse hits / number of forward hits above score threshold of current reverse hit
                q.val[(rev.idx[rev.count-1]+1):(rev.idx[rev.count]-1)] <- ((rev.count-1)/(rev.idx[rev.count]-rev.count)) # * (N.fwd/N.rev)

                #q.val[(rev.idx[rev.count-1]+1):(rev.idx[rev.count]-1)] <- rev.count/(N-rev.count)

                # if its the last reverse hit, add q-values for target hits coming after the last reverse hit
                if(r == rev.idx[length(rev.idx)] ){
                    if(r < N)
                          q.val[(rev.idx[rev.count] + 1): length(q.val)] <- rev.count/N.fwd
                }
            }

            rev.count <- rev.count +1

        }

        # set the q-values of reverse hits to 1
        q.val[rev.idx] <- 1

    }

    # append the q-values to the table
    tab <- cbind(tab, q.val)

    # restore the original ordering
    tab <- tab[ord.org, ]

    return(tab)

}


#############################################################################################
#                                    getProteinState
#
#  given the proteinGroups table of an SILAC experiment processed with an 'experimentalDesign' file
#  this function determines for each experiment and for each protein group id whether it was detected and/or
#  quantified in the particular experiment
#
#  arguments:
#    proteinGroups    - dataframe of characters, the proteinGroups table
#    plot             - logical
#    minPep           - minimal number of peptides such that a peptide is DETECTED in
#                       the particular experiment
#
#
# changelog:   20091019 implementation
#              20091028 paramater 'minPep'
#              20101118 added 'pep.col' to be more flexible; version 1.1.1.26 the column names
#                       changed
#############################################################################################
getProteinState <- function( proteinGroups, evidence=NULL, plot=T, minPep=1, pep.col="Peptides..seq..", ...  ){

    ####################################################
    #           get the experiments
    # the columns 'Experiment.X' indicates the number of
    # peptides used for QUANTITATION from experiment X
    ####################################################
    ex.idx <- grep( "^Experiment" , colnames(proteinGroups) )
    if( length(ex.idx) == 0  ) stop("\n\nNo experiments defined. Did you specify an 'experimentalDesign.txt' when running 'Identify.exe'??\n\n\n")

    ex <- colnames(proteinGroups)[ ex.idx  ]
    ex <- sub("^Experiment\\.", "", ex)

    ####################################################
    #   determine the type of the SILAC experiment,
    #   i.e.
    #      doublets: silac.type = 0
    #      triplets: silac.type = 1
    ####################################################
    silac.type <- ifelse(length(  grep( "^Ratio.H.M.", colnames(proteinGroups) ) > 0 ), 1, 0 )


    ####################################################
    #
    #  now loop over all protein groups and determine
    #  whether the protein was detected and/or quantified
    #
    #  detected:   column 'Peptides..seq..X' > minPep
    #  quantified: column 'Ratio.H.L.X' != NA
    #
    ####################################################
    detected <- quantifiedHL  <- vector("list", length(ex) )
    names(detected) <- names(quantifiedHL) <- ex

    if(silac.type == 1){
         quantifiedML <- quantifiedHM  <- vector("list", length(ex) )
         names(quantifiedML) <- names(quantifiedHM) <- ex
    }

    if( !is.null(evidence) ){
         detectedButNotQuantified <- vector("list", length(ex))
         names(detectedButNotQuantified) <- ex
    }

    for(e in ex){

        detected[[e]] <- ifelse( proteinGroups[, paste(pep.col, e , sep="") ] >= minPep , 1, 0  )
        quantifiedHL[[e]] <- ifelse( is.na(proteinGroups[, paste("Ratio.H.L.",e , sep="") ]) , 0, 1  )


        # protein group ids as vector names
        names( detected[[e]] ) <- names( quantifiedHL[[e]]  ) <- rownames(proteinGroups)

        ########################################
        #  check for non-quantified protein the
        #  SILAC state of all evidences
        #
        ########################################
        if(!is.null(evidence)){

            # protein groups that were detected but not quantified
            pgIds <- rownames(proteinGroups)[ which( (quantifiedHL[[e]] == 0) &  (detected[[e]] == 1) ) ]


            # list to store the evidences for each protein group
            evidencesOfProteinGroups <- vector("list", length(pgIds))
            names(evidencesOfProteinGroups) <- pgIds


            # loop over those protein groups
            for(p in pgIds){
                # get all associated evidences
                evidencesOfProteinGroups[[p]] <- evidence[intersect(grep("^(.*;)*0(;.*)*$", evidence[, "Protein.Group.IDs"] ) ,grep(paste("^",e,"$", sep=""),evidence[, "Experiment"])), c("Protein.Group.IDs", "Experiment", "Sequence", "SILAC.State", "Type")]
            }

            detectedButNotQuantified[[e]] <- evidencesOfProteinGroups
        }

        # SILAC triplets
        if(silac.type == 1){

            quantifiedML[[e]] <- ifelse( is.na( proteinGroups[, paste("Ratio.M.L.",e , sep="") ]) , 0, 1  )

            columns <- c("Protein.Descriptions", "Gene.Names")

            #specific.prot.mat <- matrix("", ncol=length(columns) + 1, nrow=0)
            #colnames(specific.prot.mat) <- c( "Tissue", columns )

            #for(ti in tissue){

             #   tmp <- pg[names(which(protDetectMatSpecific[, ti]==1)),  ]
             #   tmp <- cbind(rep(ti, dim(tmp)[1]), tmp  )

              #  specific.prot.mat <- rbind(specific.prot.mat, tmp)

            #}
            quantifiedHM[[e]] <- ifelse( is.na(proteinGroups[, paste("Ratio.H.M.",e , sep="") ]) , 0, 1  )

            names(quantifiedML[[e]]) <- names(quantifiedHM[[e]]) <- rownames(proteinGroups)
        }
    }


    ###############################################################
    # make a barplot: detected and quantified proteins per tissue
    ###############################################################
    if(plot){

        detectedBar <- unlist(lapply( detected, sum, na.rm=T  ))

        if(silac.type == 0){

            quantifiedBar <- unlist( lapply( quantifiedHL, sum, na.rm=T ) )
            par(mar=c(8,4,4,2))
            barplot( rbind(detectedBar, quantifiedBar  ), beside=T, main="Detected and quantified protein groups per experiment", col=c("darkgrey", "lightgrey"), ylim=c(0, max(detectedBar)+(max( detectedBar ))*0.15    ), axes=T, las=2,...  )

            legend("topright", legend=c(paste("detected (at least",minPep,"peptide(s))"), "quantified"), fill=c("darkgrey", "lightgrey"))
        } else {

            quantifiedBarHL <- unlist( lapply( quantifiedHL, sum ) )
            quantifiedBarML <- unlist( lapply( quantifiedML, sum ) )
            quantifiedBarHM <- unlist( lapply( quantifiedHM, sum ) )

            par(mar=c(8,4,4,2))
            barplot( rbind(detectedBar, quantifiedBarHL, quantifiedBarHM, quantifiedBarML  ), beside=T, main="Detected and quantified protein groups per experiment", col=c("black", "darkgrey", "grey" ,"lightgrey"), ylim=c(0, max(detectedBar)+(max( detectedBar ))*0.3    ), ... )

            legend("topright", legend=c(paste("detected (at least",minPep,"peptide(s))"), "quantified H/L", "quantified H/M", "quantified M/L"), fill=c("black", "darkgrey", "grey", "lightgrey"))


        }
    }


    #########################
    #    output
    #########################
    if(silac.type == 0){
       output <- vector("list", 4)
       names(output) <- c("experiments", "detected", "quantified")

       output[[1]] <- ex
       output[[2]] <- detected
       output[[3]] <- quantifiedHL

       if(!is.null(evidence)){
             output[[4]] <- detectedButNotQuantified
             names(output)[[4]] <- "detectedButNotQuantified"
         }

    } else {
        output <- vector("list", 5)
       names(output) <- c("experiments", "detected", "quantifiedHL", "quantifiedHM", "quantifiedML")

       output[[1]] <- ex
       output[[2]] <- detected
       output[[3]] <- quantifiedHL
       output[[4]] <- quantifiedHM
       output[[5]] <- quantifiedML


    }

    return(output)
}


################################################################################################
#                                significance A
#
#
#
# changelog:   20100420 implementation
#              20110531 finally changed to the formula without 0.5 ...
################################################################################################
significanceA <- function( r, log =T, log.base=2 ){

    ###################################################
    #  tranform to log scale
    ###################################################
    if(log)
        r <- log(r, log.base)

    quant.idx <- which(!is.na(r))
    ##################################################
    # robust and asymetrical estimate of the sd
    #
    ##################################################
    estimate.sd <- quantile(r[quant.idx], c(0.1587, 0.5, 0.8413))

    # the notation used in the MQ supplement
    r.minus1 <- estimate.sd[1]
    r.zero <- estimate.sd[2]
    r.plus1 <- estimate.sd[3]

    ##################################################
    # z tranformation
    ##################################################
    r.z <- ifelse(r >= r.zero, (r-r.zero)/(r.plus1-r.zero),  (r.zero-r)/(r.zero-r.minus1) )


    ##################################################
    # complementary error function
    ##################################################
    #sigA <- 0.5*erfc( r.z/sqrt(2) )
    sigA <- erfc( r.z/sqrt(2) )

return(sigA)
}


#############################################################
#                 Significance B
# r        - ratios
# i        - corresponding intensities
# nbin     - numeric, number of proteins in an intensity bin
#          - each bin contains at least that number of proteins
#
# changelog:  20100423 implementation
#             20110616 if length(r) < nbin sigA is applied
#############################################################
significanceB <- function(  r, i, nbin = 300, log=T, log.base=2  ){

    names(r) <- names(i) <- 0:(length(r)-1)

    # original order
    org.order <- names(r)

    ###################################################
    #  tranform to log scale
    ###################################################
    if(log)
        r <- log(r, log.base)

    ###################################################
    # order according to intensities
    ###################################################
    order.idx <- order(i, decreasing=T)
    i <- i[ order.idx ]
    r <- r[ order.idx ]

    ###################################################
    # index of all proteins having a ratio
    ###################################################
    quant.idx <- which(!is.na(r))

    ###################################################
    # get the ratios
    ###################################################
    r.quant <- r[quant.idx]
    i.quant <- i[quant.idx]

   ###################################################
    # bin the intensities
    ###################################################
    #bin <- round(seq(0, length(i[quant.idx]), length.out=ceiling(length(i[quant.idx])/nbin)  ))
    bin <- floor(seq(0, length(i.quant), length.out=ceiling(length(i.quant)/nbin)  ))
    #bin <- round(seq(0, length(i[quant.idx]), by=nbin  ))
    #bin <- round(seq(0, length(i), length.out=ceiling(length(i)/nbin)  ))


    ###################################################
    # loop over the bins and apply significance A
    ###################################################
    # vector to store the values
    sigB <- rep(NA, length(r))
    names(sigB) <- names(r)

    ###################################################
    #
    ###################################################
    if(length(r.quant)  <= nbin){

      sigB[ names(r.quant) ] <-  significanceA(r.quant, log=F)

      sigB <- sigB[ org.order ]

      return(list(sigB))
    }



    # store some informations on the bins
    info <-vector("list", length(bin)-1)
    names(info) <- paste("bin", 1:(length(bin)-1))

    # loop over the bins
    for(b in 1:(length(bin)-1)){

        # get the ratios in the current bin
        r.quant.bin <-  r.quant[ (bin[b]+1) : bin[b+1] ]

        # apply sigA
        sigB.tmp <- significanceA ( r.quant.bin, log=F)
        sigB[names(r.quant.bin)] <- sigB.tmp

        # some additional information
        info.tmp <- vector("list", 4)
        names(info.tmp ) <- c("proteins", "bin", "ratio", "intensity")
        info.tmp[[1]] <- length(r.quant.bin)
        info.tmp[[2]] <- names(r.quant.bin)
        info.tmp[[3]] <- r[names(r.quant.bin)]
        info.tmp[[4]] <- i[names(r.quant.bin)]

        info[[b]] <- info.tmp
    }

    ####################################################
    #   restore the original ordering
    ####################################################
    sigB <- sigB[org.order]

    out <- list()
    out[[1]] <- sigB
    out[[2]] <- info

    return(out)

}



#########################################################
#
#        complementary error function assuming
#        a normal distribution
#
#########################################################
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)


#################################################################################################
#                               check site redundancy
# - for all unlocalized sites, check whether there are more sites reported than the number
#   of sites seen on the peptide
# - remove all redundant sites, i.e. report only those with best localization score
#
# value:
#   - list, for each peptide id the site ids are returned
#
# changelog:  20120625 - version 1.2.7.5:
#                      - 'Number of Phospho STY' contains now the number of modifications
#                        for all peptide versions carrying the site -> made compatible again
#             20120706 - remove site swith prob 0
#
#################################################################################################
getNRSites <- function(tab, mod="Phospho..STY."){

    # check whether there are sites with loc prob of zero
    if(min( tab[,grep("^Localization.Prob$", colnames(tab), ignore.case=T) ], na.rm=T ) == 0 ){

        zero.idx <- which(tab[,grep("^Localization.Prob$", colnames(tab), ignore.case=T) ] == 0)

        tab <- tab[-zero.idx, ]

    }

    # all non-redundant peptide ids
    pep.unique <- unique(tab[, grep( "^Peptide.IDs$", colnames(tab), ignore.case=T ) ])

    # list to store the ids of filtered sites
    nr.sites <- vector("list", length(pep.unique))
    names(nr.sites) <- pep.unique

    # loop over the non-redundant peptide ids
    for(pep in pep.unique){

        # get all phosphosites associated with the current peptide id
        tab.pep <- tab[ which(tab[, grep( "^Peptide.IDs$", colnames(tab), ignore.case=T ) ] == pep),   ]

        # sort this table according to localization probability
        tab.pep <- tab.pep[order(tab.pep[, grep("^Localization.Prob$", colnames(tab), ignore.case=T)], decreasing=T),  ]

        # maximal number of modifications that peptide was seen, i.e. # mods the peptide was measured
        # problem: a site can be measured on different peptides (missed cleavages) with different
        # numbers of modifications -> take the highest number
        nMod.tmp <- tab.pep[, grep( paste("^Number.of.", mod, sep=""), colnames(tab.pep), ignore.case=T) ]

        if( sum(nchar(nMod.tmp), na.rm=T) == 0 ){
            nMod=dim(tab.pep)[1]
        } else {
             nMod <- max(unlist( lapply( strsplit( as.character(nMod.tmp[nchar(nMod.tmp) > 0 ]), ";"), function(x) max( as.numeric( x)) )) )
        }

        # if there are more sites reported than the peptide was seen with, choose the sites with highest localization prob
        if(dim(tab.pep)[1] > nMod){

            # get the 'nMod' best localized sites
            tab.pep.filt <- tab.pep[1:nMod, ]

        } else{
            tab.pep.filt <- tab.pep
        }

        # store the site ids in the output list
        nr.sites[[pep]] <- rownames(tab.pep.filt)
    }

    return(nr.sites)
}


##################################################################################################
#          check the overlap of modified and quantified sites between different experiments
#
# changelog: 20100214 implementation
#            20120625 compatibility to 1.2.7.5
##################################################################################################
siteOverlap <- function(siteTable, experiments, phPEP = 0.01, phLoc =.75, quant=T, quantCol.prefix="Ratio.H.L.Normalized"){

    #if(length(experiments %in% colnames(siteTable)) < length(experimens)  ) stop("Check your experiments!\n")

    ##########################################
    # a site is identified in experiment X if
    # there is a PEP value below 'phPEP' for
    # experiment X
    #########################################
    pepExp.idx <- locExp.idx <- vector("character", length(experiments))
    for(ex.idx in 1:length(experiments)){
        pepExp.idx[ex.idx] <- colnames(siteTable)[grep(paste("^PEP.",experiments[ex.idx], "$", sep=""), colnames(siteTable), ignore.case=T)]
        locExp.idx[ex.idx] <- colnames(siteTable)[grep( paste("^Localization.Prob.",experiments[ex.idx], "$", sep=""), colnames(siteTable), ignore.case=T)]
    }

    pepExp <- siteTable[, pepExp.idx ]
    locExp <- siteTable[, locExp.idx ]

    if(quant) quantExp <- siteTable[, paste(quantCol.prefix, experiments, sep=".")]


    olPEP.all <- vector("list", length(experiments))
    names(olPEP.all) <- experiments

    olPEP.filt <- olLocPEP.filt <- olLoc.filt <- olQuant <- olQuantLoc.filt <- olQuantLocPEP.filt <- olPEP.all

    for(ex in experiments){

        olPEP.all[[ex]] <- rownames(siteTable)[which( !is.na(pepExp[, grep( paste("^PEP.",ex, "$", sep=""), colnames(pepExp), ignore.case=T)]) )]
        olPEP.filt[[ex]] <- rownames(siteTable)[ which(pepExp[, grep( paste("^PEP.",ex, "$",sep=""),  colnames(pepExp), ignore.case=T) ] <= phPEP)  ]

        olLoc.filt[[ex]] <- rownames(siteTable)[ which(locExp[,  grep( paste( "^Localization.Prob.",ex, "$",sep="") , colnames(locExp), ignore.case=T)] >= phLoc)  ]
        olLocPEP.filt[[ex]] <- rownames(siteTable)[ (locExp[,  grep( paste("^Localization.Prob.",ex, "$", sep=""), colnames(locExp), ignore.case=T )] >= phLoc) & (pepExp[, grep( paste("^PEP.", ex, "$", sep=""), colnames(pepExp), ignore.case=T)] <= phPEP)  ]

        if(quant) {
            olQuant[[ex]] <- rownames(siteTable)[ which(!is.na( quantExp[, paste(quantCol.prefix, ex, sep=".") ]))  ]
            olQuantLoc.filt[[ex]] <- rownames(siteTable)[ ( !is.na( quantExp[, paste(quantCol.prefix, ex, sep=".") ])) &  (locExp[, grep( paste("^Localization.Prob.",ex, "$",sep=""), colnames(locExp), ignore.case=T )] >= phLoc)  ]
            olQuantLocPEP.filt[[ex]] <- rownames(siteTable)[(!is.na( quantExp[, paste(quantCol.prefix, ex, sep=".") ]) & (locExp[, grep( paste("^Localization.Prob.",ex, "$", sep=""), colnames(locExp), ignore.case=T )] >= phLoc)  ) & (pepExp[, grep(paste("^PEP.",ex, "$",sep=""), colnames(pepExp), ignore.case=T) ] <= phPEP) ]
        }
    }
    ##############
    # output
    ##############
    out <- vector("list",7)
    names(out) <- c("all sites", "PEP sites", "localized sites", "localized PEP sites", "all quantified sites", "quantified localized sites", "quantified localized PEP sites")    # names are used in summary-function!!

    out[[1]] <- olPEP.all
    out[[2]] <- olPEP.filt
    out[[3]] <- olLoc.filt
    out[[4]] <- olLocPEP.filt
    out[[5]] <- olQuant
    out[[6]] <- olQuantLoc.filt
    out[[7]] <- olQuantLocPEP.filt
    return(out)
}



##################################################################################################
#
#                                  compareExperiments
#
# The function determines all identified peptides/protein groups identified in each experiment
# defined in the 'experimentalDesign.txt'. Furthermore, it determines the number of peptides
# that were only found in 1, 2, ..., N experiments as well as  the overlap between the experiments
#
# arguments
#  evidence            - dataframe containing the 'evidence.txt' file
#  file                - NULL or a string specifying the filename of the produced pdf file
#  experiments         - NULL or character vector specifying a subset of the experiments
#  splitPDF            - logical, if TRUE each plot will result in a single pdf file
#                        defined in 'experimentalDesign'
# value
#                      a list containing the following elements
#                       'experiment'
#                            - list containing for each experiment the corresponding
#                              'evidence' rows as well as the numbers of identified
#                              peptides and proteins
#                       '# peptides'
#                            - vector containing for each number of experiments the
#                              number of identified peptide, e.g. '2' contains the number
#                              of peptides that were identified in two experiments
#                       '# proteins'
#                            - same as '# peptides' but for identified protein groups
#                       'overlap peptides'
#                            - matrix containing for each combinantion of experiments
#                              the overlap of identified peptides
#                       'overlap proteins'
#                            - same as 'overlap proteins' but for identified protein groups
#                              based on UNIQUE evidences
#
# changed: 20090319    - implementation
#          20090321    - parameter 'experiments'
#          20090322    - parameter 'splitPDF'
#          20090323    - solved a bug that causes that the protein group overlap between
#                        the experiments was not determined correctly ( I forgot to split
#                        the protein group ids by ";" in cases evidences were assigned to
#                        multiple protein groups)
#          20090507    - replaced two for-loops by C functions
#          20090512    - updated the function such that no experimentalDesign file is needed anymore
#
#          20101020    - when counting protein groups only unique evidences are used now!
#          20120626    - compatible with 1.2.7.5
##################################################################################################
compareExperiments <- function(evidence,  file="overlapExperiments.pdf", experiments=c("GeLCsample", "GeLCpellet"), splitPDF=T)
{


    ###########################################
    # check whether  Experiments were defined
    ###########################################
    if( length(grep("^Experiment$", colnames(evidence))) == 0  )
        stop("\nNo experiments defined!\n")


    #########################################
    # determine the different experiments
    #########################################
    expAll <- sort(unique(as.character(evidence[, grep("^Experiment$", colnames(evidence), ignore.case=T) ])))
    if(!is.null(experiments))
    {
        # check if the experiments were defined in 'evidence.txt'
        if( sum( expAll %in% experiments  ) != length(experiments)  )
            stop("\nExperiments you defined are not contained in 'evidence.txt'!\n")

        # extract the corresponding evidences
        exCount = 1
        for(ex in experiments)
        {
             if(exCount > 1) tmp <- rbind(tmp, evidence[which( as.character(evidence[, grep("^Experiment$", colnames(evidence), ignore.case=T) ]) == ex )  ,])
             else tmp <- evidence[which( as.character(evidence[, grep("^Experiment$", colnames(evidence), ignore.case=T) ]) == ex )  ,]

             exCount = exCount+1
        }

        evidence <- tmp
    }
    else
        experiments <- expAll

    #########################################
    # check whether SILAC
    #########################################
    if(length(grep("^SILAC", colnames(evidence))) > 0 ){

        if(length(grep("Medium", names(table(evidence[ ,'SILAC.State'])))) > 0){
            ##########################
            # triple SILAC
            ##########################
            silac.trip = T
            silac.doub = F

            # get L/H/M evidences: SILAC state is either written as heavy/light(-> MSMS level) or there is a heavy intensity (i.e. this
            # form of the peptide was not sequenced)
            evidence.H.idx <- union( grep("Heavy", evidence$SILAC.State), which(evidence$Intensity.H > 0) )
            evidence.L.idx <- union( grep("Light", evidence$SILAC.State),  which(evidence$Intensity.L > 0) )
            evidence.M.idx <- union( grep("Medium", evidence$SILAC.State), which(evidence$Intensity.M > 0))

        } else {
            ##########################
            # double SILAC
            ##########################
            silac.doub = T
            silac.trip = F

            # get L/H evidences
            evidence.H.idx <- union( grep("Heavy", evidence$SILAC.State), which(evidence$Intensity.H > 0))
            evidence.L.idx <- union( grep("Light", evidence$SILAC.State), which(evidence$Intensity.L > 0))

        }

    } else{
        silac.trip = silac.doub = F
    }

    #########################################
    # table containing only unique evidences
    #########################################
    not.unique.idx <- grep(";", evidence[ , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T) ])

    if(length(not.unique.idx) > 0 ){
        evidence.unique <- evidence[-not.unique.idx,  ]

        if(silac.doub){
            evidence.unique.L.idx <- setdiff(evidence.L.idx, not.unique.idx)
            evidence.unique.H.idx <- setdiff(evidence.H.idx, not.unique.idx)
        }
        if(silac.trip){
            evidence.unique.L.idx <- setdiff(evidence.L.idx, not.unique.idx)
            evidence.unique.M.idx <- setdiff(evidence.M.idx, not.unique.idx)
            evidence.unique.H.idx <- setdiff(evidence.H.idx, not.unique.idx)
        }

    } else {
        evidence.unique <- evidence
        if(silac.doub){
            evidence.unique.L.idx <- evidence.L.idx
            evidence.unique.H.idx <- evidence.H.idx
        }
        if(silac.trip){
            evidence.unique.L.idx <- evidence.L.idx
            evidence.unique.M.idx <- evidence.M.idx
            evidence.unique.H.idx <- evidence.H.idx
        }
    }

    #########################################
    # get the peptides per experiment
    #########################################
    pepPerExp <- vector("list", length(experiments))
    names(pepPerExp) <- experiments

    # output lists for double SILAC
    if(silac.doub){
        pepPerExpL <- pepPerExpH <- vector("list", length(experiments))
        names(pepPerExpL) <- names(pepPerExpH) <- experiments
    }
    # output lists for triple SILAC
    if(silac.trip){
        pepPerExpL <- pepPerExpM <- pepPerExpH <- vector("list", length(experiments))
        names(pepPerExpL) <- names(pepPerExpM) <- names(pepPerExpH) <- experiments
    }

    ##########################################################
    #  loop over the experiments
    ##########################################################
    for( ex in experiments)
    {
        #####################
        # double SILAC
        #####################
        if(silac.doub){
           exOut <- vector("list", 13)
           names(exOut) <- c("evidence", "evidence L", "evidence H", "unique evidence",  "unique evidence L", "unique evidence H","# peptides", "# nr peptides", "# nr peptides L", "# nr peptides H", "# protein groups", "# protein groups L", "# protein groups H")

        #######################
        # triple SILAC
        #######################
        } else if(silac.trip) {

           exOut <- vector("list", 17)
           names(exOut) <- c("evidence", "evidence L", "evidence M", "evidence H", "unique evidence",  "unique evidence L", "unique evidence M", "unique evidence H", "# peptides", "# nr peptides", "# nr peptides L", "# nr peptides M", "# nr peptides H", "# protein groups", "# protein groups L", "# protein groups M", "# protein groups H")

        #######################
        # no SILAC
        #######################
        } else {
           exOut <- vector("list", 5)
           names(exOut) <- c("evidence", "unique evidence", "# peptides", "# nr peptides", "# protein groups")
        }
        #######################################
        # get the evidences for each experiment
        #######################################
        evidence.ex.idx <- grep(paste("^", ex, "$", sep=""), evidence[, grep("^Experiment$", colnames(evidence), ignore.case=T)])
        if(silac.doub){
            evidence.ex.idx.L <- intersect(evidence.ex.idx, evidence.L.idx)
            evidence.ex.idx.H <- intersect(evidence.ex.idx, evidence.H.idx)
        }
        if(silac.trip){
            evidence.ex.idx.L <- intersect(evidence.ex.idx, evidence.L.idx)
            evidence.ex.idx.M <- intersect(evidence.ex.idx, evidence.M.idx)
            evidence.ex.idx.H <- intersect(evidence.ex.idx, evidence.H.idx)
        }

        #exOut[["evidence"]] <- evidence[ which(evidence[, "Experiment"] == ex),  ]
        exOut[["evidence"]] <- evidence[ evidence.ex.idx,  ]

        if(silac.doub){
            exOut[["evidence L"]] <- evidence.ex.idx.L
            exOut[["evidence H"]] <- evidence.ex.idx.H
        }
        if(silac.trip){
            exOut[["evidence L"]] <- evidence.ex.idx.L
            exOut[["evidence M"]] <- evidence.ex.idx.M
            exOut[["evidence H"]] <- evidence.ex.idx.H
        }

        ######################################
        # get unique evidences for each exp
        ######################################
        exOut[["unique evidence"]] <- evidence.unique[ which(evidence.unique[, grep("^Experiment$", colnames(evidence), ignore.case=T)] == ex),  ]

        if(silac.doub){
            exOut[["unique evidence L"]] <- intersect(evidence.unique.L.idx, evidence.ex.idx)
            exOut[["unique evidence H"]] <- intersect(evidence.unique.H.idx, evidence.ex.idx)
        }
        if(silac.trip){
            exOut[["unique evidence L"]] <- intersect(evidence.unique.L.idx, evidence.ex.idx)
            exOut[["unique evidence M"]] <- intersect(evidence.unique.M.idx, evidence.ex.idx)
            exOut[["unique evidence H"]] <- intersect(evidence.unique.H.idx, evidence.ex.idx)
        }

        ######################################
        # get the number of evidences
        ######################################
        exOut[["# peptides"]] <- dim(exOut[["evidence"]])[1]


        ######################################
        # get the number of nr peptides
        ######################################
        exOut[["# nr peptides"]] <- length( unique( exOut[["evidence"]][, grep("^Sequence$", colnames( exOut[["evidence"]] ), ignore.case=T) ]  )  )
        if(silac.doub){
            exOut[["# nr peptides L"]] <- length( unique( evidence[ exOut[["evidence L"]] , grep("^Sequence$", colnames( exOut[["evidence L"]] ), ignore.case=T)]  )  )
            exOut[["# nr peptides H"]] <- length( unique( evidence[ exOut[["evidence H"]] , grep("^Sequence$", colnames( exOut[["evidence H"]] ), ignore.case=T)]  )  )
        }
        if(silac.trip){
            exOut[["# nr peptides L"]] <- length( unique( evidence[ exOut[["evidence L"]] , grep("^Sequence$", colnames( exOut[["evidence L"]] ), ignore.case=T)]  )  )
            exOut[["# nr peptides M"]] <- length( unique( evidence[ exOut[["evidence M"]] , grep("^Sequence$", colnames( exOut[["evidence M"]] ), ignore.case=T)]  )  )
            exOut[["# nr peptides H"]] <- length( unique( evidence[ exOut[["evidence H"]] , grep("^Sequence$", colnames( exOut[["evidence H"]] ), ignore.case=T)]  )  )
        }

        ######################################
        # number of protein groups (unique evidences)
        ######################################
        exOut[["# protein groups"]] <- length( unique( exOut[["unique evidence"]][, grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence"]] ), ignore.case=T) ] ) )

        if(silac.doub){
            exOut[["# protein groups L"]] <- length( unique( evidence[ exOut[["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence L"]] ), ignore.case=T)]  )  )
            exOut[["# protein groups H"]] <- length( unique( evidence[ exOut[["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence H"]] ), ignore.case=T)]  )  )
        }
        if(silac.trip){
            exOut[["# protein groups L"]] <- length( unique( evidence[ exOut[["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence L"]] ), ignore.case=T)]  )  )
            exOut[["# protein groups M"]] <- length( unique( evidence[ exOut[["unique evidence M"]] , grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence M"]] ), ignore.case=T)]  )  )
            exOut[["# protein groups H"]] <- length( unique( evidence[ exOut[["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames( exOut[["unique evidence H"]] ), ignore.case=T)]  )  )
        }

        ######################################
        # store all results for current exp
        ######################################
        pepPerExp[[ex]] <- exOut

    } # end loop over exp


    ###################################################################################
    #   number of identified, non-redundant peptides/protein groups per experiment
    ###################################################################################

    # number of non-redundant peptides per experiment
    numbPeps <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides"]]  ))
    # number of all non-redundant peptides
    numbPeps <- c(numbPeps, length(unique(evidence[, grep("^Sequence$", colnames(evidence), ignore.case=T )])  ))

    # same for H/L/M
    if(silac.doub){
        numbPeps.L <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides L"]]  ))
        numbPeps.H <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides H"]]  ))
    }
    if(silac.trip){
        numbPeps.L <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides L"]]  ))
        numbPeps.M <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides M"]]  ))
        numbPeps.H <- unlist(lapply(pepPerExp, function(x) x[["# nr peptides H"]]  ))
    }

    ###################################################################################
    # number of protein groups per experiment
    ###################################################################################
    numbProts <- unlist( lapply(pepPerExp, function(x) x[["# protein groups"]]))

    # number of all protein groups
    numbProts <- c(numbProts, length( unique(unlist(strsplit( as.character(evidence[, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T) ]), ";")))))
    names(numbPeps)[length(numbPeps)] <- names(numbProts)[length(numbProts)] <- "# total"

    if(silac.doub){
        numbProts.L <- unlist( lapply(pepPerExp, function(x) x[["# protein groups L"]]))
        numbProts.H <- unlist( lapply(pepPerExp, function(x) x[["# protein groups H"]]))
    }
    if(silac.trip){
        numbProts.L <- unlist( lapply(pepPerExp, function(x) x[["# protein groups L"]]))
        numbProts.M <- unlist( lapply(pepPerExp, function(x) x[["# protein groups M"]]))
        numbProts.H <- unlist( lapply(pepPerExp, function(x) x[["# protein groups H"]]))
    }


    ####################################################################################
    # number of non-redundant peptides/protein groups identified in
    # 1, 2, ..., n experiments
    ####################################################################################

    # determine for each non-redundant peptide in which experiments it was found
    nrSequence <-  unique( as.character( evidence[, grep("^Sequence$", colnames(evidence), ignore.case=T)] ))

    # list to store the number of experiments per non-redundnat sequence
    pepsInExperiments <- vector("numeric", length(nrSequence) )
    names(pepsInExperiments) <- nrSequence

    #################################################
    # count the number of peptides per experiment
    #################################################
    pepsInExperiments <- .C( "peptidesInExperiments", as.character(nrSequence), as.integer(length(nrSequence)), as.character(evidence[, grep( "^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(dim(evidence)[1]), as.character(evidence[, grep("^Experiment$", colnames(evidence), ignore.case=T) ] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments)  )[[8]]

    # now count  the peptides for each number of experiments
    numbPepsNumbExperiments <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments == x  )  ))
    names(numbPepsNumbExperiments) <- 1:length(experiments)

    ####################################################################################
    #
    # count the number of heavy/medium/light peptides identified in 1,2,3,...,n exp
    #
    ###################################################################################
    if(silac.doub){
        ############
        # light
        ############
        # all all nr peptides
        nrSequence.L <-  unique( as.character( evidence[evidence.L.idx , grep("^Sequence$", colnames(evidence), ignore.case=T)] ))
        pepsInExperiments.L <- vector("numeric", length(nrSequence.L) )
        names(pepsInExperiments.L) <- nrSequence.L
        # count the peptides
        pepsInExperiments.L <- .C( "peptidesInExperiments", as.character(nrSequence.L), as.integer(length(nrSequence.L)), as.character(evidence[evidence.L.idx, grep("^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(length(evidence.L.idx)), as.character(evidence[evidence.L.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments.L)  )[[8]]

        # now count  the peptides for each number of experiments
        numbPepsNumbExperiments.L <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments.L == x  )  ))
        names(numbPepsNumbExperiments.L) <- 1:length(experiments)

        ############
        # heavy
        ############
        # all all nr peptides
        nrSequence.H <-  unique( as.character( evidence[evidence.H.idx , grep("^Sequence$", colnames(evidence), ignore.case=T)] ))
        pepsInExperiments.H <- vector("numeric", length(nrSequence.H) )
        names(pepsInExperiments.H) <- nrSequence.H
        # count the peptides
        pepsInExperiments.H <- .C( "peptidesInExperiments", as.character(nrSequence.H), as.integer(length(nrSequence.H)), as.character(evidence[evidence.H.idx, grep("^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(length(evidence.H.idx)), as.character(evidence[evidence.H.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments.H)  )[[8]]

        # now count  the peptides for each number of experiments
        numbPepsNumbExperiments.H <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments.H == x  )  ))
        names(numbPepsNumbExperiments.H) <- 1:length(experiments)

    }
     if(silac.trip){
        ############
        # light
        ############
        # all all nr peptides
        nrSequence.L <-  unique( as.character( evidence[evidence.L.idx , grep("^Sequence$", colnames(evidence), ignore.case=T)] ))
        pepsInExperiments.L <- vector("numeric", length(nrSequence.L) )
        names(pepsInExperiments.L) <- nrSequence.L
        # count the peptides
        pepsInExperiments.L <- .C( "peptidesInExperiments", as.character(nrSequence.L), as.integer(length(nrSequence.L)), as.character(evidence[evidence.L.idx, grep("^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(length(evidence.L.idx)), as.character(evidence[evidence.L.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments.L)  )[[8]]

        # now count  the peptides for each number of experiments
        numbPepsNumbExperiments.L <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments.L == x  )  ))
        names(numbPepsNumbExperiments.L) <- 1:length(experiments)


        ############
        # medium
        ############
        # all all nr peptides
        nrSequence.M <-  unique( as.character( evidence[evidence.M.idx , grep("^Sequence$", colnames(evidence), ignore.case=T)] ))
        pepsInExperiments.M <- vector("numeric", length(nrSequence.M) )
        names(pepsInExperiments.M) <- nrSequence.M
        # count the peptides
        pepsInExperiments.M <- .C( "peptidesInExperiments", as.character(nrSequence.M), as.integer(length(nrSequence.M)), as.character(evidence[evidence.M.idx, grep("^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(length(evidence.M.idx)), as.character(evidence[evidence.M.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments.M)  )[[8]]

        # now count  the peptides for each number of experiments
        numbPepsNumbExperiments.M <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments.M == x  )  ))
        names(numbPepsNumbExperiments.M) <- 1:length(experiments)


        ############
        # heavy
        ############
        # all all nr peptides
        nrSequence.H <-  unique( as.character( evidence[evidence.H.idx , grep("^Sequence$", colnames(evidence), ignore.case=T)] ))
        pepsInExperiments.H <- vector("numeric", length(nrSequence.H) )
        names(pepsInExperiments.H) <- nrSequence.H
        # count the peptides
        pepsInExperiments.H <- .C( "peptidesInExperiments", as.character(nrSequence.H), as.integer(length(nrSequence.H)), as.character(evidence[evidence.H.idx, grep("^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer(length(evidence.H.idx)), as.character(evidence[evidence.H.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(pepsInExperiments.H)  )[[8]]

        # now count  the peptides for each number of experiments
        numbPepsNumbExperiments.H <- unlist(lapply( 1:length(experiments), function(x) sum( pepsInExperiments.H == x  )  ))
        names(numbPepsNumbExperiments.H) <- 1:length(experiments)

    }


    ##################################################################################
    # do the same with protein groups
    # - consider only evidences that can
    #   be uniquly mapped to a single
    #   protein group
    ##################################################################################
    nrProtIds <- unique( as.character( evidence.unique[, grep("^Protein.Group.IDs$", colnames(evidence.unique), ignore.case=T) ]) )

    protsInExperiments <- vector( "numeric", length(nrProtIds) )
    names(protsInExperiments) <- nrProtIds

    ###########################
    # again a fu*** 'for' loop...
    ###########################
    protsInExperiments <- .C("peptidesInExperiments", as.character(nrProtIds), as.integer(length(nrProtIds)), as.character(evidence.unique[, grep("^Protein.Group.IDs$", colnames(evidence.unique), ignore.case=T)]), as.integer(dim(evidence.unique)[1]), as.character(evidence.unique[, grep("^Experiment$", colnames(evidence.unique), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments) )[[8]]

    numbProtsNumbExperiments <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments == x)  ) )
    names(numbProtsNumbExperiments) <- 1:length(experiments)

    if(silac.doub){
        #####################
        # light
        #####################
        nrProtIds.L <- unique( as.character( evidence[ evidence.unique.L.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
        protsInExperiments.L <- vector( "numeric", length(nrProtIds.L) )
        names(protsInExperiments.L) <- nrProtIds.L

        protsInExperiments.L <- .C("peptidesInExperiments", as.character(nrProtIds.L), as.integer(length(nrProtIds.L)), as.character(evidence[evidence.unique.L.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]), as.integer(length(evidence.unique.L.idx)), as.character(evidence[evidence.unique.L.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments.L) )[[8]]

        numbProtsNumbExperiments.L <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments.L == x)  ) )
        names(numbProtsNumbExperiments.L) <- 1:length(experiments)

        #####################
        # heavy
        #####################
        nrProtIds.H <- unique( as.character( evidence[ evidence.unique.H.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
        protsInExperiments.H <- vector( "numeric", length(nrProtIds.H) )
        names(protsInExperiments.H) <- nrProtIds.H

        protsInExperiments.H <- .C("peptidesInExperiments", as.character(nrProtIds.H), as.integer(length(nrProtIds.H)), as.character(evidence[evidence.unique.H.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]), as.integer(length(evidence.unique.H.idx)), as.character(evidence[evidence.unique.H.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments.H) )[[8]]

        numbProtsNumbExperiments.H <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments.H == x)  ) )
        names(numbProtsNumbExperiments.H) <- 1:length(experiments)
    }
    if(silac.trip){
        #####################
        # light
        #####################
        nrProtIds.L <- unique( as.character( evidence[ evidence.unique.L.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
        protsInExperiments.L <- vector( "numeric", length(nrProtIds.L) )
        names(protsInExperiments.L) <- nrProtIds.L

        protsInExperiments.L <- .C("peptidesInExperiments", as.character(nrProtIds.L), as.integer(length(nrProtIds.L)), as.character(evidence[evidence.unique.L.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]), as.integer(length(evidence.unique.L.idx)), as.character(evidence[evidence.unique.L.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments.L) )[[8]]

        numbProtsNumbExperiments.L <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments.L == x)  ) )
        names(numbProtsNumbExperiments.L) <- 1:length(experiments)

        #####################
        # medium
        #####################
        nrProtIds.M <- unique( as.character( evidence[ evidence.unique.M.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
        protsInExperiments.M <- vector( "numeric", length(nrProtIds.M) )
        names(protsInExperiments.M) <- nrProtIds.M

        protsInExperiments.M <- .C("peptidesInExperiments", as.character(nrProtIds.M), as.integer(length(nrProtIds.M)), as.character(evidence[evidence.unique.M.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]), as.integer(length(evidence.unique.M.idx)), as.character(evidence[evidence.unique.M.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments.M) )[[8]]

        numbProtsNumbExperiments.M <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments.M == x)  ) )
        names(numbProtsNumbExperiments.M) <- 1:length(experiments)

        #####################
        # heavy
        #####################
        nrProtIds.H <- unique( as.character( evidence[ evidence.unique.H.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
        protsInExperiments.H <- vector( "numeric", length(nrProtIds.H) )
        names(protsInExperiments.H) <- nrProtIds.H

        protsInExperiments.H <- .C("peptidesInExperiments", as.character(nrProtIds.H), as.integer(length(nrProtIds.H)), as.character(evidence[evidence.unique.H.idx, grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]), as.integer(length(evidence.unique.H.idx)), as.character(evidence[evidence.unique.H.idx, grep("^Experiment$", colnames(evidence), ignore.case=T)] ), as.character(experiments), as.integer(length(experiments)), as.integer(protsInExperiments.H) )[[8]]

        numbProtsNumbExperiments.H <- unlist( lapply(1:length(experiments), function(x) sum(protsInExperiments.H == x)  ) )
        names(numbProtsNumbExperiments.H) <- 1:length(experiments)

    }


    ##############################################################################################################
    #
    #                       determine the overlap between the different experiments
    #
    #
    ##############################################################################################################
    overlapPepMat <- matrix( 0, nrow=length(experiments) ,ncol=length(experiments), dimnames=list(experiments, experiments))
    overlapProtMat <- overlapProtMatRel <- overlapPepMatRel <- overlapPepMat

    if(silac.doub)
        overlapPepMat.L <- overlapPepMat.H <- overlapPepMatRel.L <- overlapPepMatRel.H <- overlapProtMat.L <- overlapProtMat.H <- overlapProtMatRel.L <- overlapProtMatRel.H <- overlapPepMat
    if(silac.trip)
        overlapPepMat.L <- overlapPepMat.H <- overlapPepMatRel.L <- overlapPepMatRel.H <- overlapProtMat.L <- overlapProtMat.H <- overlapProtMatRel.L <- overlapProtMatRel.H <- overlapPepMat.M <- overlapPepMatRel.M <- overlapProtMat.M <- overlapProtMatRel.M <- overlapPepMat

    ##########################
    # loop over matrix
    ##########################
    for(i in experiments){
          for(j in experiments){
                 # peptides: intersection of peptide ids of corresponding evidences
                 overlapPepMat[i,j] <- length( intersect(pepPerExp[[i]][["evidence"]][, grep( "^Peptide.ID$", colnames(pepPerExp[[i]][["evidence"]]), ignore.case=T ) ],   pepPerExp[[j]][["evidence"]][, grep( "^Peptide.ID$", colnames(pepPerExp[[j]][["evidence"]]), ignore.case=T )] )   )
                 # protein groups: intersection of protein group ids of UNIQUE evidences (unique to a protein group, i.e. there
                 # must be only one protein group id for each evidence...)
                 overlapProtMat[i,j] <- length( intersect(pepPerExp[[i]][["unique evidence"]][, grep("^Protein.Group.IDs$", colnames(pepPerExp[[i]][["unique evidence"]]), ignore.case=T)], pepPerExp[[j]][["unique evidence"]][, grep("^Protein.Group.IDs$", colnames(pepPerExp[[j]][["unique evidence"]]), ignore.case=T)]) )


                 ##########################
                 # double SILAC
                 ##########################
                 if(silac.doub){
                     # peptides
                     overlapPepMat.L[i,j] <- length( intersect( evidence[ pepPerExp[[i]][["evidence L"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["evidence L"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)]   )   )
                     overlapPepMat.H[i,j] <- length( intersect( evidence[ pepPerExp[[i]][["evidence H"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["evidence H"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)]   )   )

                     # protein groups
                     overlapProtMat.L[i,j] <- length( intersect(evidence[ pepPerExp[[i]][["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
                     overlapProtMat.H[i,j] <- length( intersect(evidence[ pepPerExp[[i]][["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
                 }
                 ##########################
                 # triple SILAC
                 ##########################
                 if(silac.trip){
                     # peptides
                     overlapPepMat.L[i,j] <- length( intersect( evidence[ pepPerExp[[i]][["evidence L"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["evidence L"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)]   )   )
                     overlapPepMat.H[i,j] <- length( intersect( evidence[ pepPerExp[[i]][["evidence H"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["evidence H"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)]   )   )
                     overlapPepMat.M[i,j] <- length( intersect( evidence[ pepPerExp[[i]][["evidence M"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["evidence M"]], grep("^Peptide.ID$", colnames(evidence), ignore.case=T)]   )   )

                     # protein groups
                     overlapProtMat.L[i,j] <- length( intersect(evidence[ pepPerExp[[i]][["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["unique evidence L"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
                     overlapProtMat.H[i,j] <- length( intersect(evidence[ pepPerExp[[i]][["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["unique evidence H"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
                     overlapProtMat.M[i,j] <- length( intersect(evidence[ pepPerExp[[i]][["unique evidence M"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)], evidence[ pepPerExp[[j]][["unique evidence M"]] , grep("^Protein.Group.IDs$", colnames(evidence), ignore.case=T)]) )
                 }


             }
    }
    ############################################
    # normalize to 1
    ############################################
    for(i in experiments){

        overlapPepMatRel[i, ] <- overlapPepMatRel[,i] <- overlapPepMat[i,]/max(overlapPepMat[i,], overlapPepMat[,i]  )
        overlapProtMatRel[i, ] <- overlapProtMatRel[,i] <- overlapProtMat[i,]/max(overlapProtMat[i,], overlapProtMat[,i]  )

        if(silac.doub){
            overlapPepMatRel.L[i, ] <- overlapPepMatRel.L[,i] <- overlapPepMat.L[i,]/max(overlapPepMat.L[i,], overlapPepMat.L[,i]  )
            overlapPepMatRel.H[i, ] <- overlapPepMatRel.H[,i] <- overlapPepMat.H[i,]/max(overlapPepMat.H[i,], overlapPepMat.H[,i]  )

            overlapProtMatRel.L[i, ] <- overlapProtMatRel.L[,i] <- overlapProtMat.L[i,]/max(overlapProtMat.L[i,], overlapProtMat.L[,i]  )
            overlapProtMatRel.H[i, ] <- overlapProtMatRel.H[,i] <- overlapProtMat.H[i,]/max(overlapProtMat.H[i,], overlapProtMat.H[,i]  )
        }
        if(silac.trip){
            overlapPepMatRel.L[i, ] <- overlapPepMatRel.L[,i] <- overlapPepMat.L[i,]/max(overlapPepMat.L[i,], overlapPepMat.L[,i]  )
            overlapPepMatRel.M[i, ] <- overlapPepMatRel.M[,i] <- overlapPepMat.M[i,]/max(overlapPepMat.M[i,], overlapPepMat.M[,i]  )
            overlapPepMatRel.H[i, ] <- overlapPepMatRel.H[,i] <- overlapPepMat.H[i,]/max(overlapPepMat.H[i,], overlapPepMat.H[,i]  )

            overlapProtMatRel.L[i, ] <- overlapProtMatRel.L[,i] <- overlapProtMat.L[i,]/max(overlapProtMat.L[i,], overlapProtMat.L[,i]  )
            overlapProtMatRel.M[i, ] <- overlapProtMatRel.M[,i] <- overlapProtMat.M[i,]/max(overlapProtMat.M[i,], overlapProtMat.M[,i]  )
            overlapProtMatRel.H[i, ] <- overlapProtMatRel.H[,i] <- overlapProtMat.H[i,]/max(overlapProtMat.H[i,], overlapProtMat.H[,i]  )
        }

    }



    ##################################################################
    #                          the ouptut
    ##################################################################
    if(silac.doub){
        output <- vector("list",13 )
        names(output) <- c( "experiment" , "# peptides", "# peptides L", "# peptides H", "# proteins", "# proteins L", "# proteins H", "overlap peptides", "overlap proteins", "overlap peptides L", "overlap peptides H", "overlap proteins L", "overlap proteins H")

        output[["# peptides L"]] <- numbPepsNumbExperiments.L
        output[["# peptides H"]] <- numbPepsNumbExperiments.H

        output[["# proteins L"]] <- numbProtsNumbExperiments.L
        output[["# proteins H"]] <- numbProtsNumbExperiments.H

        output[["overlap peptides L"]] <- overlapPepMat.L
        output[["overlap peptides H"]] <- overlapPepMat.H

        output[["overlap proteins L"]] <- overlapProtMat.L
        output[["overlap proteins H"]] <- overlapProtMat.H

    } else if(silac.trip){
        output <- vector("list",17 )
        names(output) <- c( "experiment" , "# peptides", "# peptides L", "# peptides M", "# peptides H", "# proteins", "# proteins L", "# proteins M", "# proteins H", "overlap peptides", "overlap proteins", "overlap peptides L", "overlap peptides M", "overlap peptides H", "overlap proteins L", "overlap proteins M", "overlap proteins H")

        output[["# peptides L"]] <- numbPepsNumbExperiments.L
        output[["# peptides M"]] <- numbPepsNumbExperiments.M
        output[["# peptides H"]] <- numbPepsNumbExperiments.H

        output[["# proteins L"]] <- numbProtsNumbExperiments.L
        output[["# proteins M"]] <- numbProtsNumbExperiments.M
        output[["# proteins H"]] <- numbProtsNumbExperiments.H


        output[["overlap peptides L"]] <- overlapPepMat.L
        output[["overlap peptides M"]] <- overlapPepMat.M
        output[["overlap peptides H"]] <- overlapPepMat.H

        output[["overlap proteins L"]] <- overlapProtMat.L
        output[["overlap proteins M"]] <- overlapProtMat.M
        output[["overlap proteins H"]] <- overlapProtMat.H


    } else{
        output <- vector("list",5 )
        names(output) <- c( "experiment" , "# peptides", "# proteins", "overlap peptides", "overlap proteins")
    }


    output[["experiment"]] <- pepPerExp
    output[["# peptides"]] <- numbPepsNumbExperiments
    output[["# proteins"]] <- numbProtsNumbExperiments
    output[["overlap peptides"]] <- overlapPepMat
    output[["overlap proteins"]] <- overlapProtMat

    return(output)
}

##################################################################
#               helper function
#
#  all multiple protein group IDs matching the same
#  evidence are splited and appendend to the end of
#  'evidence'´such that each row corresponds
#  to a single protein group ID
#
# changed:    20090507  implementation
#
##################################################################
splitMultipleProteinIDsInEvidence <- function(evidence)
{
    # get the number of all redundant protein group IDs
    pgNumb <-  .C("NumberProteinIDsEvidence", as.character(evidence[, "Protein.Group.IDs"]), as.integer(dim(evidence)[1] ), as.integer(0)  )[[3]]

    # now split
    res <- .C("splitProteinIDsInEvidence", as.character(evidence[, "Protein.Group.IDs"]), as.integer(dim( evidence )[1]), as.character( evidence[, "Experiment"] ), as.character(rep("", pgNumb)), as.character(rep("", pgNumb))  )


    return( cbind(res[[4]], res[[5]]) )
}

########################################################################################################################
#
#                                              checkSeparation
#
# The function determines the number of non-redundant peptides in each LC-MS run.
#
# arguments
#   - evidence                  - dataframe containing the 'evidence.txt'
#
# value
#   - list containing the following elements:
#
#
#
# changed: 20090319      - implementation of basic functionality
#          20090507      - replaced the for-loop by a C function
#          20101020      - switched to function 'fancyBarplot'
#          20120710      - compatible to 1.2.7.5 (case of column names doesn't matter anymore)
########################################################################################################################
checkSeparation <- function(evidence, file="checkFractionation.pdf" )
{

    # determine the number of MS runs (based
    # on the filenames )
    rawFiles <- unique( as.character(  evidence[, grep( "^Raw.File$", colnames(evidence), ignore.case=T)] ))

    ##############################################################
    #  determine the number of identified peptides (NOT UNIQUE)
    #  within each fraction
    ##############################################################
    numbPepPerFrac <- vector("numeric", length(rawFiles))
    names(numbPepPerFrac) <- rawFiles

    # list to store the peptide sequences of each well
    seqPepPerFrac<- vector("list", length(rawFiles))
    names(seqPepPerFrac) <- rawFiles

    for(rF in rawFiles)
    {
      # get the sequences
      seqPepPerFrac[[rF]] <- as.character(evidence[ which(as.character(evidence[, grep( "^Raw.File$", colnames(evidence), ignore.case=T)]) == rF ) ,  grep( "^Sequence$", colnames(evidence), ignore.case=T)])

      # count the peptides
      #numbPepPerFrac[rF] <- sum( as.character(evidence[, "Raw.File"]) == rF )
      numbPepPerFrac[[rF]] <- length( seqPepPerFrac[[rF]] )
    }

    ###########################################################################
    #  determine for each peptide in how many fractions it has been identified
    #
    ###########################################################################

    # non-redundant peptide sequences
    nrSeq <- unique(as.character(evidence[, "Sequence"]))

    pepInFrac <- vector("numeric", length(nrSeq)  )
    names(pepInFrac) <- nrSeq


    pepInFrac <- .C("peptidesInExperiments", as.character(nrSeq), as.integer(length(nrSeq)), as.character( evidence[, grep( "^Sequence$", colnames(evidence), ignore.case=T)] ), as.integer( dim(evidence)[1] ), as.character(evidence[, grep( "^Raw.File$", colnames(evidence), ignore.case=T)]), as.character(rawFiles), as.integer(length(rawFiles)), as.integer(pepInFrac)  )[[8]]


    # determine the number
    numbPepInFrac <- unlist( lapply(1:length(rawFiles), function(x) sum(pepInFrac == x))  )
    names(numbPepInFrac) <- 1:length(rawFiles)


    #############################################################################
    #                        plot
    #############################################################################

    fancyBarplot(numbPepInFrac, main="Number of Non-Redundant Peptides per Number of Raw Files", xlab="# raw files", ylab="# non-redundant peptides")


}



#########################################################################################################################
## Plot a 2-Way, 3-Way or 4-Way Venn Diagram ##
###############################################
## Author: Thomas Girke
## Last update: Nov 6, 2008
## Utility: Plots a non-proportional 2-, 3- or 4-way venn diagram based on overlaps among data sets (vectors)
## Detailed instructions for running this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn

## Define venndiagram function
venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE, ...) {
	## Remove duplicates and NA fields in x, y, z and w
	if(unique==T) {
		x <- unique(x); x <- as.vector(na.omit(x))
		y <- unique(y); y <- as.vector(na.omit(y))
		if(!missing("z")) {
			z <- unique(z); z <- as.vector(na.omit(z))
		}
		if(!missing("w")) {
			w <- unique(w); w <- as.vector(na.omit(w))
		}
	}

	## Check valid type selection
	if(!type %in% c("2", "2map", "3", "3map", "4", "4map", "4el", "4elmap")) {
		return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
	}

	## Plot a 2-way venn diagram
	if(type=="2") {
		# Define ovelap queries
		q1 <- x[x %in% y]
		q2 <- x[!x %in% y]
		q3 <- y[!y %in% x]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep="")} else {mysub <- ""}
		if(plot==T) {
			## Plot the 2-way venn diagram
			symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 2-way mapping venn diagram
	if(type=="2map") {
		olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
	}

	## Plot a 3-way venn diagram
	if(type=="3") {
		## Define ovelap queries
		q1 <- x[x %in% y & x %in% z]
		q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
		q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
		q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
		q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep="")} else { mysub <- "" }
		if(plot==T) {
			## Plot the 3-way venn diagram
			symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 3-way mapping venn diagram
	if(type=="3map") {
		olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
		symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
	}

	## Overlap queries for 4-way venn diagram
	if(type=="4" | type=="4el" | type=="4elmap") {
		## Define ovelap queries
		xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
		q1 <- xy[xy %in% zw]
		q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
		q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
		q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
		q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
		q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w]
		q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
		q10 <- x[!x %in% c(y,z,w)]
		q11 <- y[!y %in% c(x,z,w)]
		q12 <- z[!z %in% c(x,y,w)]
		q13 <- w[!w %in% c(x,y,z)]
		q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
		q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)

		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep="") } else { mysub <- "" }

	## Plot 4-way venn diagram as circles
		if(plot==T & type=="4") {
			symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
			text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
			text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
		}

	## Plot 4-way venn diagram as ellipses
	if(plot==T & (type=="4el" | type=="4elmap")) {
		olDF <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=countDF$count)
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(lines=lines, olDF, title=title, labels=labels, sub=mysub, main, lcol=lcol, tcex=1.3, ...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=title, sub=mysub, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, ...)
			text(olDF[1:15,1], olDF[1:15,2], olDF[1:15,3], col=tcol, ...)
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels, col=lcol, ...)
			close.screen(all=TRUE)
		}
		## Plot 4-way ellipse venn diagram
		if(type=="4el") {
			ellipseVenn(olDF=olDF, lcol=lcol, lines=lines, labels=labels, title=title, ...)
		}

		## Plot 4-way ellipse mapping venn diagram
		if(type=="4elmap") {
			olDFdebug <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=paste("q", 1:15, sep=""), ...)
			ellipseVenn(olDF=olDFdebug, lcol=lcol, lines=lines, labels=paste(labels, "=", c("x","y","z","w")), title="Mapping Venn Diagram", ...)
		}
	}

	## Return query list
	return(qlist)
	}

	## Plot 4-way circle mapping venn diagram
	if(type=="4map") {
		olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
		text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
	}

}

## Generate overlap reports
olReport <- function(qlist=qlist, missing=".", type=1) {
	## Check valid type selection
	if(!type %in% c(1, 2, 3, 4)) {
		return("Error: the 'type' argument can only be set to the values: 1, 2 or 3.")
	}

	## Output data frame with overlap keys in separate columns
	if(type==1) {
		ids <- sort(unique(as.vector(unlist(qlist))))
		qDF <- matrix(ids, nrow=length(ids), ncol=length(qlist), dimnames=list(1:length(ids), names(qlist)))
		lqDF <- as.data.frame(lapply(names(qlist), function(x) qDF[,x] %in% qlist[[x]]))
		colnames(lqDF) <- colnames(qDF)
		lqDF <- as.matrix(lqDF)
		qDF[!lqDF] <- missing
		qDF <- data.frame(IDs=ids, qDF)
		return(qDF)
	}

	## Output data frame with overlap section numbers (qNo) in one column
	if(type==3) {
		collapsedDF <- data.frame(IDs=as.vector(unlist(qlist)), qNo=rep(names(qlist), sapply(qlist, length)))
		collapsedDF <- collapsedDF[order(collapsedDF$IDs), ]
		rownames(collapsedDF) <- 1:length(collapsedDF[, 1])
		return(collapsedDF)
	}

	## Output data frame with overlap counts
	if(type==2) {
		qStat <- data.frame(count=sapply(qlist, length))
		return(qStat)
	}

	## Output presence-absence matrix
	if(type==4) {
		ids <- sort(unique(as.vector(unlist(qlist))))
		qDF <- matrix(ids, nrow=length(ids), ncol=length(qlist), dimnames=list(1:length(ids), names(qlist)))
		lqDF <- as.data.frame(lapply(names(qlist), function(x) qDF[,x] %in% qlist[[x]]))
		colnames(lqDF) <- colnames(qDF)
		lqDF <- as.matrix(lqDF)
		lqDF[lqDF] <- 1
		lqDF[!lqDF] <- 0
		rownames(lqDF) <- ids
		lqDF <- lqDF[names(rev(sort(rowSums(lqDF)))),]
		return(lqDF)
	}
}

##########################################################################################################
#
#
# changelog:  20100929 implementation
#
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){

    out <- vector( "character", length(color) )

    for(col in 1:length(color)){

        col.rgb <- col2rgb(color[col])

        out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)

    }
    return(out)
}


##########################################################################################################
#
#                  capitalize all words in 's'
#
# function taken from R-help of 'tolower'/'toupper'
#
# changelog: 20120329 implentation
##########################################################################################################
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                  {s <- substring(s,2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

##########################################################################################################
#
#
#
# changelog: 20120329 implementation
##########################################################################################################
fixColumnNames <- function(colnames){

    # capitalize each part separated by '.'
    fix <- unlist(lapply(colnames, function(x) paste(capwords(unlist(strsplit(x, '\\.'))), collapse='.') ))

    # remove trailing dots
    fix <- sub("\\.*$", "", fix)

    return(fix)
}
