#!/bin/bash

# run chopchop with command similar to web-server
# https://bitbucket.org/valenlab/chopchop/src/master/

genome=${1}
gene_name=${2}
outdir=${3:-"chopchop_results"}

# defaults
limitPrintResults=100


# create result folder
mkdir -p ${outdir}

#######################
# standard run command

function standard_run (){
	chopchop.py \
		--jsonVisualize \
		--padSize \
		--MODE 1 \
		--PAM NGG \
		--maxMismatches 3 \
		--guideSize 20 \
		--scoringMethod DOENCH_2016 \
		--fivePrimeEnd NN \
		--backbone AGGCTAGTCCGT \
		--replace5P GG \
		--target CODING \
		--enzymeCo N \
		--minResSiteLen 4 \
		--primer3options 'PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60' \
		--primerFlanks 290 \
		--guidePadding 20 \
		--rm1perfOff \
		--outputDir ${outdir} \
		--genome ${genome} \
		${gene_name} \
		> ${outdir}/results.txt \
		2> ${outdir}/python.err
	}

#######################
# ONT enrichment command
#    {"fastaInput":"","forSelect":"nanoporE","geneInput":"ICAM1","isIsoform":false,"opts":["-J","-BED","-GenBank","-G","hg38","-filterGCmin","20","-filterGCmax","80","-filterSelfCompMax","0","-consensusUnion","-t","WHOLE","-n","N","-R","4","-3","PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60","-P","-A","290","-DF","300","-a","20","-T","1","-g","20","-scoringMethod","DOENCH_2016","-f","NN","-v","3","-M","NGG","-BB","AGGCTAGTCCGT","--nonOverlapping"]}

function enrich_run() {
	chopchop.py \
		--jsonVisualize \
		--MODE 1 \
		--PAM NGG \
		--maxMismatches 3 \
		--guideSize 20 \
		-scoringMethod DOENCH_2016 \
		--fivePrimeEnd NN \
		--backbone "AGGCTAGTCCGT" \
		--target WHOLE \
		--enzymeCo N \
		--minResSiteLen 4 \
		--primer3options 'PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60' \
		--primerFlanks 290 \
		--guidePadding 20 \
		--BED \
		--GenBank \
		--offtargetsTable \
		--filterGCmin 20 \
		--filterGCmax 80 \
		--filterSelfCompMax 0 \
		--consensusUnion \
		--makePrimers \
		--displaySeqFlanks 300 \
		--nonOverlapping" \
		--limitPrintResults ${limitPrintResults} \
		--outputDir ${outdir} \
		--genome ${genome} \
		--targets ${gene_name} \
		> ${outdir}/results.txt \
		2> ${outdir}/python.err
	}

# choose method to run
standard_run
# enrich_run


########################################################################################
#       -Target TARGETS, --targets TARGETS
#                         Target genes or regions
#   -r, --gRVD            Use RVD 'NN' instead of 'NH' for guanine nucleotides.
#                         'NH' appears to be more specific than 'NN' but the
#                         choice depends on assembly kit.
#   -D DATABASE, --database DATABASE
#                         Connect to a chopchop database to retrieve gene:
#                         user_name:passwd@host/database
#   -e EXON_NUMBER, --exon EXON_NUMBER
#                         Comma separated list of exon indices. Only find sites
#                         in this subset.
#   -TDP TARGETDOWNSTREAMPROMOTER, --targetDownstreamPromoter TARGETDOWNSTREAMPROMOTER
#                         how many bp to target downstream of TSS
#   -TUP TARGETUPSTREAMPROMOTER, --targetUpstreamPromoter TARGETUPSTREAMPROMOTER
#                         how many bp to target upstream of TSS
#   -G GENOME, --genome GENOME
#                         The genome to search.
#   -g GUIDE_SIZE, --guideSize GUIDE_SIZE
#                         The size of the guide RNA.
#   -c, --scoreGC         Score GC content. True for CRISPR, False for TALENs.
#   -SC, --noScoreSelfComp
#                         Do not penalize self-complementarity of CRISPR.
#   -BB BACKBONE, --backbone BACKBONE
#                         Penalize self-complementarity versus backbone regions
#                         (comma-separated list, same strand as guide). Requires
#                         -C.
#   -R5 REPLACE_5P, --replace5P REPLACE_5P
#                         Replace bases from 5' end (with e.g. 'GG')
#   -t TARGETREGION, --target TARGETREGION
#                         Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE.
#                         Default is CODING.
#   -T {1,2,3,4}, --MODE {1,2,3,4}
#                         Set mode (int): default is Cas9 = 1, Talen = 2, Cpf1 =
#                         3, Nickase = 4
#   -taleMin TALEMIN, --taleMin TALEMIN
#                         Minimum distance between TALENs. Default is 14.
#   -taleMax TALEMAX, --taleMax TALEMAX
#                         Maximum distance between TALENs. Default is 20.
#   -nickaseMin NICKASEMIN, --nickaseMin NICKASEMIN
#                         Minimum distance between TALENs. Default is 10.
#   -nickaseMax NICKASEMAX, --nickaseMax NICKASEMAX
#                         Maximum distance between TALENs. Default is 31.
#   -offtargetMaxDist OFFTARGETMAXDIST, --offtargetMaxDist OFFTARGETMAXDIST
#                         Maximum distance between offtargets for Nickase.
#                         Default is 100.
#   -f FIVEPRIMEEND, --fivePrimeEnd FIVEPRIMEEND
#                         Specifies the requirement of the two nucleotides 5'
#                         end of the CRISPR guide: A/C/G/T/N. Default: NN.
#   -n ENZYME_CO, --enzymeCo ENZYME_CO
#                         The restriction enzyme company for TALEN spacer.
#   -R MINRESSITELEN, --minResSiteLen MINRESSITELEN
#                         The minimum length of the restriction enzyme.
#   -v MAX_MISMATCHES, --maxMismatches MAX_MISMATCHES
#                         The number of mismatches to check across the sequence.
#   -m MAX_HITS, --maxOffTargets MAX_HITS
#                         The maximum number of off targets allowed.
#   -M PAM, --PAM PAM     The PAM motif.
#   -o OUTPUT_DIR, --outputDir OUTPUT_DIR
#                         The output directory. Default is the current
#                         directory.
#   -F, --fasta           Use FASTA file as input rather than gene or genomic
#                         region.
#   -p PADSIZE, --padSize PADSIZE
#                         Extra bases searched outside the exon. Defaults to the
#                         size of the guide RNA for CRISPR and TALEN + maximum
#                         spacer for TALENS.
#   -P, --makePrimers     Designes primers using Primer3 to detect mutation.
#   -3 PRIMER3OPTIONS, --primer3options PRIMER3OPTIONS
#                         Options for Primer3. E.g. 'KEY1=VALUE1,KEY2=VALUE2'
#   -A PRIMERFLANKS, --primerFlanks PRIMERFLANKS
#                         Size of flanking regions to search for primers.
#   -DF DISPLAYSEQFLANKS, --displaySeqFlanks DISPLAYSEQFLANKS
#                         Size of flanking regions to output sequence into
#                         locusSeq_.
#   -a GUIDEPADDING, --guidePadding GUIDEPADDING
#                         Minimum distance of primer to target site.
#   -O LIMITPRINTRESULTS, --limitPrintResults LIMITPRINTRESULTS
#                         The number of results to print extended information
#                         for. Web server can handle 4k of these.
#   -w, --uniqueMethod_Cong
#                         A method to determine how unique the site is in the
#                         genome: allows 0 mismatches in last 15 bp.
#   -J, --jsonVisualize   Create files for visualization with json.
#   -nonO, --nonOverlapping
#                         Will not produce overlapping guides, saves time, and
#                         recommended for permissive PAMs (e.g. Cas13d).
#   -scoringMethod {XU_2015,DOENCH_2014,DOENCH_2016,MORENO_MATEOS_2015,CHARI_2015,G_20,KIM_2018,ALKAN_2018,ZHANG_2019,ALL}, --scoringMethod {XU_2015,DOENCH_2014,DOENCH_2016,MORENO_MATEOS_2015,CHARI_2015,G_20,KIM_2018,ALKAN_2018,ZHANG_2019,ALL}
#                         Scoring used for Cas9 and Nickase. Default is G_20. If
#                         a method fails to give scores, CHOPCHOP will output 0
#                         instead of terminating.
#   -repairPredictions {mESC,U2OS,HEK293,HCT116,K562}, --repairPredictions {mESC,U2OS,HEK293,HCT116,K562}
#                         Use inDelphi from Shen et al 2018 to predict repair
#                         profiles for every guideRNA, this will make
#                         .repProfile and .repStats files
#   -rm1perfOff, --rm1perfOff
#                         For fasta input, don't score one off-target without
#                         mismatches.
#   -isoforms, --isoforms
#                         Search for offtargets on the transcriptome.
#   -filterGCmin FILTERGCMIN, --filterGCmin FILTERGCMIN
#                         Minimum required GC percentage. Default is 0.
#   -filterGCmax FILTERGCMAX, --filterGCmax FILTERGCMAX
#                         Maximum allowed GC percentage. Default is 100.
#   -filterSelfCompMax FILTERSELFCOMPMAX, --filterSelfCompMax FILTERSELFCOMPMAX
#                         Maximum acceptable Self-complementarity score. Default
#                         is -1, no filter.
#   -consensusUnion, --consensusUnion
#                         When calculating consensus sequence from multiple
#                         isoforms default uses intersection. This option
#                         specifies union of isoforms.
#   -BED, --BED           Create results as BED file, can be used for
#                         integration with UCSC.
#   -GenBank, --GenBank   Create results as GenBank file, sequence of targeted
#                         region with introns is included.
#   -offtargetsTable, --offtargetsTable
#                         Create .tsv table with off-targets. Not all off-
#                         targets will be reported when early stopping will work
#                         on a guide! Limited also to CRISPR mode only and
#                         limited by --limitPrintResults option.