
require('argparser')

# Create a parser
p=arg_parser("Run LOLA")

# Add command line arguments
p=add_argument(p, "--database", help="database", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/baf_perturbseq/results/2020-03-26/data/nm/t1/resources/regions/")

p=add_argument(p, "--LOLAdir", help="usually LOLACore/hg38", type="character",
default="LOLACore/hg38")

p=add_argument(p, "--set", help="set of regions of interest", type="character")

p=add_argument(p, "--universe", help="set of regions of interest", type="character")

p=add_argument(p, "--outpref", help="", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/baf_perturbseq/results/2020-03-26/results/2020-04-06_NMF_genes/test.LOLA2")

p=add_argument(p, "--direction", help="can be enrichment of depletion", type="character",
default="enrichment")

# Parse the command line arguments
args=parse_args(p)

#imports
#=======
require(LOLA)
require(simpleCache)

#set database
#============
setwd(args$database)
regionDB = loadRegionDB(args$LOLAdir)

#the total set
regionSet = readBed(args$universe)

#region subset
regionSet_here=readBed(args$set)

require(GenomicRanges)
universeSets = GRangesList(regionSet)
userSets=GRangesList(regionSet_here)
locResults = runLOLA(userSets, regionSet,regionDB,cores=1,direction=args$direction)

locResults=as.data.frame(locResults)
rownames(locResults)=locResults[,'filename']

#now, save the results, and the locResults, so we can make nice names after
write.table(locResults,paste(args$outpref,'.LOLA.results.txt',sep=''),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)

write.table(locResults[,c('filename','cellType','antibody','description')],
            paste(args$outpref,'.LOLA.anno.txt',sep=''),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
