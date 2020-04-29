
require('argparser')

# Create a parser
p=arg_parser("Run LOLA")

# Add command line arguments
p=add_argument(p, "--database", help="database", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/baf_perturbseq/results/2020-03-26/data/nm/t1/resources/regions/")

p=add_argument(p, "--LOLAdir", help="usually LOLACore/hg38", type="character",
default="LOLACore/hg38")

p=add_argument(p, "--zscore", help="z score threshold", type="numeric",default=3)

p=add_argument(p, "--genecoords", help="", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/data/genomes/hg38/promoters/hg38.TSS.size5000.bed.gz")

p=add_argument(p, "--outpref", help="", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/baf_perturbseq/results/2020-03-26/results/2020-04-06_NMF_genes/test.LOLA2")

p=add_argument(p, "--subsetfilelist", help="", type="character",
default="/ahg/regevdata/projects/Cell2CellCommunication/baf_perturbseq/results/2020-03-26/results/2020-04-06_NMF_genes/nmf_filelist.txt")

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

#parameters
#===============
full_set_f=args$genecoords
subset_file_list_f=args$subsetfilelist
ZCORE_THRESHOLD=args$zscore

#the total set
regionSet = readBed(full_set_f)

f_df=read.table(subset_file_list_f)
colnames(f_df)=c('full_path')

f_df['name']=t(as.data.frame(strsplit(as.character(f_df[,1]),'NMF_genes/')))[,2]
print(head(f_df))

i=1
cur_path=f_df[i,'full_path']
cur_sample=f_df[i,'name']

cur_df=read.table(as.character(cur_path),header=TRUE)
print(head(cur_df))
genes=as.character(cur_df[which(as.numeric(as.character(cur_df[,'zscore']))>=ZCORE_THRESHOLD),
                                'gene'])
print(genes)

regionSet_here=regionSet[genes]
print(regionSet_here)

require(GenomicRanges)
universeSets = GRangesList(regionSet)
userSets=GRangesList(regionSet_here)
locResults = runLOLA(userSets, regionSet,regionDB,cores=1,direction=args$direction)

locResults=as.data.frame(locResults)
rownames(locResults)=locResults[,'filename']

results=data.frame(matrix(0,ncol = dim(f_df)[1], nrow = dim(locResults)[1]))
colnames(results)=f_df$name
rownames(results)=locResults$filename
print(dim(results))
results[,f_df[i,'name']]=locResults[rownames(results),'pValueLog']

for (i in c(2:(dim(results)[2]))){
    cur_path=f_df[i,'full_path']
    cur_sample=f_df[i,'name']

    cur_df=read.table(as.character(cur_path),header=TRUE)
    print(head(cur_df))
    genes=as.character(cur_df[which(as.numeric(as.character(cur_df[,'zscore']))>=ZCORE_THRESHOLD),
                                    'gene'])
    print(genes)

    regionSet_here=regionSet[genes]
    print(regionSet_here)

    require(GenomicRanges)
    universeSets = GRangesList(regionSet)
    userSets=GRangesList(regionSet_here)
    locResults = runLOLA(userSets, regionSet,regionDB,cores=1)

    locResults=as.data.frame(locResults)
    rownames(locResults)=locResults[,'filename']
    
    results[,f_df[i,'name']]=locResults[rownames(results),'pValueLog']
}

#now, save the results, and the locResults, so we can make nice names after
write.table(results,paste(args$outpref,'.LOLA.results.txt',sep=''),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)

write.table(locResults[,c('filename','cellType','antibody','description')],
            paste(args$outpref,'.LOLA.anno.txt',sep=''),sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
