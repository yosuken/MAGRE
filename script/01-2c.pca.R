
## args
nargs = 2
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != nargs){
	stop(paste0("num of argument should be ", nargs, " (", length(args), " given)"))
}

odir = args[1] ## output directory
type = args[2] ## type (coverge or tetranucleotide)

## variables
N = 10 ## parse PC1 to PC10


### files
pref  = paste(odir, "/", type, sep="")
fin   = paste(pref, ".tsv", sep="")
fout1 = paste(pref, ".pca_coordinates", sep="")
fout2 = paste(pref, ".pca_importance", sep="")
fout3 = paste(pref, ".pc1_zscore", sep="")

### matrix
m  = as.matrix(read.delim(fin, row.names=1, as.is = T))
if (nrow(m) == 1) { break } ### [!!!] if n_contig == 1 then NO output files will be generated

### extract only non-zero columns
### ref: https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
m   = m[, colSums(m != 0, na.rm = T) > 0]
## [!!!] m[, 1] returns 'vector', not 'matrix'

if (is.vector(m)) { break } ### [!!!] if nonzero sample is only one (e.g., SRS1016432_bin.52) then NO output files will be generated

### normalize
if (type == "tetranucleotide"){
	### normalize tetranucleotide (4mer count --> 4mer fraction in each contig)
	m  = t(apply(m, 1, function(i){ return(i/sum(i)) })) ## scale rowwise (to be sum of row is 1), then TRANSPOSE
} else if (type == "coverage"){
	### normalize coverage (coverage --> coverage fraction in each sample)
	m  =   apply(m, 2, function(i){ return(i/sum(i)) })  ## scale colwise (to be sum of col is 1)
}

### [!!!] error with scale=T of prcomp()
# m1  = prcomp(m, scale=T) ## --> error when fin == ERS492926_bin.70.tsv (Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd')
### error in ERS492926 if `scale=T` in prcomp() --> Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
# Calls: prcomp -> prcomp.default -> svd -> La.svd
# Execution halted

### PCA by prcomp()
m1  = prcomp(m) ## default: scale=F, center=T

### extract primary dimentions
max = dim(m1$x)[[2]]
n   = min(N, max)
o1  = m1$x[,1:n]
o2  = summary(m1)$importance[,1:n]

### convert to z-score
o3  = scale(m1$x[,1:1])
colnames(o3) = c("PC1_zscore")

write.table(o1, fout1, col.names=NA, quote=F, sep="\t")
write.table(o2, fout2, col.names=NA, quote=F, sep="\t")
write.table(o3, fout3, col.names=NA, quote=F, sep="\t")
