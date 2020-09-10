library(data.table)
library(bigsnpr)

######################
#Paths
######################
bgen_path = "~/ukbb/v3/imputation/ukb_imp_chr22_v3.bgen"
bgi_dir = "~/ukbb/v3/imputation/"
bk_file = "~/ukb_hm_in_icbp_groups_chr22"

snp_list_path = "~/ukb_hm_snps_in_sumstats_chr22.csv"
#individuals_path = "~/group_eids_index_shared.csv" #based on .fam, contains indices outside the range of the .sample, leads to incomplete data
individuals_path = "~/group_eids_index_shared_from_sample.csv" #based on .sample, better!

######################
#Individuals and SNPs
######################

E = fread(individuals_path)
snp_list = fread(snp_list_path,header=FALSE)$V1

######################
#Make the RDS from the BGEN
######################
(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

if(file.exists(paste0(bk_file,".bk"))){
  print("Removing existing, old bk file:")
  print(bk_file)
  system(paste0("rm ",paste0(bk_file,".bk")))
}

system.time(
  rds_path <- bigsnpr::snp_readBGEN(
    bgenfiles = bgen_path,
    list_snp_id = list(snp_list),
    ind_row = E$index,
    backingfile = bk_file,
    bgi_dir = bgi_dir,
    ncores = NCORES
  )
)

######################
#Load the RDS
######################
obj.bigSNP <- snp_attach(rds_path)
G   <- obj.bigSNP$genotypes

######################
#Check for missing data
######################

y = big_counts(G,ind.col=1:ncol(G),byrow=TRUE)
individual_NAs = table(y[4,])
N_missing = sum(individual_NAs[names(individual_NAs)!="0"])
print("Number of individuals with missing dosages:")
print(N_missing)
