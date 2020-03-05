DB = "hg38"
#DB = "hg19"
TMPDIR = "/scratch/tmp"
# set to None to not use cluster
#CLUSTER=None
CLUSTER="ku"
TARGETPERC = 99 # percentile, for top50 consensus building
rmskTable = "rmskHmmer2020" # the rmsk table, usually it's rmsk, but you can specify a different one here
RMTORB = "data/rmToDfam.tab"
RBTORM = "data/rbToRm.tab"
