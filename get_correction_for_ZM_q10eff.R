# get Tcorr in the first year for ZM to calculate a compensation factor

dat <- readLines('data/log.txt')

t_init <- grep('Time\\: 0.000000e\\+00', dat)[1]
t_end <- grep('Time\\: 3.660000e\\+02', dat)[1]-1

year1 <- dat[t_init:t_end]
albi <- year1[grep('ALBI', year1)]

albi_vals <- as.numeric(gsub(', Tcorr:.*', '', gsub('.*Tscalar:', '', albi)))

mean(albi_vals) # 0.82377

hist(albi_vals)

# so correction for mum should be 1/0.82377 = 1.21 (sounds really low)


10*0.82377*1.213931
