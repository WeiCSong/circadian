trait=alldisturbance
your_genomewide_sumstats="/dssg/home/acct-bioxsyy/bioxsyy-user2/ukb/gwas/allbart.disturbance.summary"

tissue="Liver"
for chr in {1..22}
do
Rscript fusion_twas-master/FUSION.assoc_test.R --sumstats $your_genomewide_sumstats --weights TCSC/weights/allEUR_tissues/v8_allEUR_${tissue}_blup.pos --weights_dir TCSC_weight_files/weights --ref_ld_chr LDREF/1000G.EUR. --chr $chr --out $(tissue}/${trait}.${chr}.dat
done
