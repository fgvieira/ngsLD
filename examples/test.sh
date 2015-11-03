SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd


##### Clean up
rm -f testLD*


N_IND=24
N_SITES=10000


##### TRUE genotypes
cat $SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testLD_T.geno.gz
zcat testLD_T.geno.gz | awk '{pos+=int(rand()*10000+1); print $1"\t"pos}' > testLD.pos
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_T.geno.gz --pos testLD.pos --max_dist 10 | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_T.ld

##### Genotypes
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 2 -postCutoff 0.95 -out testLD_2
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_2.geno.gz --pos testLD.pos --max_dist 10 | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_2.ld

##### Genotypes' likelihood and posterior probabilities
# Binary, log-scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 32 -out testLD_32
gunzip -f testLD_32.geno.gz
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_dist 10                                              | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_32.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_dist 10 --call_geno                                  | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_32-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_dist 10 --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_32-CGf.ld

# Text, normal scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 8 -out testLD_8
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_dist 10                                              | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_8.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_dist 10 --call_geno                                  | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_8-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_dist 10 --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1 -k 2,2g -k 3,3 -k 4,4g > testLD_8-CGf.ld



##### Check MD5
rm -f *.arg
md5sum testLD* | fgrep -v '.gz' | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsDist: All tests OK!"
else
    echo "ngsDist: test(s) failed!"
fi
