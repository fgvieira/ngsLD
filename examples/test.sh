SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd


##### Clean up
rm -f testLD*


N_IND=24
N_SITES=10000


##### Genotypes
cat $SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testLD_T.geno.gz
zcat testLD_T.geno.gz | awk 'BEGIN{cnt=1} pos > 10000 {pos=0; cnt++} {pos+=int(rand()*1000+1); print $1"_"cnt"\t"pos}' > testLD.pos
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_T.geno.gz --pos testLD.pos --max_kb_dist 20 --min_maf 0.05                               | sort -k 1,1V -k 2,2V > testLD_T.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_T.geno.gz --pos testLD.pos --max_kb_dist 20 --min_maf 0.05 --rnd_sample 0.5 --seed 12345 | sort -k 1,1V -k 2,2V > testLD_Tr.ld

##### Genotypes' likelihood and posterior probabilities
# Binary, normal-scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 32 -out testLD_32
gunzip -f testLD_32.geno.gz
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05                                              | sort -k 1,1V -k 2,2V > testLD_32.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --call_geno                                  | sort -k 1,1V -k 2,2V > testLD_32-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_32.geno   --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1V -k 2,2V > testLD_32-CGf.ld

# Text, normal scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGeno 8 -out testLD_8
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05                                              | sort -k 1,1V -k 2,2V > testLD_8.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --call_geno                                  | sort -k 1,1V -k 2,2V > testLD_8-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_8.geno.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1V -k 2,2V > testLD_8-CGf.ld

# LD prunning
../scripts/prune_graph.pl --i testLD_8.ld --max_dist 5 --min_weight 0.95 --weight_field 6 --print_excl testLD_pruned.id.gz | sort -k 1,1V > testLD_unlinked.id
gunzip testLD_pruned.id.gz


##### Check MD5
rm -f *.arg
TMP=`mktemp`
md5sum testLD* | fgrep -v '.gz' | sort -k 2,2 > $TMP
if diff $TMP test.md5 > /dev/null
then
    echo "ngsLD: All tests OK!"
else
    echo "ngsLD: test(s) failed!"
fi
