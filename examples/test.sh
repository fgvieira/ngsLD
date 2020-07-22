SIM_DATA=../../ngsSim/examples
ANGSD=../../angsd


##### Clean up
rm -f testLD*


N_IND=24
N_SITES=10000


##### Genotypes
cat $SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testLD_T.geno.gz
zcat testLD_T.geno.gz | perl -an -e 'BEGIN{srand(12345)} if($pos > 10000) {$pos=0; $cnt++}; $pos += int(rand()*1000+1); print $F[0]."_".($cnt+1)."\t".$pos."\n"' > testLD.pos
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_T.geno.gz --pos testLD.pos --max_kb_dist 20 --min_maf 0.05 --extend_out                               | sort -k 1,1V -k 2,2V > testLD_T.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_T.geno.gz --pos testLD.pos --max_kb_dist 20 --min_maf 0.05 --extend_out --rnd_sample 0.5 --seed 12345 | sort -k 1,1V -k 2,2V > testLD_Tr.ld

##### Genotypes' likelihood and posterior probabilities
# Binary, normal-scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 3 -out testLD_3
gunzip testLD_3.glf.gz
rm testLD_3.glf.pos.gz
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_3.glf   --log_scale --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out                                              | sort -k 1,1V -k 2,2V > testLD_3.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_3.glf   --log_scale --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out --call_geno                                  | sort -k 1,1V -k 2,2V > testLD_3-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_3.glf   --log_scale --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1V -k 2,2V > testLD_3-CGf.ld

# Text, normal scale
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fas.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out testLD_2
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_2.beagle.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out                                              | sort -k 1,1V -k 2,2V > testLD_2.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_2.beagle.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out --call_geno                                  | sort -k 1,1V -k 2,2V > testLD_2-CG.ld
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_2.beagle.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out --call_geno --N_thresh 0.3 --call_thresh 0.9 | sort -k 1,1V -k 2,2V > testLD_2-CGf.ld
# Ignoring missing data
../ngsLD --n_threads 10 --verbose 1 --n_ind $N_IND --n_sites $N_SITES --geno testLD_2.beagle.gz --probs --pos testLD.pos --max_kb_dist 10 --min_maf 0.05 --extend_out --ignore_miss_data | sort -k 1,1V -k 2,2V > testLD_2.no_miss.ld
gunzip testLD_2.beagle.gz

# LD prunning
../scripts/prune_graph.pl --in_file testLD_2.ld --max_kb_dist 5 --field_weight 7 --min_weight 0.5 --print_excl testLD_pruned.id | sort -k 1,1V > testLD_unlinked.id



##### Check MD5
rm -f *.arg
TMP=`mktemp --suffix .ngsLD`
md5sum testLD* | fgrep -v '.gz' | sort -k 2,2 > $TMP
if diff $TMP test.md5 > /dev/null
then
    echo "ngsLD: All tests OK!"
else
    echo "ngsLD: test(s) failed!"
fi
