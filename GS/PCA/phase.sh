for i in {1..22}
do
shapeit/bin/shapeit --input-bed chr_$i --input-map chr_$i'_map.txt' --output-max haplotypes/chr_$i'.haps' haplotypes/chr_$i'.sample' --thread 48 --duohmm -W 5
done