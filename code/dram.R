module load conda/25.1.1
source activate /projects/ehi/data/0_Environments/conda/DRAM
for i in *.fa.gz; do DRAM.py annotate -i $i -o ${i/.fa.gz/_annotate} --threads 8 --min_contig_size 1500; done
