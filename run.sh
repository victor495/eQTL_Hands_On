##Víctor Lorenzana Culebradas 
##eQTL Hands-On
cd eQTL
##task1
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g
##task2:
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt 
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz
tabix -p vcf input/processed/genotypes.chr22.vcf.gz
##Q1:
##Las opciones de BCFtools:
##          -v: comma-separated list of variant types to select.
##          -m: print sites with at least INT alleles listed in REF and ALT columns
## 	    -M: print sites with at most INT alleles listed in REF and ALT columns
##	    -q: minimum allele frequency (INFO/AC / INFO/AN) of sites to be printed.
##	    -S: File of sample names to include or exclude if prefixed with "^". 
##          -O: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v).
##          -o: FILE: output file name.
##	    -d: If a record is present multiple times, output only the first instance
##Q2:con esta funcion descomprimimos, eliminamos los comentarios y vemos el numero de lineas que corresponde con el numero de genotipos
zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l
##74656 genotipos
##Q3: obtengo un "text file" de los archivos VCF usando bcftools stats i luego miro el contenido usando la función less
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > tmp/stats.before
less tmp/stats.before 
##2504
bcftools stats input/processed/genotypes.chr22.vcf.gz > tmp/stats.after
less tmp/stats.after
##445
##task3:
##Q1: GENCODE v12
##Q2: GRCh37 (mirandolo desde la pagina de genecode)
##Q3: 19948
##Q4: PATH=$PATH:~/eQTL/bin
release=12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed
##Q5: con la funcion awk
##Q6: las coordenadas BED empiezan en 0 mientras que GTF empieza en 1, por lo tando el TSS 10 en BED corresponde a TSS 9 en GTF.
##Q7: Porque no se puede leer i escribir el mismo archivo a la vez
##task4:
join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed
##Q1: Protein-coding genes
##Q2: si les dades no estan normalitzades no es poden fer estudis estadistics
##Q3-4: realitzaria diferents plots per veure la distribucio de les dades
normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
bgzip tmp/genes.chr22.norm.bed 
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed
##trask5:

check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf
##Q1:Els plots de les dades normalizades mostran una distribució normal a diferencia de les dades previes, es a dir, la normalitzacio es correcte.
##task6:
head -1 input/processed/genes.chr22.norm.bed | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
##Q1:Characteristics and factor values of the strain and population  
##task7:
##Q1:--vcf i --bed demanen la inut data, --scale i --center delimiten les dimesions del grafic. --maf es la frequencia minima per alel.
##Q2: conte coordenades individuals dels principals components i el percentatge de variacio per a cada component
QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes
pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf
##Q3:En el PCA de la expressio el punts estan mes o meny generalitzats en tot el grafic. En el pca del genotype en canvi es diferencien 2 grups diferents seguint l'eix de les x
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf
##Q4:population, laboratory, sexe
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf
##Q5: laboratory
##task8:
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
gzip input/processed/covariates.tsv
##Q1: Tenen més variabilitat de mitja
##task9:
QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 0.01 --out result/nominals.txt
##Q1:
##pvalor:
root@2ea1a1f74577:/home/eqtl/eQTL# wc -l <(awk '{print $12}' result/nominals.txt | uniq)
##31725 /dev/fd/63
##beta size:
root@2ea1a1f74577:/home/eqtl/eQTL# wc -l <(awk '{print $13}' result/nominals.txt | uniq)
##31770 /dev/fd/63
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf
##Q2: la majoria dels valors de p corresponen a menys del 0.002
plink --ld <snp_id1> <snp_id2> --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2
##Q3:snps rs36084991 rs5746938. 

#Haplotype     Frequency    Expectation under LE
#   ---------     ---------    --------------------
#         AGT      0.246067                0.060549
#         TGT     -0                       0.185518
#          AG     -0                       0.185518
#          TG      0.753933                0.568414

#   In phase alleles are AGT/TG

##task10:He fet un script.sh(amb la seq. que dones) i le executat amb la terminal
./permutation.sh
## realitzo un altre script amb la sequencia de R anomenat Rchek.R
Rscript Rchek.R
##task11:
mtc.R -n result/nominals.txt -p result/permutation.txt --method "bonferroni" --alpha 0.05 --out tmp/bonferroni.txt
mtc.R -n result/nominals.txt -p result/permutation.txt --method "fdr" --alpha 0.05 --out tmp/fdr.txt
mtc.R -n result/nominals.txt -p result/permutation.txt --method "perm-fdr" --alpha 0.05 --out result/eqtls.tsv
#Q1:
wc -l result/eqtls.tsv
wc -l result/nominals.txt
##task12:
eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose
##task13:
rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed
for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
  bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed
sed -i "s/^chr//" input/processed/ERB.collapsed.bed
for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do 
  QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt
plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf
##Q1:H3K36me3, H3K79me2, H3K4me1, H4K20me1, etc. Son histones.
##Q2:Els valors per sota 1 de odds ratio son valors no esperats.
##task14:
sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv
## Com que no funciona el ENSEMBL amb tots els SNPs, agago al atzar fins a trobar alguna variant amb high impact.
##Q1:intron_variant: 47%
##downstream_gene_variant: 14%
##upstream_gene_variant: 15%
##non_coding_transcript_variant: 11%
##3_prime_UTR_variant: 1%
##NMD_transcript_variant: 6%
##regulatory_region_variant: 2%
##non_coding_transcript_exon_variant: 3%
##synonymous_variant: 2%
##intergenic_variant: 2%
##Others
##Q2: 23 high impacts, consequence: stop_gained, frameshift_variants,splice_acceptor_variant, non_coding_transcript_variant
##Q3:6
##task15:
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt
##Q1: response to lipopolysaccharide: Es qualsevol proces que canvii de estat o activitat degut a un lipopolisacarid. Es un component principal de la paret celular de bacteries gram negatives
##response to molecule of bacterial origin: Es qualsevol proces que canvii degut a una molecula de origen bacteria. 
##task16:
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt
wc -l <(awk '$20 > "0.9"' result/rtc.txt )
##Q1:38
awk '$1==$2' result/rtc.txt
##Q2:rs909685, es un snp relacionat amb la artritis reumatoide.
##Q3:intron_variant: 75%
##NMD_transcript_variant: 13%
##non_coding_transcript_variant: 13%
##task17:
gene=
compZscore.R --gene ENSG00000237438.1 --nominal result/nominals.txt -k 50 --output tmp/ENSG00000237438.1.rs_z
plink --r square --snps $(cut -f1 tmp/ENSG00000237438.1.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ENSG00000237438.1
CAVIAR -z tmp/ENSG00000237438.1.rs_z -l tmp/ENSG00000237438.1.ld -o result/ENSG00000237438.1
##Q1: 2 variants, 0.97869 and 0.934763
##Q2: 
grep "rs2845395" result/nominal.txt
##pvalue=7.76821e-09
grep "rs2845402" result/nominals.txt
##pvalue=4.92356e-08
##Els pvalors comparant amb altres variants son molt més petits.
grep "rs3016112" result/nominals.txt
##pvalue=0.00276116
##Q3:intron_variant: 38%
##non_coding_transcript_variant: 38%
##downstream_gene_variant: 14%
##upstream_gene_variant: 10%
##task18:
gene= 
cat <(echo "MarkerName P.value") <(grep rs909685 result/nominals.txt | cut -d " " -f8,12) > tmp/metal.rs909685
##he guardat el plot el result/plots
##task19:
##Fet.
