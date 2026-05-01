# QIIME2 pipeline for microbiome analysis
# Author: Anna Kindberg
# Description: Processes raw sequences and computes diversity measures
# Environment: rachis-qiime2-2026.4     

# Import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path subset_manifest.txt \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Visualize demultiplexed sequence quality
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# Quality filtering and denoising
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 240 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 200 \
  --o-representative-sequences asv-seqs.qza \
  --o-table asv-table.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

# Summarize feature-table
qiime feature-table summarize \
  --i-table asv-table.qza \
  --m-metadata-file BNZ_Fire.16S.metadata.tsv \
  --o-summary asv-table.qzv \
  --o-sample-frequencies sample-frequencies.qza \
  --o-feature-frequencies asv-frequencies.qza

# Visualize denoised sequence quality
qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --m-metadata-file asv-frequencies.qza \
  --o-visualization asv-seqs.qzv

# Import SILVA data base for training the classifier
qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

# Reverse transcribe SILVA sequences
qiime rescript reverse-transcribe \
—I-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
—o-dna-sequences silva-138.2-ssu-nr99-dna-seqs.qza

# Train the Naive Bayes taxonomy classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138.2-ssu-nr99-dna-seqs.qza \
  --i-reference-taxonomy silva-138.2-ssu-nr99-tax.qza \
  --o-classifier silva-16S-rRNA-classifier.qza

# Classify ASV sequences
qiime feature-classifier classify-sklearn \
  --i-classifier silva-16S-rRNA-classifier.qza \
  --i-reads asv-seqs.qza \
  --o-classification silva-taxonomy.qza

# Summarize and visualize taxonomy
qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --i-taxonomy silva-taxonomy.qza \
  --m-metadata-file asv-frequencies.qza \
  --o-visualization asv-seqs-tax-sum.qzv

# Make taxa barplots
qiime taxa barplot \
—I-table asv-table.qza\
—i-taxonomy silva-taxonomy.qza \
—-m-metadata-file /mnt/drownlab/projects/BNZ_data/metadata/BNZ_Fire.16S.metadata.tsv \
—o-visualization taxa-bar-plots.qzv

# Generate Diversity Measures at specified sampling-depth
qiime boots kmer-diversity \
  --i-table asv-table.qza \
  --i-sequences asv-seqs.qza \
  --m-metadata-file BNZ_Fire.16S.metadata.tsv \
  --p-sampling-depth 3000 \
  --p-n 10 \
  --p-replacement \
  --p-alpha-average-method median \
  --p-beta-average-method medoid \
  --output-dir boots-kmer-diversity



