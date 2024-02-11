### Log in to server ###

#Create parkinsons directory in data directory and navigate to it
mkdir /data/parkinsons
cd /data/parkinsons

#Import data and demultiplex
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/parkinsons/parkinsons_manifest.txt \
  --output-path ./demux_seqs.qza

#Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux_seqs.qzv

### Open local command line (not in server) ###

#Navigate to Qiime directory in git repository folder
cd Desktop/MICB-475-Team-2/qiime_files

#Import server demux file into git repository
scp root@10.19.139.167:/data/parkinsons/demux_seqs.qzv .

### Visualize demux_seqs.qzv at view.qimme2.org ###

### Denoising and clustering while in /data/parkinsons directory ###

#Determine ASVs with DADA2
  qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 251 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

# Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

# Visualize ASVs stats
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/parkinsons/parkinsons_metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

### Taxonomic Analysis ###

qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

# Taxonomy barplots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /mnt/datasets/project_2/parkinsons/parkinsons_metadata.txt \
  --o-visualization taxa-bar-plots.qzv

### Open local command line (not in server) ###

#Navigate to Qiime directory in git repository folder
cd Desktop/MICB-475-Team-2/qiime_files

#Import server taxonomy files into git repository
scp root@10.19.139.167:/data/parkinsons/taxonomy.qzv .

scp root@10.19.139.167:/data/parkinsons/taxa-bar-plots.qzv .

### Visualize taxonomy files at view.qimme2.org ###

### Filtering while in /data/parkinsons directory ###
#Remove mitochondria and chloroplasts
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/parkinsons/parkinsons_metadata.txt

### Open local command line (not in server) ###

#Navigate to Qiime directory in git repository directory
cd Desktop/MICB-475-Team-2/qiime_files

#Import server filtered table into git repository
scp root@10.19.139.167:/data/parkinsons/table-no-mitochondria-no-chloroplast.qzv .

### While in /data/parkinsons directory ###

#Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

### Exporting files to git repository for R analysis ###

#Navigate to data directory
cd /data

#Make new export directory and navigate to it
mkdir parkinsons_export
cd parkinsons_export

#Export files
qiime tools export \
--input-path ../parkinsons/table.qza \
--output-path table_export 

qiime tools export \
--input-path ../parkinsons/taxonomy.qza \
--output-path taxonomy_export 

qiime tools export \
--input-path ../parkinsons/rooted-tree.qza \
--output-path rooted_tree_export

#Navigate to table_export and change feature-table.biom to feature-table.txt
cd table_export

biom convert \
-i feature-table.biom \
--to-tsv \
-o feature-table.txt

### Open local command line (not in server) to export parkinsons_export###

#Navigate to Qiime directory in git repository directory
cd Desktop/MICB-475-Team-2/qiime_files

#Import parkinsons_export to git repository
scp -r root@10.19.139.167:~/data/parkinsons_export .

#Navigate to parkinsons_export and import parkinsons_metadata.txt to git repository

cd Desktop/MICB-475-Team-2/qiime_files/parkinsons_export

scp root@10.19.139.167:/mnt/datasets/project_2/parkinsons/parkinsons_metadata.txt .
