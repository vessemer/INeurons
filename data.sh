DATASET_DIR='../datasets'


# Lin HC He Z Ebert S Schörnig M Santel M Weigert A Hevers W Nadif Kasri N 
# Taverna E Camp JG Treutlein B (2020) Mendeley Data scRNAseq dataset.
# https://doi.org/10.17632/y3s4hnyvg6

cd $DATASET_DIR
wget 'https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/y3s4hnyvg6-1.zip' --output-document 'MendeleyData.zip'
mkdir MENDELEY; mv MendeleyData.zip MENDELEY; cd MENDELEY; unzip MendeleyData.zip; cd ..
for p in `ls MENDELEY/*/*.gz`; do echo $p; gzip -d $p; done;


# Ju X Schörnig M Ebert S Treutlein B Taverna E 
# (2020) ArrayExpress ID E-MTAB-9233. scRNAseq dataset.
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9233/

wget 'https://www.ebi.ac.uk/biostudies/files/E-MTAB-9233/zip' --output-document 'E-MTAB.zip'
unzip E-MTAB.zip

