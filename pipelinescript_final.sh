#!/bin/bash

echo "OMG, script is starting!"
echo "csv contig file is being made"
#The config file for the germline-DNA pipeline will be created in this script.

echo "makecsvfile.sh - hello"
read -p "Enter the first FASTQ.gz file here: " FastqR1file
read -p "(Optional) Enter the second FASTQ.gz file here: " FastqR2file
read -p "(Optional) Enter the VCF file u want to compare with later: " VCFfile2
read -p "Enter a sample name. Please do not use a sample name that has already been used: " Sample
read -p "Enter a name for the config.csv file here: " config
read -p "Enter an output directory here: " outputdirectory
echo $FastqR1file, $FastqR2file, $Sample, $config, $outputdirectory

sudo mkdir $outputdirectory
#Checking if there are R1 and R2 files. If there is only a R1 file the config file for only R1 will be made with the MD5 code..
#If both R1 and R2 files are in the folder, the config file for both files will be made with the MD5 codes..
#If there is no fq file, the script will stop.
if [[ -z $FastqR1file ]]
then
   echo "No Fastq R1 file given or input name is incorrect!"
   exit
elif [[ ! -z $FastqR1file ]] && [[ -z $FastqR2file ]]
then
   echo "No Fastq R2 file given, pipeline will run with only R1."
   md5R1=($(md5sum $FastqR1file))
   sudo touch files/$config
   echo '"sample","library","readgroup","R1","R1_md5"
"'$Sample'","lib1","rg2","'$FastqR1file'","'$md5R1'"' | sudo tee files/$config

elif [[ ! -z $FastqR2file ]] && [[ ! -z $FastqR1file ]]
then
   echo $FastqR1file
   echo $FastqR2file
   md5R1=($(md5sum $FastqR1file))
   md5R2=($(md5sum $FastqR2file))
   sudo touch files/$config
   echo '"sample","library","readgroup","R1","R1_md5","R2","R2_md5"
"'$Sample'","lib1","rg2","'$FastqR1file'","'$md5R1'","'$FastqR2file'","'$md5R2'"' | sudo tee files/$config 

fi

echo $md5R1
echo $md5R2
sudo chmod +x files/$config
echo "config file has been made"

#The input.json file for the germline-DNA pipeline will be made here. Change the fasta and index files to the correct files.
sudo touch files/input-germline-$Sample.json
echo '{
  "Somatic.bwaIndex": {
    "fastaFile": "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta",
    "indexFiles": [
      "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.sa",
      "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.amb",
      "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.ann",
      "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.bwt",
      "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.pac"
    ]
  },
  "Somatic.dbsnpVCF": "DBSNP/homo_sapiens_somatic.vcf.gz",
  "Somatic.dbsnpVCFIndex": "DBSNP/homo_sapiens_somatic.vcf.gz.tbi",
  "Somatic.sampleConfigFile": "files/'$config'",
  "Somatic.referenceFasta": "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta",
  "Somatic.referenceFastaFai": "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.fai",
  "Somatic.referenceFastaDict": "files/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
  "Somatic.dockerImagesFile": "germline-DNA/dockerImages.yml"
}' | sudo tee files/input-germline-$Sample.json

#This is the command to run the BioWDL germline-DNA pipeline.

sudo java -jar germline-DNA/cromwell-54.jar run -o germline-DNA/cromwell_options_no_user.json -i files/input-germline-$Sample.json \
--imports germline-DNA/somatic_v4.0.0.zip germline-DNA/somatic_v4.0.0.wdl

echo "germline-DNA DONE"
Bamfile=""
Baifile=""

#Here the script checks for BAM and BAI files in all the output files the germline-DNA pipeline generates.

for file in outputfiles/samples/$Sample/*; do
    if [[ $file == *bqsr.bam ]]
    then
       Bamfile2="${file//\.bqsr/}"
       sudo mv "$file" "$Bamfile2";
       Bamfile="$Bamfile2"
       echo Bamfile
       echo $Bamfile
    elif [[ $file == *bqsr.bai ]]
    then 
       Baifile2="${file//\.bqsr/}"
       sudo mv "$file" "$Baifile2";
       Baifile="$Baifile2"
       echo Baifile
       echo $Baifile
    fi
done;

#Here the scripts makes the cromwell_options.json file so the outputdirectory is the given outputdirectory.

sudo touch gatk-variantcalling/cromwell_options_$Sample.json
echo '{
"final_workflow_outputs_dir": "'$outputdirectory'",
"use_relative_output_paths": true,
"default_runtime_attributes": {
  "docker_user": "$EUID"
  }
}' | sudo tee gatk-variantcalling/cromwell_options_$Sample.json

#Here the input JSON file for BioWDL GATK-variantcalling will be made. Change the fasta file and index files to the correct files.

sudo touch files/inputfile-GATK-$Sample.json
echo '{
  "MultisampleCalling.referenceFasta": "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta",
  "MultisampleCalling.referenceFastaFai": "files/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.fai",
  "MultisampleCalling.referenceFastaDict": "files/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
  "MultisampleCalling.bamFilesAndGenders": [
    {"file": "'$Bamfile'",
     "index":  "'$Baifile'",
     "gender":"NULL"}
   ]
}' | sudo tee files/inputfile-GATK-$Sample.json

sudo chmod -R 777 files/inputfile-GATK-$Sample.json
sudo chmod -R 777 gatk-variantcalling/cromwell_options_$Sample.json

#This is the command to run the BioWDL GATK-variantcalling pipeline:

sudo java -jar gatk-variantcalling/cromwell-54.jar run -o gatk-variantcalling/cromwell_options_$Sample.json -i files/inputfile-GATK-$Sample.json \
gatk-variantcalling/multisample-variantcalling.wdl

#Checks for the generated vcf.gz file and compares it with the VCF.gz file from the other party.

VCFfile=""
for file in $outputdirectory/*; do
    if [[ $file == *.vcf.gz ]]
    then
       VCFfile="$file"
    fi
done;

if [[ -z $VCFfile2 ]]
then
   bcftools isec $VCFfile $VCFfile2 -p isec
fi

#Looks for the vcf.gz file and gunzip the file to use it in the VEP tool.

for file in $outputdirectory/*; do
    if [[ $file == *.vcf.gz ]]
    then
       sudo gunzip $file
    fi
done;

VCFfileVEP=""
for file in $outputdirectory/*; do
    if [[ $file == *.vcf ]]
    then
       VCFfileVEP="$file"
    fi
done;

#The command below is to run the VEP tool. Change the version to the correct version.

./ensembl-vep/vep -i $VCFfileVEP --cache --cache_version 104 --offline --force_overwrite
