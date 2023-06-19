#!/bin/bash

## Step 1: set up ENV variable and activate miniconda from CSD Python API
## 
echo "[*] set up CCD mapping enviroment"
## change the path below to py-rcsb_ccmodels folder
source /path_rcsb_ccmodels_installation/py-rcsb_ccmodels/ccdc-api-env.sh

## Step 2: run CCD external mapping
## 
## set CLI command path
cli_command="/path_CSD_installation/Python_API_2022/miniconda/bin/cc_models_cli"
if [ ! -f ${cli_command} ]
then
    echo "[*] cli command does not exist, exit"
    exit 1
else
    echo "[*] run CLI command of ${cli_command}"
fi

## set final data path, the data path should be permanent
data_folder="/path_data_to_save_to"
if [ ! -d ${data_folder} ]
then
    echo "[*] destination data path does not exist, exit"
    exit 1
else
    echo "[*] final mapping file will be copied to ${data_folder}"
fi

## set run path, the run path should be at local disk, or else the mapping runs slow
run_folder="/path_to_run"
if [ ! -d ${run_folder} ]
then
    echo "[*] run path does not exist, exit"
    exit 1
else
    cd ${run_folder}
    echo "[*] run CCD mapping at ${run_folder}"
    if [ ! -d CACHE ]
    then
	mkdir CACHE
	echo "[*] create CACHE folder"
    else
	echo "[*] use existing CACHE folder, data will be overwritten"
    fi
fi

echo "[**] generarate search targets for both CSD and COD searches"
${cli_command} --generate --num_proc 4 --cache_path ${run_folder}/CACHE 

echo "[**] run CSD-CCD mapping"
echo "[***] CSD search"
${cli_command} --search_ccdc --num_proc 4 --cache_path ${run_folder}/CACHE 

echo "[***] CSD model build on matched target"
${cli_command} --build_ccdc --num_proc 4 --cache_path ${run_folder}/CACHE --build_align_type graph-relaxed-stereo-sdeq 

echo "[**] run COD-CCD mapping"
echo "[***] COD search"
${cli_command} --search_cod --num_proc 4 --chunk_size 10 --cache_path ${run_folder}/CACHE 

echo "[***] COD sdf files fetch"
${cli_command} --fetch_cod --num_proc 4 --chunk_size 10 --cache_path ${run_folder}/CACHE 

echo "[***] COD model build on matched target"
${cli_command} --build_cod --num_proc 4 --chunk_size 20 --cache_path ${run_folder}/CACHE --build_align_type graph-relaxed-stereo-sdeq --build_cod_timeout 120.0 

echo "[**] assembly both/either CSD and COD matches"
${cli_command} --assemble  --cache_path ${run_folder}/CACHE 


## Step 3: verify and copy file to permanent location, and send email notification
## 
result_folder=${run_folder}/CACHE/cc-model-files
file_date=$(date +'%Y-%m-%d')
filename="chem_comp_models-${file_date}.cif"
filepath=${result_folder}/${filename}
echo "[*] looking for result file at ${filepath}"

## Update email below to yours
if [ -f ${filepath} ]
then
    cp ${filepath} ${data_folder}
    echo "[*] copy result file to ${data_folder}/${filename}"
    echo -e "Subject:Succeed to generate CCD mapping of ${data_folder}/${filename}" | /usr/sbin/sendmail your_email@domain > mail.log
else
    echo "[*] error, no result file found at ${result_folder}"
    echo -e "Subject:Fail to generate CCD mapping at ${result_folder}" | /usr/sbin/sendmail your_email@domain > mail.log
fi
