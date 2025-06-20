mkdir ./data/tumor_data
./gdc-client download -m ./data/manifests/tumor_gdc_manifest.txt -d ./data/tumor_data > ./results/tumor_download.log 2>&1
mkdir ../data/normal_data
./gdc-client download -m ./data/manifests/normal_gdc_manifest.txt -d ./data/normal_data > ./results/normal_download.log 2>&1
