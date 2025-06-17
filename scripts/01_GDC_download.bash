mkdir ./data/tumor_data2
./gdc-client download -m ./data/manifests/tumor_gdc_manifest2.txt -d ./data/tumor_data2 > ./results/tumor2_download.log 2>&1
mkdir ../data/normal_data2
./gdc-client download -m ./data/manifests/normal_gdc_manifest2.txt -d ./data/normal_data2 > ./results/normal2_download.log 2>&1
