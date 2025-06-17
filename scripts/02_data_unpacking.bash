find ./data/normal_data2 -name "*.tsv" -exec mv {} ./data/normal_data2 \;
find ./data/normal_data2 -name "*.parcel" -type f -delete
find ./data/normal_data2 -name "annotations.txt" -type f -delete
find ./data/normal_data2 -type d -empty -delete
find ./data/tumor_data2 -name "*.tsv" -exec mv {} ./data/tumor_data2 \;
find ./data/tumor_data2 -name "*.parcel" -type f -delete
find ./data/tumor_data2 -name "annotations.txt" -type f -delete
find ./data/tumor_data2 -type d -empty -delete
