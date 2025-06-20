find ./data/normal_data -name "*.tsv" -exec mv {} ./data/normal_data \;
find ./data/normal_data -name "*.parcel" -type f -delete
find ./data/normal_data -name "annotations.txt" -type f -delete
find ./data/normal_data -type d -empty -delete
find ./data/tumor_data -name "*.tsv" -exec mv {} ./data/tumor_data \;
find ./data/tumor_data -name "*.parcel" -type f -delete
find ./data/tumor_data -name "annotations.txt" -type f -delete
find ./data/tumor_data -type d -empty -delete
