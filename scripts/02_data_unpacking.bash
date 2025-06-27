echo "🔄 Unpacking normal data"
echo "	✅ Putting all gene counts tables in ./data/normal_data"
find ./data/normal_data -name "*.tsv" -exec mv {} ./data/normal_data \;

echo "	✅ Deleting all .parcel files"
find ./data/normal_data -name "*.parcel" -type f -delete

echo "	✅ Deleting annotations"
find ./data/normal_data -name "annotations.txt" -type f -delete

echo "	✅ Deleting empty folders"
find ./data/normal_data -type d -empty -delete

echo "🔄 Unpacking tumor data"
echo "	✅ Putting all gene counts tables in ./data/normal_data"
find ./data/tumor_data -name "*.tsv" -exec mv {} ./data/tumor_data \;

echo "	✅ Deleting all .parcel files"
find ./data/tumor_data -name "*.parcel" -type f -delete

echo "	✅ Deleting annotations"
find ./data/tumor_data -name "annotations.txt" -type f -delete

echo "	✅ Deleting empty folders"
find ./data/tumor_data -type d -empty -delete
