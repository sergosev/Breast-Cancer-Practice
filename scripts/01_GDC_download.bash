T_DIR="./data/tumor_data"

if [ ! -d "$T_DIR" ]; then
	echo "🔧 Creating tumor_data directory in ./data"
	mkdir ./data/tumor_data
fi

echo "📥 Downloading tumor data..."
./gdc-client download -m ./data/manifests/tumor_gdc_manifest.txt -d ./data/tumor_data > ./results/tumor_download.log


N_DIR="./data/normal_data"

if [ ! -d "$N_DIR" ]; then
	echo "🔧 Creating normal_data directory in ./data"
	mkdir ./data/normal_data
fi

echo "📥 Downloading normal data..."
./gdc-client download -m ./data/manifests/normal_gdc_manifest.txt -d ./data/normal_data > ./results/normal_download.log
