library="$1"
ref="$2"

cd && cd artic/rampart

source activate artic-rampart

yarn --offline

npm run build

node rampart.js --basecalledDir /home/artic/Desktop/LIBRARIES/"$library"/BASECALLED_FILES/ --demuxedDir /home/artic/Desktop/RAMPART_OUTPUT/ --referencePanelPath /home/artic/Desktop/RAMPART_CONFIG_AND_REFERENCE_FILES/arbovirus_reference_genomes.fasta --referenceConfigPath /home/artic/Desktop/RAMPART_CONFIG_AND_REFERENCE_FILES/CHIKV_AsianECSA.json

