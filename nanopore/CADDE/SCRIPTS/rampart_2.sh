# version edited by Filipe M.

library="$1"		# Complete path for the  BASECALLED_FILES for the specific library
ref="$2"

cd && cd artic/rampart

source activate artic-rampart

#yarn --offline

#Npm run build

rampart --basecalledPath $library --referencesPath /home/artic/Desktop/RAMPART_CONFIG_AND_REFERENCE_FILES/arbovirus_reference_genomes.fasta
