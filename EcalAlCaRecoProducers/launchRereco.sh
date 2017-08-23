#!/bin/bash
trap "kill 0" SIGINT
json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt
jsonName=271036_275783-Prompt

json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt
jsonName=271036_276811-ICHEP

json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-279588_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt
jsonName=271036_279588-Prompt

json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
jsonName=271036_284044-23Sep2016

json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-299042_13TeV_PromptReco_Collisions17_JSON.txt
jsonName=294927_299042-Prompt

json25ns=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt
jsonName=294927-300575_Prompt_v1

##
PERIOD=LEGACY2016
PERIOD=RUN2017
#
#tags=( config/reRecoTags/Cal_Mar2017_ICcomb_v5.py)
tags=( config/reRecoTags/Cal_Aug2017_ref.py config/reRecoTags/Cal_Aug2017_Ped_v1.py config/reRecoTags/Cal_Aug2017_Ped_v2.py )
for tagfile in ${tags[@]}
do
	echo
#	./scripts/removeRereco.sh -t $tagfile -f alcarereco_datasets.dat
#	./scripts/removeRereco.sh -t $tagfile -f ntuple_datasets.dat --json_name=$jsonName
#	continue

	for CHECK in  '' --check
	do
		case $tagfile in 
			*/Cal_*_ref*.py)
				#./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json25ns --json_name="noJSON" ${CHECK} --alcarerecoOnly  --singleEle --weightsReco
		#		./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json25ns --json_name="noJSON" ${CHECK} --alcarerecoOnly  --singleEle
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json25ns --json_name="noJSON" ${CHECK} --alcarerecoOnly 
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json25ns --json_name="noJSON" ${CHECK} --alcarerecoOnly  --weightsReco
				;;
			*)
				echo 
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json25ns --json_name=$jsonName ${CHECK} --alcarerecoOnly
				;;
		esac
	done
	for CHECK in --check
	do
		case $tagfile in 
			*/Cal_*_ref*.py)
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json25ns --json_name=$jsonName --ntupleOnly  $CHECK --singleEle --weightsReco
		#		./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json25ns --json_name=$jsonName --ntupleOnly  $CHECK --singleEle
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json25ns --json_name=$jsonName --ntupleOnly  $CHECK 
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json25ns --json_name=$jsonName --ntupleOnly  $CHECK --weightsReco #| grep 'root;//' |sort |uniq > tmp/`basename $tagfile .py`.dat
				;;
			*)
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json25ns --json_name=$jsonName --ntupleOnly  $CHECK #| grep 'root;//' |sort |uniq > tmp/`basename $tagfile .py`.dat
				;;
		esac
	done

done
exit 0
