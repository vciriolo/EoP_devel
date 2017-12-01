#!/bin/bash
trap "kill 0" SIGINT

json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt
jsonName=271036_276811-ICHEP

json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-279588_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt
jsonName=271036_279588-Prompt


json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON.txt
jsonName=294927-305364_Prompt_v1

json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
jsonName=271036_284044-23Sep2016

json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306126_13TeV_PromptReco_Collisions17_JSON.txt
jsonName=294927-306126_Prompt_v1

##
PERIOD=LEGACY2016
PERIOD=RUN2017F
#
tags=(  config/reRecoTags/Cal_Nov2017_{ref,IC}_v1.py  )
tags=(  config/reRecoTags/Cal_Oct2017_cand_v5.py )

if  git status --porcelain -uno | grep -v launch | grep -v ZFitter | grep -q -v _datasets  ; then
	echo "You have uncommitted changes, please commit everything before making a production" 
	exit 1
else
	GITCOMMIT=`git rev-parse HEAD`
	if [ "`git rev-parse HEAD`" != "`git rev-parse origin/master`" ];then
		echo "[ERROR] You are not allowed to make any production if all commits are propagated to the master branch of the repository" >> /dev/stderr
#		exit 2
	fi
fi


for tagfile in ${tags[@]}
do
	echo
#	./scripts/removeRereco.sh -t $tagfile -f alcarereco_datasets.dat
#	./scripts/removeRereco.sh -t $tagfile -f ntuple_datasets.dat --json_name=$jsonName
#	continue

	for CHECK in --createOnly #''  --check
	do
		case $tagfile in 
			*/Cal_*_cand*.py )
				#./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json --json_name="noJSON" ${CHECK} --alcarerecoOnly  --singleEle --weightsReco
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json --json_name="noJSON" ${CHECK} --alcarerecoOnly  --singleEle
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json --json_name="noJSON" ${CHECK} --alcarerecoOnly 
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json --json_name="noJSON" ${CHECK} --alcarerecoOnly  --weightsReco
				;;
			*)
				echo 
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile  --json=$json --json_name=$jsonName ${CHECK} --alcarerecoOnly 
				;;
		esac
	done
continue
	for CHECK in   --check
	do
		case $tagfile in 
			*/Cal_*ref*.py)
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json --json_name=$jsonName --ntupleOnly  $CHECK --singleEle --weightsReco
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json --json_name=$jsonName --ntupleOnly --singleEle  $CHECK 
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json --json_name=$jsonName --ntupleOnly  $CHECK 
#				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json --json_name=$jsonName --ntupleOnly  $CHECK --weightsReco #| grep 'root;//' |sort |uniq > tmp/`basename $tagfile .py`.dat
				;;
			*)
				./scripts/RerecoQuick.sh -p ${PERIOD} -t $tagfile --json=$json --json_name=$jsonName --ntupleOnly  $CHECK  #| grep 'root;//' |sort |uniq > tmp/`basename $tagfile .py`.dat
				;;
		esac
	done

done
exit 0
