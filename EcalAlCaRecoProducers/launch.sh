#!/bin/bash

#CHECK=--check
#CHECK=--createOnly
#CHECK=--submitOnly

jsonDCS=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/DCSOnly/json_DCSONLY.txt
jsonNameDCS=DCSonly

jsonRereco=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
jsonNameRereco=271036-284044_23Sep2016_v1

jsonPrompt=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON.txt
jsonNamePrompt=294927-301141_Prompt_v1

scheduler=caf
tag_Prompt=config/reRecoTags/92X_dataRun2_Prompt_v8.py
tag_Rereco=config/reRecoTags/80X_dataRun2_2016LegacyRepro_v3.py
tag_Moriond=config/reRecoTags/80X_dataRun2_2016SeptRepro_v7.py
tag_PromptH=config/reRecoTags/80X_dataRun2_Prompt_v16.py
tag_MC=config/reRecoTags/80X_mcRun2_asymptotic_2016_TrancheIV_v7.py

fileList=alcareco_datasets.dat
PERIOD=RUN2017
#PERIOD=LEGACY2016
#PERIOD=MORIOND2017
#PERIOD=MORIOND17 # MC


if  git diff-index --quiet HEAD; then
	GITCOMMIT=`git rev-parse HEAD`
	if [ "`git rev-parse HEAD`" != "`git rev-parse origin/master`" ];then
		echo "[ERROR] You are not allowed to make any production if all commits are propagated to the master branch of the repository" >> /dev/stderr
		exit 2
	fi
else
	echo "You have uncommitted changes, please commit everything before making a production" 
	exit 1
fi

extraName=$GITCOMMIT

# set IFS to newline in order to divide using new line the datasets
IFS=$'\n'
datasetsData=(`./scripts/parseDatasetFile.sh $fileList | grep VALID | sed 's|$|,|' | grep "${PERIOD}," | grep -v SIM`)
datasetsMC=(`./scripts/parseDatasetFile.sh $fileList | grep VALID | sed 's|$|,|' | grep "${PERIOD}," | grep SIM`)

for dataset in ${datasetsMC[@]} ${datasetsData[@]} #
do
	datasetName=`echo $dataset | awk '{print $6}'`
#	echo $datasetName
#	echo $dataset
	case $datasetName in
		
		*H-03Feb*)
			json=$jsonRereco
			jsonName=$jsonNameRereco
			#extraName=regressionMoriond17v2
##			./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_PromptH} -s noSkim --scheduler=${scheduler}   --json=$json --json_name=$jsonName  --extraName=${extraName} ${CHECK} $dataset
			;;
		
		*03Feb*)
			json=$jsonRereco
			jsonName=$jsonNameRereco
			#extraName=regressionMoriond17v2
##			./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_Moriond} -s noSkim --scheduler=${scheduler}   --json=$json --json_name=$jsonName  --extraName=${extraName} ${CHECK} $dataset
			;;
		*18Apr2017*)
			json=$jsonRereco
			jsonName=$jsonNameRereco
##			./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_Rereco} -s noSkim --scheduler=${scheduler}   --json=$json --json_name=$jsonName --extraName=${extraName} ${CHECK} $dataset
			;;
		*miniAODv2) #MC
			./scripts/prodNtuples.sh --type=MINIAOD --isMC -t ${tag_MC} -s noSkim --scheduler=${scheduler}    --extraName=${extraName} ${CHECK} $dataset
			;;
		*miniAOD) #data prompt 2017
			case $PERIOD in
				*_DCS)
					json=$jsonDCS
					jsonName=$jsonNameDCS
					;;
				*)
					json=$jsonPrompt
					jsonName=$jsonNamePrompt
					;;
			esac
			./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_Prompt} -s noSkim --scheduler=${scheduler}   --json=$json --json_name=$jsonName --extraName=${extraName} ${CHECK} $dataset
			;;
	esac

done
exit 0
for dataset in ${datasetsData[@]}
  do
	
	if [ "`echo $dataset | grep -c SingleElectron`" != "0" -a "`echo $dataset | grep -c ZElectron`" != "0" ];then continue; fi
	datasetName=`echo $dataset | awk '{print $4}'`
#	echo $dataset
#	continue
	case $datasetName in
		*PromptReco*)
			./scripts/prodNtuples.sh  --type=MINIAOD -t ${tag_Prompt} -s noSkim  --json=${json} --json_name=${jsonName} --scheduler=${scheduler} --doEleIDTree --extraName=newNtuples ${CHECK} ${dataset} || exit 1
			;;
		*)
			./scripts/prodNtuples.sh  --type=MINIAOD -t ${tag_Rereco} -s noSkim  --json=${json} --json_name=${jsonName} --scheduler=${scheduler} --doEleIDTree --extraName=newNtuples ${CHECK} ${dataset} || exit 1
			;;
	esac
	
done
	
exit 0
while [ "1" == "1" ];do
sleep 15m
done
exit 0

#./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep Run2016B | grep -v MINIAOD | head -1 | tail -1` --type ALCARECO -t ${tag_Data} -s ZSkim  --json=${json} --json_name=${jsonName} --scheduler=${scheduler} #${CHECK}
#./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep Run2016B | grep -v MINIAOD | head -2 | tail -1` --type ALCARECO -t ${tag_Data} -s ZSkim  --json=${json} --json_name=${jsonName} --scheduler=${scheduler} #${CHECK}
#./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep Run2016C | grep -v  MINIAOD | head -2 | tail -1` --type ALCARECO -t ${tag_Data} -s ZSkim --json=${json} --json_name=${jsonName} --scheduler=${scheduler} #${CHECK}

#



# while [ ! -e "prod_alcareco/DYJets_madgraph-RunIISpring16/allRange/res/finished" ] || [ ! -e "prod_alcareco/DYJets_amcatnlo-RunIISpring16/allRange/res/finished" ]
# do
# sleep 5m	
# ./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep RunIISpring16 | grep -v MINIAODSIM | head -1` --isMC --type ALCARECO -t ${tag_MC} -s ZSkim  --scheduler=${scheduler} ${CHECK}
	
# ./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep RunIISpring16 | grep -v MINIAODSIM | head -2` --isMC --type ALCARECO -t ${tag_MC} -s ZSkim  --scheduler=${scheduler} ${CHECK}

# done


# while [ "1" == "1" ];
# do
#    sleep 5m
# ./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep RunIISpring16 | grep -v MINIAODSIM | head -1` --isMC --type ALCARECO -t ${tag_MC} -s ZSkim  --scheduler=${scheduler} ${CHECK}
# ./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep RunIISpring16 | grep -v MINIAODSIM | head -2` --isMC --type ALCARECO -t ${tag_MC} -s ZSkim  --scheduler=${scheduler} ${CHECK}
# done




# #./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep RunIISpring16-miniAODv1 | head -1` --isMC --type MINIAOD -t ${tag_MC} -s noSkim  --scheduler=${where} ${CHECK}

# #if [[ $CHECK == "--createOnly" ]]; then
# #    crab -c prod_ntuples/MINIAODNTUPLE/80X_mcRun2_asymptotic_2016_v3/DYJets_amcatnlo-RunIISpring16-miniAODv1/allRange/RunII-2016/ -submit 1
# #fi

