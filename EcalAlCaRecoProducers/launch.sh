#!/bin/bash
jsonDCS=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/DCSOnly/json_DCSONLY.txt
jsonDCS=/eos/project/c/cms-ecal-calibration/data/json/300576-301532_DCSonly.json

#json=$jsonDCS
jsonRereco=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
#Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt
jsonPrompt=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt

jsonNamePrompt=294927-300575_Prompt_v1
jsonNameRereco=271036-284044_23Sep2016_v1
jsonNameDCS=300576-301532_DCSonly

CHECK=--check
#CHECK=--createOnly
#CHECK=--submitOnly

#where=remoteGlidein
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

IFS=$'\n'
datasetsData=(`./scripts/parseDatasetFile.sh $fileList | grep VALID | sed 's|$|,|' | grep "${PERIOD}," | grep -v SIM`)
datasetsMC=(`./scripts/parseDatasetFile.sh $fileList | grep VALID | sed 's|$|,|' | grep "${PERIOD}," | grep SIM`)
# set IFS to newline in order to divide using new line the datasets


extraName=19Jul2017
extraName=1d296a55622a0222533deaa7b763d6d834e1bcec
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
			json=$jsonPrompt
			jsonName=$jsonNamePrompt
			./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_Prompt} -s noSkim --scheduler=${scheduler}   --json=$json --json_name=$jsonName --extraName=${extraName} ${CHECK} $dataset
			case $datasetName in 
				*Run2017C-noSkim-Prompt-v2*|*Run2017C-noSkim-Prompt-v3*)
					./scripts/prodNtuples.sh --type=MINIAOD -t ${tag_Prompt} -s noSkim --scheduler=${scheduler}   --json=$jsonDCS --json_name=$jsonNameDCS --extraName=${extraName} ${CHECK} $dataset
					;;
			esac
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

