#!/bin/bash
source $CMSSW_BASE/src/Calibration/EcalAlCaRecoProducers/scripts/prodFunctions.sh

############################### OPTIONS
#------------------------------ default
SCHEDULER=caf
QUEUE=cmscaf1nd
STORAGE_ELEMENT=caf
isMC=0
UI_WORKING_DIR=prod_ntuples
USER_REMOTE_DIR_BASE=group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples
LUMIS_PER_JOBS=100
EVENTS_PER_JOB=300000

DOTREE=1
DOEXTRACALIBTREE=0
DOEXTRASTUDYTREE=0
DOELEIDTREE=1
CREATE=y
SUBMIT=y

SKIM=""
OUTFILES="ntuple.root"
if [ "$DOELEIDTREE" != "0" ]; then
	OUTFILES="${OUTFILES},eleIDTree.root"
fi

CRABVERSION=2
crab2File=tmp/ntuple.cfg
crab3File=tmp/ntuple.py

case ${CRABVERSION} in
	2) crabFile=${crab2File};;
	3) crabFile=${crab3File};;
	*) 
		echo "[`basename $0` ERROR] CRABVERSION = ${CRABVERSION} not valid" >> /dev/stderr
		exit 1
		;;
esac

JOBINDEX=""
ISPRIVATE=0
BUNCHSPACING=0

usage(){
    echo "`basename $0` {parseDatasetFile options} --type={type} -t tagFile [options]"
    echo "---------- provided by parseDatasetFile (all mandatory)"
    echo "    -r runRange"
    echo "    -d, --datasetpath path"
    echo "    -n, --datasetname name"
    echo "    --store dir"
    echo "    --remote_dir dir: origin files remote dir"

    echo "---------- provided by command-line (mandatory)"
    echo "    --type TYPE: one of the following"
    echo "           ALCARECOSIM: alcareco produced on MC"
    echo "           ALCARERECO: alcareco format after rereco on ALCARAW"
    echo "                       supports only CAF as scheduler and CRAB2!"
    echo "           MINIAOD: ntuple production from miniAOD"
	echo "     -t tagFile: config file with the global tag. It's one of the files in config/reRecoTags/"
    echo "    *** for MC ***"
    echo "    --isMC"
#    echo "    --isParticleGun: redundant, --skim=partGun is the same"
    echo "    *** for DATA ***"
    echo "    --json_name jsonName: additional name in the folder structure to keep track of the used json"
    echo "    --json jsonFile.root"
	echo "    --weightsReco: only if producing from ALCARERECO using weights for local reco"

    echo "---------- optional common"
    echo "    --doExtraCalibTree (needed for E/p calibration)"
    echo "    --doExtraStudyTree (for extra studies on single rechit quantities)"
    echo "    --doEleIDTree      (has shower shapes)"
    echo "    --noStandardTree"
    echo "    --createOnly"
    echo "    --submitOnly"
    echo "    --check"
    echo "    --isPrivate: it is a privately produced dataset (not central or prompt)"
	echo "    --useParent" 

    echo "----------"
    echo "    --tutorial: tutorial mode, produces only one sample in you user area"
    echo "    -H, --expertHelp"
}

expertUsage(){
    echo "------------------------------ Expert options"
    echo "    --ntuple_store arg: storage element for the ntuples (=${STORAGE_ELEMENT})"
    echo "    --ntuple_remote_dir dir (=${USER_REMOTE_DIR_BASE})"
    echo "    --scheduler caf,lsf,remoteGlidein (=${SCHEDULER})"
    echo "    -f, --filelist arg: produce ntuples from a list of files"
    echo "    -s, --skim arg: apply the skim (ZSkim, WSkim)"
    echo "    --ui_working_dir arg: crab task folder (=${UI_WORKING_DIR})"
    echo "    --extraName arg: additional name for folder structure (to make different versions) (='')"
    echo "    --file_per_job arg: number of files to process in 1 job (=1)"
    echo "    --develRelease: CRAB do not check if the CMSSW version is in production (only if you are sure what you are doing)"
}

    
#------------------------------ parsing

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o hHd:n:s:r:t:f: -l help,expertHelp,datasetpath:,datasetname:,skim:,runrange:,store:,remote_dir:,scheduler:,isMC,isParticleGun,ntuple_remote_dir:,json:,tag:,type:,json_name:,ui_working_dir:,extraName:,doExtraCalibTree,doExtraStudyTree,doEleIDTree,noStandardTree,createOnly,submitOnly,check,isPrivate,file_per_job:,develRelease,weightsReco,useParent -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

checkVerboseOption

set -- $options
#echo $options

while [ $# -gt 0 ]
do
    case $1 in
	-h|--help) usage; exit 0;;
	-H|--expertHelp)  expertUsage; exit 0;;
	-r|--runrange) RUNRANGE=$2; shift;;
	-d|--datasetpath) DATASETPATH=$2; shift ;;
	-n|--datasetname) DATASETNAME=$2; shift ;;
	--store) ORIGIN_STORAGE_ELEMENT=$2; 
	    case ${ORIGIN_STORAGE_ELEMENT} in
		caf) ;;
		caf.cern.ch) ;;
		*)
		    echo "[ERROR] Origin storage_element ${ORIGIN_STORAGE_ELEMENT} != caf not implemented yet" >> /dev/stderr
		    exit 1
		    ;;
	    esac
	    shift;;
	--remote_dir) ORIGIN_REMOTE_DIR_BASE=$2; shift;;

	--ntuple_remote_dir) USER_REMOTE_DIR_BASE=$2; echo ${USER_REMOTE_DIR_BASE}; shift;;

 	--type) TYPE=$2; shift;
	    case $TYPE in 
		alcareco | ALCARECO)
		    if [ "${isMC}" == "1" ]; then 
			TYPE=ALCARECOSIM; 
		    else
			TYPE=ALCARECO
		    fi
		    ;;
		ALCARECOSIM)
		    isMC=1
		    ;;
		MINIAOD| miniAOD)
				TYPE=MINIAODNTUPLE
				;;
		alcarereco | ALCARERECO)
		    TYPE=ALCARERECO
		    if [ "${isMC}" == "1" ]; then
				echo "[ERROR `basename $0`] Incompatible options: TYPE=${TYPE} and --isMC" >> /dev/stderr
				exit 1
		    fi
			CRABVERSION=2
		    ;;
		*)
		    echo "[OPTION ERROR] Type $TYPE not recognize" >> /dev/stderr
		    usage >> /dev/stderr
		    exit 1
		    ;;
	    esac
	    ;;
	--isMC) isMC=1;;
	--isParticleGun) isPARTICLEGUN="y"; SKIM=partGun;;
 	--json) JSONFILE=$2;  shift;;
	--json_name) JSONNAME=$2; shift;;

	-f|--filelist) FILELIST="$FILELIST $2"; echo ${FILELIST}; shift ;;
	-s|--skim) SKIM=$2 ; shift;;
	-t | --tag) 
			if [ -n "${TAG}" ];then
				if [ "`basename $2 .py`" != "${TAG}" ];then
					echo "[ERROR] -t options called twice with different tag!" >> /dev/stderr
					exit 1
				fi
			else
				TAGFILE=$2; 
				case $TAGFILE in
					config/reRecoTags/*) TAG=`basename ${TAGFILE} .py` ;;
					*) TAG=$TAGFILE; TAGFILE=config/reRecoTags/$TAG.py;;
				esac						
				if [ -n "${VERBOSE}" ];then	echo "[OPTION] GLOBALTAG:$TAGFILE"; fi
			fi
			shift 
			;;
	--ntuple_store) STORAGE_ELEMENT=$2; shift;;
	--ui_working_dir) UI_WORKING_DIR=$2; shift;;
	--scheduler) SCHEDULER=$2; shift;;
	#--puWeight) PUWEIGHTFILE=$2; shift;;
	--extraName) EXTRANAME=$2;shift;;
        #name of the output files is hardcoded in ZNtupleDumper
	--doExtraCalibTree)
			if [ -n "${VERBOSE}" ]; then
				echo "[OPTION] doExtraCalibTree"; 
			fi
			DOEXTRACALIBTREE=1; OUTFILES="${OUTFILES},extraCalibTree.root"
			;;
	--doExtraStudyTree)
			if [ -n "${VERBOSE}" ];then 
				echo "[OPTION] doExtraStudyTree"; 
			fi
			DOEXTRASTUDYTREE=1; OUTFILES="${OUTFILES},extraStudyTree.root"
			;;
	--doEleIDTree)
			DOELEIDTREE=1;OUTFILES="${OUTFILES},eleIDTree.root"
			;;
	--noStandardTree) DOTREE=0; OUTFILES=`echo ${OUTFILES} | sed 's|ntuple.root,||'`;;
	--createOnly) echo "[OPTION] createOnly"; unset SUBMIT;;
	--submitOnly) echo "[OPTION] submitOnly"; unset CREATE;;
	--check)      
			#echo "[OPTION] checking jobs"; 
			CHECK=y; EXTRAOPTION="--check"; unset CREATE; unset SUBMIT;;
	--isPrivate)      echo "[OPTION] private dataset"; ISPRIVATE=1;;
	--useParent) USEPARENT=y ;;
 	--file_per_job)
			#echo "[OPTION] file per job: $2"; 
			FILE_PER_JOB=$2; shift ;;
	--develRelease) echo "[OPTION] Request also CMSSW release not in production!"; DEVEL_RELEASE=y;;
	--weightsReco)    echo "[OPTION `basename $0`] using weights for local reco"; BUNCHSPACING=-1;;

	(--) shift; break;;
	(-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
	(*) break;;
    esac
    shift
done

if [ -z "$DATASETNAME" ];then 
    echo "[ERROR] DATASETNAME not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

if [ -z "$RUNRANGE" ];then 
    echo "[ERROR] RUNRANGE not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

if [ -z "$TYPE" ];then 
    echo "[ERROR] TYPE not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

if [ -z "$JSONFILE" -a "$isMC" != "1" ];then 
    echo "[ERROR] JSONFILE not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

if [ -z "$JSONNAME" -a "$isMC" != "1" ];then 
    echo "[ERROR] JSONNAME not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

if [ -z "${TAG}" ];then
	echo "[ERROR] Tag not defined" >> /dev/stderr
	usage >> /dev/stderr
	exit 1
fi

if [ "$DOEXTRACALIBTREE" != "0" -a "$TYPE" != "ALCARERECO" ]; then
	echo "[ERROR] cannot produce the extra calib tree other datasets than ALCARERECO because there are no uncalib rechits"
	exit 1
fi

if [ "$DOEXTRASTUDYTREE" != "0" -a "$TYPE" != "ALCARERECO" ]; then
	echo "[ERROR] cannot produce the extra calib tree other datasets than ALCARERECO because there are no uncalib rechits"
	exit 1
fi


case $SCHEDULER in
	CAF|caf)
		CRABVERSION=2
		;;
	LSF|lsf)
		CRABVERSION=2
		QUEUE=1nd
		;;
	*)
		CRABVERSION=3
		;;
esac

case $DATASETNAME in
#	*-50nsReco) BUNCHSPACING=50;;
#	*-25nsReco) BUNCHSPACING=25;;
	*-weightsReco) BUNCHSPACING=-1;;
esac

#Setting the ENERGY variable
setEnergy $DATASETPATH



if [ "${TYPE}" != "ALCARERECO" -a "${TYPE}" != "ALCARAW" ];then
    ORIGIN_REMOTE_DIR_BASE=`echo ${ORIGIN_REMOTE_DIR_BASE} | sed 's|alcaraw|alcareco|g'`
#elif [ "${TYPE}" != "alcarereco" ];then
#    echo ${ORIGIN_REMOTE_DIR_BASE}
#    ORIGIN_REMOTE_DIR_BASE=`echo ${ORIGIN_REMOTE_DIR_BASE} | sed 's|sandbox|alcareco|g'`
fi



#JOBINDEX=`cat crab_projects/crab_${DATASETNAME}/crab.log | grep -oh -m1 '[0-9]\{6\}_[0-9]\{6\}'` 
#echo " JOB INDEX $JOBINDEX"

#echo "[INFO] Temporarily set Crab to CRAB2 to allow CAF submission"
#echo "[INFO] You can set yourself back to CRAB3 by running iniCmsEnv.sh again"
#source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh


setStoragePath $STORAGE_ELEMENT $SCHEDULER 
###
###
# make the filelist before parsing the options and arguments
options="-d ${DATASETPATH} -n ${DATASETNAME} -r ${RUNRANGE} --remote_dir ${ORIGIN_REMOTE_DIR_BASE}"
case $TYPE in 
	ALCARERECO)
		options="$options -t ${TAGFILE}"; 
		setUserRemoteDirAlcarereco $ORIGIN_REMOTE_DIR_BASE
		ORIGIN_REMOTE_DIR=${USER_REMOTE_DIR}
		;;
	MINIAODNTUPLE)
		case $SKIM in 
			ZSkim)
				echo "[ERROR `basename $0`] currently no skim possible on miniaod: skim=$SKIM not possible" >> /dev/stderr
				exit 1
				;;
		esac
		#setUserRemoteDirMiniaod $ORIGIN_REMOTE_DIR_BASE
		#ORIGIN_REMOTE_DIR=${USER_REMOTE_DIR}
		ORIGIN_REMOTE_DIR=${ORIGIN_REMOTE_DIR_BASE}
		;;
	*)
#		TAG=""
		setUserRemoteDirAlcareco $ORIGIN_REMOTE_DIR_BASE
		ORIGIN_REMOTE_DIR=${USER_REMOTE_DIR}
		;;
esac



if [ -z "$USER_REMOTE_DIR_BASE" ];then 
    echo "[ERROR] USER_REMOTE_DIR_BASE not defined" >> /dev/stderr
    usage >> /dev/stderr
    exit 1
fi

UI_WORKING_DIR=$UI_WORKING_DIR/${TYPE}/${TAG}/${DATASETNAME}/${RUNRANGE}/${JSONNAME}/${EXTRANAME}

USER_REMOTE_DIR=$USER_REMOTE_DIR_BASE/${ENERGY}/${TYPE}
if [ -n "${TAG}" ];then USER_REMOTE_DIR=$USER_REMOTE_DIR/${TAG}; fi
USER_REMOTE_DIR=$USER_REMOTE_DIR/${DATASETNAME}/${RUNRANGE}
if [ -n "${JSONNAME}" ];then USER_REMOTE_DIR=$USER_REMOTE_DIR/${JSONNAME}; fi
if [ -n "${EXTRANAME}" ];then USER_REMOTE_DIR=$USER_REMOTE_DIR/${EXTRANAME}; fi

if [ -z "${CHECK}" ];then
	if [ "${TYPE}" == "ALCARERECO" ];then
		if [ "`cat ntuple_datasets.dat | grep ${DATASETNAME}  | grep ${JSONNAME} | grep $TAG$ | grep -c $RUNRANGE`" != "0" ];then
			echo "[WARNING] Ntuple for rereco $TAG already done for ${RUNRAGE} ${DATASETNAME}"

			for file in `eos.select ls -l $STORAGE_PATH/$USER_REMOTE_DIR/  | sed '/^d/ d' | awk '{print $9}'`
			do 
				echo "[FILE] $STORAGE_PATH/$USER_REMOTE_DIR/$file"
			done
			exit 0
		fi
#else
#    if [ "`cat ntuple_datasets.dat | grep -v ALCARERECO | grep ${DATASETNAME} | grep ${JSONNAME} | grep -c $RUNRANGE`" != "0" ];then
#	echo "[WARNING] Ntuple for $JSONNAME  already done for ${RUNRANGE} ${DATASETNAME}"
#	exit 0
#    fi
	fi
fi
USER_REMOTE_DIR=$USER_REMOTE_DIR/unmerged


#==============================



###############################


#------------------------------




#${ENERGY}/
#${DATASETNAME}/tmp-${DATASETNAME}-${RUNRANGE}
OUTFILES=`echo $OUTFILES | sed 's|^,||'`
# echo ${ORIGIN_REMOTE_DIR}
# echo
# echo ${USER_REMOTE_DIR}
# exit 0
if [ -n "${CREATE}" ];then

case ${ORIGIN_REMOTE_DIR} in
    database)
		FILELIST=""
		;;
    *)
	echo ${FILELIST}
	if [ -z "${FILELIST}" ];then
	    FILELIST=filelist/${TAG}/${DATASETNAME}-${RUNRANGE}.list
	fi

	if [ ! -e "${FILELIST}" ];then
	    #sample=tempFileList-${DATASETNAME}-${RUNRANGE}-${TAG}
			echo $USER
			echo "makefilelist.sh -f filelist/${TAG} `basename ${FILELIST} .list` $STORAGE_PATH/ /${ORIGIN_REMOTE_DIR}"
	    makefilelist.sh -f "filelist/${TAG}" `basename ${FILELIST} .list` $STORAGE_PATH/${ORIGIN_REMOTE_DIR}  || exit 1
	    #filelistDatasets.sh $options || exit 1
            # remove PUDumper files!
	    if [ -n "$TAG" ];then
			sed -i '/PUDumper/ d' filelist/*/*.list
			sed -i '/ntuple/ d'   filelist/*/*.list
	    else
			sed -i '/PUDumper/ d' filelist/*.list
			sed -i '/ntuple/ d' filelist/*.list
	    fi
	fi
	;;
esac

if [ -n "$FILELIST" ]; then
    nFiles=`cat $FILELIST | wc -l`
    if [ -n "$NJOBS" ];then
	
		let FILE_PER_JOB=$nFiles/$NJOBS
		if [ "`echo \"$nFiles%$NJOBS\" | bc`" != "0" ];then
			let FILE_PER_JOB=$FILE_PER_JOB+1
		fi
    elif [ -n "$FILE_PER_JOB" ];then
		NJOBS=`perl -w -e "use POSIX; print ceil($nFiles/${FILE_PER_JOB}), qq{\n}"`
		if [ "`echo \"${nFiles}%${FILE_PER_JOB}\" | bc -l`" != "0" ];then
			let NJOBS=$NJOBS+1
		fi
    else
		NJOBS=$nFiles
		FILE_PER_JOB=1
    fi
fi

if [ "$RUNRANGE" == "allRange" -o "`echo $RUNRANGE |grep -c -P '[0-9]+-[0-9]+'`" == "0" ];then
    unset RUNRANGE
fi

if [ ! -d "tmp" ];then mkdir tmp/; fi

#Write the crab config file
cat > ${crab2File} <<EOF
UI_WORKING_DIR=$UI_WORKING_DIR
pset=python/alcaSkimming.py
psetparams="type=${TYPE} doTree=${DOTREE} doExtraCalibTree=${DOEXTRACALIBTREE} doExtraStudyTree=${DOEXTRASTUDYTREE} doEleIDTree=${DOELEIDTREE} doTreeOnly=1  jsonFile=${JSONFILE} isCrab=1 skim=${SKIM} tagFile=${TAGFILE} isPrivate=$ISPRIVATE  bunchSpacing=${BUNCHSPACING}"

outFiles=${OUTFILES}
use_parent=${USEPARENT}
queue=${QUEUE}

user_remote_dir="$USER_REMOTE_DIR"
storage_path="$STORAGE_PATH"
EOF
case ${ORIGIN_REMOTE_DIR_BASE} in
        database)
		case $DATASETPATH in
			*USER)
				if [ -z ${DBS_URL} ];then
					DBS_URL=phys03
				fi
				;;
		esac
		if [ -n "${DBS_URL}" ];then
			echo "dbs_url=${DBS_URL}" >> ${crab2File}
		fi

		if [ "$isMC" != "1" ];then
        cat >> ${crab2File} <<EOF
dataset=${DATASETPATH}
datasetname=${DATASETNAME}
runrange=${RUNRANGE}
datasetsite=T2_CH_CERN
#dbs_url = phys03
EOF
		else
			cat >> ${crab2File} <<EOF
dataset=${DATASETPATH}
datasetname=${DATASETNAME}
runrange=allRange
datasetsite=T2_CH_CERN
EOF
		fi
        ;;
        *)
        cat >> ${crab2File} <<EOF
dataset=None
datasetname=${DATASETNAME}
runrange=${RUNRANGE}
EOF
	;;
esac


cat >> ${crab2File} <<EOF
EOF

crab3outfiles=`echo $OUTFILES | sed "s|^|\'|;s|$|'|;s|,|\',\'|"`
###################################crab 3 .py writing###############################
#I used the tool: crab2cfgTOcrab3py [crab2confgiName.cfg] [crab3configName.py]
##jsonFile=${JSONFILE}
cat > ${crab3File} <<EOF
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = '${UI_WORKING_DIR}'
config.General.requestName = '${TYPE}'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/alcaSkimming.py'
config.JobType.outputFiles = [$crab3outfiles]
config.JobType.pyCfgParams=['type=${TYPE}','doTree=${DOTREE}','doExtraCalibTree=${DOEXTRACALIBTREE}','doExtraStudyTree=${DOEXTRASTUDYTREE}','doEleIDTree=${DOELEIDTREE}','doTreeOnly=1','isCrab=1','skim=${SKIM}','tagFile=${TAGFILE}','isPrivate=$ISPRIVATE','bunchSpacing=${BUNCHSPACING}', 'MC=${isMC}']

config.Data.inputDataset = '${DATASETPATH}'
config.Data.inputDBS = 'global'
EOF

if [[ ${isMC} = "0" ]]; then
cat >> ${crab3File} <<EOF
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = '${JSONFILE}'
config.Data.runRange = '${RUNRANGE}'
EOF
elif [[ ${isMC} = "1" ]]; then
cat >> ${crab3File} <<EOF
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
EOF
fi

cat >> ${crab3File} <<EOF
config.Data.publication = False
#config.Data.outputDatasetTag = 'outputdatasetTag'
#config.Site.storageSite = '$STORAGE_PATH' #=> come si scrive su eos??
#config.Site.storageSite = 'srm-eoscms.cern.ch'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_IT_Rome'
config.Data.outLFNDirBase = '${USER_REMOTE_DIR}' 

EOF


##############At this point you wrote the crab config file in tmp/${crabFile}####################

if [[ ${CRABVERSION} = "2" ]]; then
    echo "You are using LSFsubmit tool"
	create ${crabFile} ${FILELIST} || exit $?
fi

echo "Working dir is " $UI_WORKING_DIR
#in Crab 3 you directly submit (do not create the job first)
fi  #This is the end of the option createOnly

if [ -n "$SUBMIT" -a -z "${CHECK}" ];then
    if [[ ${CRABVERSION} = "2" ]]; then
		submit ${UI_WORKING_DIR} ${DATASETNAME} ${QUEUE}
    fi

    if [[ ${CRABVERSION} = "3" ]]; then
	echo "you are using crab 3"
	echo ${crabFile}

	##?? makeArguments.sh -f $FILELIST -u $UI_WORKING_DIR -n $FILE_PER_JOB || exit 1

	##?? splittedOutputFilesCrabPatch.sh -u $UI_WORKING_DIR
	crab submit -c  ${crabFile}
    fi

    STRING="${RUNRANGE}\t${DATASETPATH}\t${DATASETNAME}\t${STORAGE_ELEMENT}\t${USER_REMOTE_DIR_BASE}\t${TYPE}\t${TAG}\t${JSONNAME}"
    echo -e $STRING >> ntuple_datasets.dat

#else
    #echo "crab -c ${UI_WORKING_DIR} -submit"
fi

OUTFILES=`echo ${OUTFILES} | sed 's|,| |g'`

if [ -n "${CHECK}" ];then
    if [ ! -e "${UI_WORKING_DIR}/res/finished" ];then
	#echo $dir >> tmp/$TAG.log 
		check $UI_WORKING_DIR
#		echo "[STATUS] Unfinished ${UI_WORKING_DIR}"
#		resubmitCrab.sh -u ${UI_WORKING_DIR}
    else
		if [ "${isMC}" == "1" -a "${TYPE}" != "MINIAODNTUPLE" ];then
			OUTFILES="$OUTFILES PUDumper"
		fi
		for file in $OUTFILES
		do
			file=`basename $file .root`
#			echo "FILE $file"
			#( mergeOutput.sh -u ${UI_WORKING_DIR} -g $file --noRemove ) &
			case $file in
				*ntuple*)
					mergeOutput.sh -u ${UI_WORKING_DIR} -g $file  || exit 1
					;;
				*eleIDTree*)
					 mergeOutput.sh -u ${UI_WORKING_DIR} -g $file  || exit 1
					;;
				*extraCalibTree*)
					 mergeOutput.sh -u ${UI_WORKING_DIR} -g $file  || exit 1
					;;
				*)
#					time mergeOutput.sh -u ${UI_WORKING_DIR} -g $file --noRemove || exit 1
					;;
			esac
		done
		wait
    fi
#    echo "mergeOutput.sh -u ${UI_WORKING_DIR} -n ${DATASETNAME} -r ${RUNRANGE}"
fi

