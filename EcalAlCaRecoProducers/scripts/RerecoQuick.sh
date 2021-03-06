#!/bin/bash
# Semplified script for rereco and ntuple of 2011 and 2012 data
##############
NTUPLE_REMOTE_DIR=group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples
ALCARERECO_REMOTE_DIR=group/dpg_ecal/alca_ecalcalib/ecalelf/alcarereco
#or defined as `echo $dataset |  sed 's|.*--remote_dir ||;s| .*$||;s|alcaraw|alcarereco|'`
SCHEDULER=caf
CREATE=y  #not used
SUBMIT=y  #not used
NTUPLE=y  #not used
RERECO=y  #not used
fileList=alcaraw_datasets.dat
crabVersion=2
#DOEXTRACALIBTREE="--doExtraCalibTree --doEleIDTree"
DOEXTRACALIBTREE=""
JSONNAME=""
JSONFILE=""

usage(){
    echo "`basename $0` -p period -t tagFile"
    echo "    -t | --tag_file tagFile"
    echo "    --scheduler caf,lsf (=$SCHEDULER)"
    echo "--------------- FOR RERECO only"
    echo "    -p | --period <test | RUN2011AB | RUN2012AB | RUN2012ABC | Cal_Nov2012 | Cal_Dic2012>"
    echo "    --singleEle: if you want to rereco also single electron datasets"
    echo "--------------- FOR NTUPLE only"
    echo "    --json_name jsonName"
    echo "    --json jsonFile.root"
    echo "    --doExtraCalibTree"
    echo "    --doEleIDTree"
    echo "    --noStandardTree"
    echo "--------------- Expert"
    echo "    --alcarerecoOnly"
    echo "    --ntupleOnly"
    echo "    --submitOnly"
    echo "    --createOnly"
    echo "    --check"
    echo "    --crabVersion (=2)"
    echo "    --file_per_job arg: number of files to process in 1 job (=1)"
    echo "--------------- ECALELF Tutorial"
    echo "    --tutorial"
    echo "--------------- Expert Option"
    echo "    --ntupleCheck: In the case that only ntuple are being run because alrerereco was done previously"
}


#------------------------------ parsing
# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o ht:p: -l help,tag_file:,period:,scheduler:,singleEle,createOnly,submitOnly,check,ntupleCheck,alcarerecoOnly,ntupleOnly,crabVersion:,json_name:,json:,tutorial,doExtraCalibTree,doEleIDTree,noStandardTree,file_per_job:,weightsReco -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

set -- $options

while [ $# -gt 0 ]
do
    case $1 in
	-h|--help) usage; exit 0;;
	-t | --tag) TAGFILE=$2; echo "[OPTION `basename $0`] TAGFILE:$TAGFILE"; TAG=`basename ${TAGFILE} .py`; shift;;
	-p | --period) PERIOD=$2; shift;;
	--tutorial) 
	    echo "[OPTION] Entering in tutorial mode"
	    TUTORIAL="--tutorial";
	    ;;
	--scheduler) SCHEDULER=$2; 
	    case $SCHEDULER in 
		lsf)
		    ;;
		caf)
		    ;;
			remoteGlidein)
				;;
		*)
		    echo "[ERROR] Scheduler ${SCHEDULER} not defined: use only lsf or caf" >> /dev/stderr
		    exit 1
		    ;;
	    esac
	    shift;;
	--singleEle) echo "[OPTION `basename $0`] singleEle"; SINGLEELE=y;;
	--createOnly) echo "[OPTION `basename $0`] createOnly"; unset SUBMIT; EXTRAOPTION="--createOnly";;
	--submitOnly) echo "[OPTION `basename $0`] submitOnly"; unset CREATE; EXTRAOPTION="--submitOnly";;
	--check)      echo "[OPTION `basename $0`] checking jobs"; CHECK=y; EXTRAOPTION="--check";;
	--ntupleCheck)      echo "[OPTION `basename $0`] checking nutple jobs"; NUTPLECHECK=y; EXTRAOPTION="--ntupleCheck";;
	--alcarerecoOnly) echo "[OPTION `basename $0`] alcarerecoOnly"; unset NTUPLE; EXTRAEXTRAOPTION="--alcarerecoOnly";;
	--ntupleOnly) echo "[OPTION `basename $0`] ntupleOnly"; unset RERECO; fileList=alcarereco_datasets.dat;; #EXTRAEXTRAOPTION="--ntupleOnly";;
 	--crabVersion) crabVersion=$2;  shift;;
 	--json) JSONFILE="--json=$2";  shift;;
	--json_name) JSONNAME="--json_name=$2"; shift;;
	--doExtraCalibTree) DOEXTRACALIBTREE="${DOEXTRACALIBTREE} --doExtraCalibTree";;
	--doEleIDTree)      DOEXTRACALIBTREE="${DOEXTRACALIBTREE} --doEleIDTree";;
	--noStandardTree)   DOEXTRACALIBTREE="${DOEXTRACALIBTREE} --noStandardTree";;
 	--file_per_job)     DOEXTRACALIBTREE="${DOEXTRACALIBTREE} --file_per_job=$2"; shift ;;
	--weightsReco)      DOWEIGHTSRECO=y;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
    (*) break;;
    esac
    shift
done
#------------------------------ checking


if [ ! -r "$TAGFILE" ];then
    echo "[ERROR] TAGFILE not found or not readable" >> /dev/stderr
    exit 1
fi

# check on readability of json
# check json name

#UI_WORKING_DIR=prod_alcarereco


periods=`grep -v '#' $fileList | sed -r 's|[ \t]+|\t|g' | cut -f 7 | sed 's|,|\n|g' | sort | uniq`
for period in $periods
  do
  if [ "$PERIOD" == "$period" ];then
      break;
  fi
  unset period
done
if [ -z "$PERIOD" ];then
    echo "[ERROR] PERIOD $PERIOD not defined"
    exit 1
fi



if [ -n "${TUTORIAL}" ];then
    case ${PERIOD} in
	RUN2012A)
	    ;;
	*)
	    echo "[ERROR] With the tutorial mode, the only permitted period is:"
	    echo "        RUN2012A"
	    #echo "        Be sure to have it in alcaraw_datasets.dat and to have selected it using the parseDatasetFile.sh"
	    exit 1
	    ;;
    esac
    if [ -n "${SINGLEELE}" ];then
	echo "[WARNING] No single electron rereco in tutorial mode" >> /dev/stderr
	unset SINGLEELE
    fi
fi

if [ -n "${RERECO}" ];then
	datasets=`./scripts/parseDatasetFile.sh $fileList | grep VALID | sed 's|$|,|' | grep "${PERIOD},"`
else
	datasets=`./scripts/parseDatasetFile.sh $fileList | grep ${TAG}$`
fi
# set IFS to newline in order to divide using new line the datasets
IFS=$'\n'



for dataset in $datasets
  do
	

#    echo "=========================================================================================="
	if [ -z "${SINGLEELE}" -a "`echo $dataset | grep -c SingleElectron`" != "0" ];then continue; fi
	if [ "`echo $dataset | grep -c SingleElectron`" != "0" -a "`echo $dataset | grep -c ZElectron`" != "0" ];then continue; fi
#	echo " [INFO] Dataset $dataset"
	# if [ -n "${DOWEIGHTSRECO}" ];then 
	# 	if [ "`echo $dataset | grep -c weightsReco`" == "0" ]; then continue; fi
	# else
	# 	if [ "`echo $dataset | grep -c weightsReco`" != "0" ]; then continue; fi
	# fi
	echo "============================================================"
	if [ -n "${RERECO}" ];then
		if [ -n "${DOWEIGHTSRECO}" ];then
			./scripts/prodAlcarereco.sh  -t ${TAGFILE} \
				--scheduler=$SCHEDULER --weightsReco --doExtraCalibTree --doEleIDTree ${DOEXTRACALIBTREE} ${EXTRAOPTION} ${EXTRAEXTRAOPTION} \
  				${JSONFILE} ${JSONNAME} --crabVersion=${crabVersion} \
				${TUTORIAL} $dataset 
		else
			./scripts/prodAlcarereco.sh  -t ${TAGFILE} \
			--scheduler=$SCHEDULER --doEleIDTree --doExtraCalibTree ${DOEXTRACALIBTREE} ${EXTRAOPTION} ${EXTRAEXTRAOPTION} \
  				${JSONFILE} ${JSONNAME} --crabVersion=${crabVersion} \
				${TUTORIAL} $dataset 
		fi
	else
		./scripts/prodNtuples.sh  -t ${TAGFILE} --type=ALCARERECO \
  			${JSONFILE} ${JSONNAME} \
			--scheduler=$SCHEDULER --file_per_job=10   --doEleIDTree --doExtraStudyTree --doExtraCalibTree --extraName=pedNoise  ${DOEXTRACALIBTREE} ${EXTRAOPTION} ${EXTRAEXTRAOPTION} \
			${TUTORIAL} $dataset 
#--scheduler=$SCHEDULER --file_per_job=10  --doExtraCalibTree --doEleIDTree --doExtraStudyTree --extraName=withExtraStudy  ${DOEXTRACALIBTREE} ${EXTRAOPTION} ${EXTRAEXTRAOPTION} \
		
	fi
	
done


exit 0

