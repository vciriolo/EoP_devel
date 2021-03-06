#!/bin/bash

cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="M_{ee}"
set yrange  [84:92]
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_PS_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat PSv2 \
	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_Ped_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat Pedv1 \
#	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat refv1 \
#	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_IC_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat ICv1 \


cp stability.pdf stability-invMass.pdf
exit 0

cat > v.gpl <<EOF
peakLabelUnit=""
VARIABLE="R_9"
#set yrange  [0.9:0.98]
EOF

gnuplot -c macro/R9stability.gpl \
	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_PS_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat PSv2 \
	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_Ped_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat Pedv1 \
#	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat refv1 \
#	testNew/dato/rereco/Cal_Nov2017/Cal_Nov2017_IC_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat ICv1 \


cp stability.pdf stability-R9Ele.pdf 
exit 0
### REF vs cand_v1
cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="M_{ee}"
set yrange  [88:92]
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat candv1 \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat candv2 \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v3-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat candv3 \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v4-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat candv4 

#	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat ref  \

cp stability.pdf stability-invMass-detailed.pdf


gnuplot -c macro/stability_validation.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v4-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat 

cp stability.pdf stability-invMass-validation.pdf


pdfunite stability-invMass.pdf stability-invMass-detailed.pdf stability-invMass-validation.pdf stability-R9Ele.pdf stability.pdf
exit 0

### REF vs cand_v1
cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="ES"
#set yrange  [88:92]
unset yrange 
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/esEnergySCEle-stability_runNumber.dat PSv1 \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ES_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/esEnergySCEle-stability_runNumber.dat ESv2 \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_cand_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/esEnergySCEle-stability_runNumber.dat PSEoPv1 \

cp stability.pdf stability-ES-detailed.pdf


### Ped_v2 vs Ped_v3
cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="M_{ee}"
set yrange  [89:92]
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat Pedv2  \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v3-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat Pedv3
cp stability.pdf PCLvalidation.pdf


cat > v.gpl <<EOF
peakLabelUnit=" "
VARIABLE="R9"
#set yrange  [0.9:0.95]
EOF

gnuplot -c macro/R9stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat Pedv2  \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v3-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat Pedv3
cp stability.pdf PCLvalidation-R9Ele.pdf 

pdfunite stability-invMass.pdf stability-invMass-detailed.pdf stability-invMass-validation.pdf stability-R9Ele.pdf PCLvalidation.pdf PCLvalidation-R9Ele.pdf stability.pdf 

