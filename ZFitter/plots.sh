#!/bin/bash

cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="M_{ee}"
set yrange  [88:92]
EOF

gnuplot -c macro/stability.gpl testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat ref testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat pedv1 testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_Ped_v2-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat pedv2 testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_PS_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat PSv1

cp stability.pdf stability-invMass.pdf


### REF vs PS_v1
cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="M_{ee}"
set yrange  [88:92]
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat ref  \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_PS_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat PSv1
cp stability.pdf stability-invMass-detailed.pdf


gnuplot -c macro/stability_validation.gpl  testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_PS_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/invMass_ECAL_ele-stability_runNumber.dat 
cp stability.pdf stability-invMass-validation.pdf


cat > v.gpl <<EOF
peakLabelUnit=""
VARIABLE="R_9"
set yrange  [0.9:0.95]
EOF

gnuplot -c macro/R9stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat ref \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_PS_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/R9Ele-stability_runNumber.dat PSv1
cp stability.pdf stability-R9Ele.pdf 

### REF vs PS_v1
cat > v.gpl <<EOF
peakLabelUnit="[GeV]"
VARIABLE="ES"
#set yrange  [88:92]
EOF

gnuplot -c macro/stability.gpl \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_ref-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/esEnergySCEle-stability_runNumber.dat ref  \
	testNew/dato/rereco/Cal_Oct2017/Cal_Oct2017_PS_v1-Z/loose25nsRun22016Moriond/invMass_ECAL_ele/anyVar/fitres/esEnergySCEle-stability_runNumber.dat PSv1

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

