#!/bin/bash

ref_map.pl -T 8 -m 5 -S -b 1 -D "all samples, MAPQ >= 10, ref align to GigaDB G. fortis genome" -o ./output/ \
-s ./aligned_reads/ARK_CAH143.bam \
-s ./aligned_reads/ARK_CAH144.bam \
-s ./aligned_reads/ARK_CAH156.bam \
-s ./aligned_reads/ARK_CAH163.bam \
-s ./aligned_reads/FLA_FL28W.bam \
-s ./aligned_reads/FLA_FL30W.bam \
-s ./aligned_reads/FLA_FLWN26.bam \
-s ./aligned_reads/FLA_FLWN27.bam \
-s ./aligned_reads/FLA_FLWN29.bam \
-s ./aligned_reads/FLA_FLWN33.bam \
-s ./aligned_reads/FLA_PB52307.bam \
-s ./aligned_reads/FLA_PB52314.bam \
-s ./aligned_reads/FLA_PB52315.bam \
-s ./aligned_reads/FLA_PB52329.bam \
-s ./aligned_reads/FLA_PB52330.bam \
-s ./aligned_reads/FLA_PB52335.bam \
-s ./aligned_reads/GUA_DHB4503.bam \
-s ./aligned_reads/GUA_DHB4507.bam \
-s ./aligned_reads/GUA_DHB4551.bam \
-s ./aligned_reads/GUA_JK02010.bam \
-s ./aligned_reads/GUA_JK02015.bam \
-s ./aligned_reads/HON_JK00031.bam \
-s ./aligned_reads/JAL_KSW3233.bam \
-s ./aligned_reads/JAL_KSW3254.bam \
-s ./aligned_reads/JAL_KSW3273.bam \
-s ./aligned_reads/KNS_JK04545.bam \
-s ./aligned_reads/KNS_SAR7994.bam \
-s ./aligned_reads/KNS_SAR7995.bam \
-s ./aligned_reads/KNS_SAR7996.bam \
-s ./aligned_reads/LOU_CAH097.bam \
-s ./aligned_reads/LOU_CAH098.bam \
-s ./aligned_reads/LOU_CAH100.bam \
-s ./aligned_reads/LOU_CAH140.bam \
-s ./aligned_reads/LOU_JK04540.bam \
-s ./aligned_reads/NCL_PB65637.bam \
-s ./aligned_reads/NCL_PB65639.bam \
-s ./aligned_reads/NCL_PB65648.bam \
-s ./aligned_reads/NCL_PB65655.bam \
-s ./aligned_reads/NCL_PB65658.bam \
-s ./aligned_reads/OAX_ASJ90.bam \
-s ./aligned_reads/OAX_DHB5578.bam \
-s ./aligned_reads/OAX_DHB5582.bam \
-s ./aligned_reads/OAX_DHB5715.bam \
-s ./aligned_reads/OAX_JMD069.bam \
-s ./aligned_reads/OKL_CAH081.bam \
-s ./aligned_reads/OKL_CAH090.bam \
-s ./aligned_reads/OKL_JK04554.bam \
-s ./aligned_reads/OKL_JMD369.bam \
-s ./aligned_reads/QUI_GLS10.bam \
-s ./aligned_reads/QUI_GLS172.bam \
-s ./aligned_reads/QUI_GLS23.bam \
-s ./aligned_reads/QUI_MGL19.bam \
-s ./aligned_reads/QUI_PEP2969.bam \
-s ./aligned_reads/SCL_PB65537.bam \
-s ./aligned_reads/SCL_PB65545.bam \
-s ./aligned_reads/SCL_PB65556.bam \
-s ./aligned_reads/SCL_PB65557.bam \
-s ./aligned_reads/SCL_PB65564.bam \
-s ./aligned_reads/SIN_CSW7732.bam \
-s ./aligned_reads/SIN_CSW7736.bam \
-s ./aligned_reads/SIN_CSW7737.bam \
-s ./aligned_reads/SIN_CSW7763.bam \
-s ./aligned_reads/SIN_EAG016.bam \
-s ./aligned_reads/SIN_EAG024.bam \
-s ./aligned_reads/SIN_EAG032.bam \
-s ./aligned_reads/SIN_RCF2610.bam \
-s ./aligned_reads/TAB_CAM318.bam \
-s ./aligned_reads/TEX_BTS05071.bam \
-s ./aligned_reads/TEX_CAH087.bam \
-s ./aligned_reads/TEX_CAH146.bam \
-s ./aligned_reads/TEX_CAH149.bam \
-s ./aligned_reads/TEX_CAH153.bam \
-s ./aligned_reads/TEX_CAH155.bam \
-s ./aligned_reads/TEX_CAH161.bam \
-s ./aligned_reads/TEX_CAH162.bam \
-s ./aligned_reads/TEX_CAH170.bam \
-s ./aligned_reads/TEX_CAH177.bam \
-s ./aligned_reads/TEX_JK04518.bam \
-s ./aligned_reads/TEX_JK04519.bam \
-s ./aligned_reads/TEX_JK04563.bam \
-s ./aligned_reads/VCZ_TUX1094.bam \
-s ./aligned_reads/VCZ_TUX1107.bam \
-s ./aligned_reads/VCZ_TUX1402.bam \
-s ./aligned_reads/VCZ_TUX204.bam \
-s ./aligned_reads/YUC_BRB849.bam \
-s ./aligned_reads/YUC_BRB875.bam \
-s ./aligned_reads/YUC_BRB878.bam \
-s ./aligned_reads/YUC_BRB883.bam \
-s ./aligned_reads/YUC_BRB893.bam \
-s ./aligned_reads/YUC_BRB928.bam \
-s ./aligned_reads/YUC_BRB941.bam \
-s ./aligned_reads/YUC_BRB942.bam \
-s ./aligned_reads/YUC_CAM013.bam \
-s ./aligned_reads/YUC_CAM109.bam \
-s ./aligned_reads/YUC_CAM31.bam \