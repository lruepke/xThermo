ind_rhov=8
ind_rhol=17
ind_hl=20
ind_hv=11
ind_col=$ind_rhov
tol=100
Tmax=1000
mmc4_csv=../../../gallery_H2ONaCl/pT/mmc4_IAPS84.csv
awk -F ',' -v col=$ind_col 'NR==1{{col1=col-2;col2=col-1;};print $1, $2, $col1, $col2, $col}' $mmc4_csv

awk -F ',' -v col=$ind_col -v tol=$tol -v Tmax=$Tmax '{{col1=col-2;col2=col-1;};if(sqrt($col*$col)>tol && $1<Tmax)print $1, $2, $col1, $col2, $col, NR}' $mmc4_csv