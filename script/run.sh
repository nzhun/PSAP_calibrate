i=1;
for file in /home/local/ARCS/nz2274/PSAP_Setup/data/temp/*.txt
do
fp=$(echo $file|awk 'BEGIN{FS="/";}{print $NF}')
if [ $(( i % 70 )) == 0 ];
then
echo "wait v nohup.$fp.txt"
 nohup Rscript /home/local/ARCS/nz2274/PSAP_Setup/scirpt/R_opt_popscore.R $file > nohup.$fp.nohup.log
else
echo "no wait   nohup.$fp.txt"
 nohup Rscript /home/local/ARCS/nz2274/PSAP_Setup/scirpt/R_opt_popscore.R $file  > nohup.$fp.nohup.log &

fi
i=$((i+1))
done
