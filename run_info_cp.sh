

smk_log=$(ls -lht .snakemake/log/* | head -n1 | awk '{print $9}')
runtime=$(basename $smk_log | cut -d. -f1)

#echo $smk_log
#echo $runtime

config=$1
outdir=$2
out_config=$runtime"."$(basename $config)
git_commit=$(git log --pretty=format:'%h' -n 1)
out_git=$runtime".git_commit.txt"

echo "cp $smk_log $outdir"
echo "cp $config $outdir/$out_config"
echo "echo $git_commit > $outdir/$out_git"

#user prompt to confirm whether to perform all moves
#https://stackoverflow.com/questions/226703/how-do-i-prompt-for-yes-no-cancel-input-in-a-linux-shell-script/27875395#27875395
echo 'Do you wish to proceed with all the above commands (put 1 or 2 to select yes or no)?'
select yn in "yes" "no"; do
  case $yn in
    yes ) echo 'you selected to perform these commands'; cp $smk_log $outdir && cp $config $outdir/$out_config;
     break;;
    no ) echo 'You have opted to not perform the above commands. Aborting...'; exit;;
  esac
done
