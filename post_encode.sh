#!/bin/bash
 
#this script makes the sample.csv file for R pacakge DiffBinds.
#also create a script to make a QC table for all samples.
 
#rm -rf make_QC.sh
#rm -rf sample.csv
 
#echo "#!/bin/bash" > make_QC.sh
#echo -n "qc2tsv --collapse-header " >> make_QC.sh
#conditions=( PAX3-13 PAX3-27 )
analysis_dir=$1
#dir in root
 
mkdir -p result
 
#croo to gather results of ${analysis_dir} pipelines
#on Franklin croo requires ENCODE pipeline to be loaded
cd result
mkdir -p file_transfer
#write header for sample.csv
echo "SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller" > sample.csv 

#for dir in `ls -d ../${analysis_dir}/*`; do
for dir in `find ${analysis_dir}/ -maxdepth 1 -type d | tail -n+2`; do
  dir=`basename $dir`
  mkdir -p $dir
  cd $dir
  find ../../${analysis_dir}/${dir} -maxdepth 3 -name "metadata.json" -exec croo {} \;
  nrep=`grep "rep" ../../${analysis_dir}/${dir}/example.json  | grep "R1" | wc -l`
  #now write files to file_transfer
  crootable=`ls croo.filetable*`

  mkdir -p ../file_transfer/${dir}

  for (( j=1; j <= $nrep; j++ )); do
    #get filtered bam file
    echo "processing replicate $j..."
    link=`grep "Alignment/Replicate $j/Filtered BAM" $crootable | cut -f2`
    target_bam=`readlink -f $link`
    bam_loc=`echo "../file_transfer/${dir}/${dir}_rep${j}.bam" | xargs -I {} readlink -f {}`
    cp $target_bam $bam_loc
    echo "moving BAM file ${dir}_rep${j}.bam..."
    #sample ID should be everything before the S1/S2
    sampleID=`basename $link | tr "_" "\n" | sed '/S[0-9]\+/Q' | head -c -1 | tr "\n" "_"`
    #get peaks
    link=`grep "Peak/Replicate $j/Blacklist-filtered narrowpeak" $crootable | grep ".narrowPeak.gz" | cut -f2`
    target_peaks=`readlink -f $link`
    bed_loc=`echo "../file_transfer/${dir}/${dir}_rep${j}.bed" | xargs -I {} readlink -f {}`
    zcat $target_peaks > $bed_loc
    echo "moving peak file ${dir}_rep${j}.bed"
    #sampleID=`basename $link | tr "_" "\n" | grep $dir`
    #sampleID=`basename $target_peaks | tr "_" "\n" | sed '/trim/q' | tr "\n" "_"`

    #FC and p-val bigwig track
    echo "moving bigwig files..."
    grep "bigwig" $crootable  | grep fold-enrichment | grep "Replicate $j" | cut -f2 | xargs -I {} cp {} ../file_transfer/${dir}/${dir}_rep${j}_FC.bigwig
    grep "bigwig" $crootable  | grep p-val | grep "Replicate $j" | cut -f2 | xargs -I {} cp {} ../file_transfer/${dir}/${dir}_rep${j}_pval.bigwig

    #find peak caller, spp or MACS2. need ignore case
    if [[ `ls croo.report* | xargs -I {} sh -c "grep -i macs {}" |  grep "Replicate-$j" | wc -l` -gt 1 ]]; then
      peakcaller=macs
    else
      peakcaller=spp
    fi

    echo "$sampleID,$dir,$j, $bam_loc, $bed_loc, $peakcaller" >> ../sample.csv
  done

    #if number of replicates > 1
  if [[ $nrep -gt 1 ]]; then
    grep "Peak/IDR reproducibility/Conservative peak" $crootable | grep "narrowPeak.gz" | cut -f2 | xargs -I {} zcat {} > ../file_transfer/${dir}/${dir}_IDR_Conservative.bed
    echo "mvoing IDR peak file ${dir}_IDR_Conservative.bed"
  fi

  #QC files
  cp qc/qc.html ../file_transfer/${dir}/${dir}_qc.html
  echo "moving QC file ${dir}_qc.html"

  echo "moving JSON file"
  #add a line for moving json
  cp ../../${analysis_dir}/${dir}/example.json ../file_transfer/${dir}/example.json
  cd ..

done

#working directory is now result

#diffbind analysis if nrep>1
if $nrep -gt 1; then
  diffbind_ID=(sbatch --job-name="diffbind.R" --wrap="R CMD BATCH --no-save diffBind.R")

#get reference genome for motif analysis
if grep -q "mm10" ../file_transfer/${dir}/example.json; 
  then 
    genome="mm10"
  else
    genome="hg38"
  fi

#motif analysis
ml load homer/4.11.1

mkdir -p motif
mkdir file_transfer/motif
cd motif
for dir in `ls -d ../file_transfer/ | grep -v motif`; do
  dir=`basename $dir`
  mkdir -p ${dir}
  bedfile=`find ../file_transfer/${dir} -name *IDR*`
  #running motif and annotation
  sbatch --job-name=$dir --output=homer_${dir}.log --export=dir=$dir --wrap="findMotifsGenome.pl ${bedfile} ${genome} ${dir} -size given
                                                                            annotatePeaks.pl ${bedfile} ${genome} >  ${dir}/${dir}_annotation.txt
                                                                            cp ${dir}/${dir}_annotation.txt ../file_transfer/motif/${dir}_annotation.txt
                                                                            cp ${dir}/knownResults.html ../file_transfer/motif/${dir}_knownResults.html"                                           
done


#motif for differntial accessibility/binding sites
sbatch --dependency=afterok:${diffbind_ID##* } --job-name="Motif_Differential_sites" --wrap

sbatch  --wrap="findMotifsGenome.pl SDHB_Naive_WT_Naive_sig.txt hg38 . -size given -bg SDHB_Naive_WT_Naive_nonsig.txt
                 annotatePeaks.pl SDHB_Naive_WT_Naive_sig.txt hg38 >  Naive_annotation.txt"
             


#unload homer
ml unload homer/4.11.1
