#export PATH="/share/backup/luowen/install/meerkat/Meerkat/bin/:$PATH"
export PATH="/share/backup/xiongteng/Software/Meerkat/bin/:$PATH"
#export LD_LIBRARY_PATH=/share/backup/luowen/install/meerkat/Meerkat/src/mybamtools/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/share/backup/xiongteng/Software/Meerkat/src/mybamtools/lib/:$LD_LIBRARY_PATH
export BLASTDIR='/share/app/blast-2.2.26/bin'
export BLASTDATADIR='/share/app/blast-2.2.26/data'

meerkat=/share/backup/xiongteng/Software/Meerkat/scripts
bwa=/share/backup/xiongteng/Software/bwa-0.6.2
samtools=/share/backup/xiongteng/Software/samtools-0.1.19
#blast=/opt/bio/ncbi/bin/
blast=/share/app/blast-2.2.26/bin/
#blast=/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/ncbi-blast-2.5.0+/bin

meerkat_DB=/share/backup/xiongteng/Software/Meerkat/meerkat_DB
hg19_bioDB=/share/backup/xiongteng/Software/Meerkat/meerkat_DB/hg19_bwa_idx/

normal=
tumor=

##step1 pre_process
perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $normal.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools && echo end step1 at:`date`

perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $tumor.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools && echo end step1 at:`date`

##step2 meerkat
mv $tumor.blacklist.gz $tumor.blacklist.real.gz
ln $normal.blacklist.gz $tumor.blacklist.gz

perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $normal.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast && echo end dc at:`date`
perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $tumor.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast  && echo end dc at:`date`

##step3 mechanism
perl $meerkat/mechanism.pl -b $normal.bam -R $meerkat_DB/hg19_rmsk.txt && echo end step3 at:`date`
perl $meerkat/mechanism.pl -b $tumor.bam -R $meerkat_DB/hg19_rmsk.txt && echo end step3 at:`date`


##step4 somatic calling

somatica=$tumor.somatica.variants
somaticb=$tumor.somaticb.variants
somaticc=$tumor.somaticc.variants
somaticd=$tumor.somaticd.variants
somatice=$tumor.somatice.variants
somaticf=$tumor.somaticf.variants
somaticg=$tumor.somaticg.variants

normal_blacklistrg=""
tumor_blacklistrg=""
normal_bam=$normal.bam
tumor_bam=$tumor.bam

perl $meerkat/somatic_sv.pl -i $tumor.variants -o $somatica -R $meerkat_DB/hg19_rmsk.txt -F re_try/ -l 1000
perl $meerkat/somatic_sv.pl -i $somatica -o $somaticb -R $meerkat_DB/hg19_rmsk.txt -n 1 -B $normal_bam -I $normal.isinfo -D 5 -Q 10 -S $samtools
perl $meerkat/somatic_sv.pl -i $somaticb -o $somaticc -R $meerkat_DB/hg19_rmsk.txt -u 1 -B $normal_bam -Q 10 -S $samtools
perl $meerkat/somatic_sv.pl -i $somaticc -o $somaticd -R $meerkat_DB/hg19_rmsk.txt -f 1 -B $normal_bam -Q 10 -S $samtools
perl $meerkat/somatic_sv.pl -i $somaticd -o $somatice -R $meerkat_DB/hg19_rmsk.txt -e 1 -B $tumor_bam -I $tumor.isinfo -D 5 -Q 10 -S $samtools
perl $meerkat/somatic_sv.pl -i $somatice -o $somaticf -R $meerkat_DB/hg19_rmsk.txt -z 1 -S $samtools
perl $meerkat/somatic_sv.pl -i $somaticf -o $somaticg -R $meerkat_DB/hg19_rmsk.txt -d 40 -t 20 -S $samtools

echo end at: `date`

#step5 germline calling
germa=$normal.germa.variants
germb=$normal.germb.variants
germc=$normal.germc.variants
germd=$normal.germd.variants
germe=$normal.germe.variants

perl $meerkat/somatic_sv.pl -i $normal.variants -o $germa -R $meerkat_DB/hg19_rmsk.txt -E 0
perl $meerkat/somatic_sv.pl -i $germa -o $germb -R $meerkat_DB/hg19_rmsk.txt -E 0 -e 1 -B $normal_bam -I $normal.isinfo  
perl $meerkat/somatic_sv.pl -i $germb -o $germc -R $meerkat_DB/hg19_rmsk.txt -E 0 -u 1 -B $normal_bam
perl $meerkat/somatic_sv.pl -i $germc -o $germd -R $meerkat_DB/hg19_rmsk.txt -E 0 -d 40 -t 40 
perl $meerkat/somatic_sv.pl -i $germd -o $germe -R $meerkat_DB/hg19_rmsk.txt -E 0 -z 1 -p 5 -P 10 

echo end at: `date`

#step6 annotation 

perl $meerkat/fusions.pl -i $somaticg -G  $meerkat_DB/hg19_refGene.sorted.txt 

