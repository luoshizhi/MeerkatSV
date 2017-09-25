#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
use lib "$Bin/lib";
use AnaMethod;

my ($i, $sample_pair, $method, $outdir, $move, $monitorOption, $help, $qsubMemory, $pymonitor);

GetOptions(
	"i:s" => \$i,
	"outdir:s" => \$outdir,
	"move:s" => \$move,
	"method" => \$method,
	"m:s" => \$monitorOption,
	"help|?" => \$help,
	"qsubMemory:s" => \$qsubMemory,
	"pymonitor:s" => \$pymonitor,

	
);

my $usage = <<USE;
Usage:
description:Calling SV
author: Luoshizhi, luoshizhi\@genomics.cn
version :beta.1
date: 2017-08-15
usage: perl $0 [options]
	Common options:
	-i*		<str>	all tumor normal bam file list. format:sample_name  normal  tumor
	-outdir		<str>	outdir.[./]
	-method		<tr>	Calling SVs method: Meekat.[Meekat]
	-move		<str>	if this parameter is set,final result will be moved to it from output dir.
	-m		<str>	monitor options. will create monitor shell while defined this option
	-qsubMemory	<str>	"25G,10G,5G,1G"
	-help|?			print help information

	Software options:
	-pymonitor	<str>	monitor path [\$Bin/bin/monitor]
	
e.g.:
	perl $0 -i bam.list  -outdir ./outdir
USE
die $usage unless ($i && $outdir);
mkdir $outdir if  !-d $outdir;
$outdir=ABSOLUTE_DIR($outdir);

$method ||="Meerkat";
$monitorOption ||="taskmonitor -P common -p meekat_test  -q bc.q";
$qsubMemory ||= "25G,5G,5G,1G";
my @qsubMemory = split /,/,$qsubMemory;
$qsubMemory[0] ||= "25G";
$qsubMemory[1] ||= "10G";
$qsubMemory[2] ||= "5G";
$qsubMemory[3] ||= "1G";
$pymonitor ||= "$Bin/bin/monitor";
my ($shell, $process, $list)=("$outdir/shell/$method", "$outdir/process/$method", "$outdir/list");
mkpath $shell;mkpath $process;mkpath $list;

my $PATH="/ifshk4/BC_PUB/biosoft/pipe/bc_tumor/newblc/Pipeline/bin/:$Bin/bin/:\$PATH";
#my $PATH="/share/backup/luowen/install/meerkat/Meerkat/bin/:\$PATH";
my $LD_LIBRARY_PATH="$Bin/src/mybamtools/lib/:\$LD_LIBRARY_PATH";
#my $LD_LIBRARY_PATH="/share/backup/xiongteng/Software/Meerkat/src/mybamtools/lib/:\$LD_LIBRARY_PATH";
my $BLASTDIR="$Bin/scr/blast-2.2.26/bin/";
#my $BLASTDIR="/share/app/blast-2.2.26/bin";
my $BLASTDATADIR="$Bin/DB/blast-2.2.26/data/";
#my $BLASTDATADIR="/share/app/blast-2.2.26/data";

my $meerkat="$Bin/scripts/";
#my $meerkat="/share/backup/xiongteng/Software/Meerkat/scripts";
my $bwa="$Bin/src/bwa-0.6.2/";
#my $bwa="/share/backup/xiongteng/Software/bwa-0.6.2";
my $samtools="$Bin/src/samtools-0.1.19/";
#my $samtools="/share/backup/xiongteng/Software/samtools-0.1.19";
my $blast="$Bin/scr/blast-2.2.26/bin";
#my $blast="/share/app/blast-2.2.26/bin/";
my $meerkat_DB="$Bin/DB";
#my $meerkat_DB="/share/backup/xiongteng/Software/Meerkat/meerkat_DB";
my $hg19_bioDB="$meerkat_DB/hg19_bwa_idx/";
my $PERL5LIB="$Bin/src/perllib:\$PERL5LIB";
#my $PERL5LIB="/share/app/perl-5.22.0/lib/site_perl/5.22.0/:\$PERL5LIB";
my $env =<<"	ENVs.";
export PATH=$PATH && \\
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH && \\
export BLASTDIR=$BLASTDIR && \\
export BLASTDATADIR=$BLASTDATADIR && \\
export PERL5LIB=$PERL5LIB &&\\
meerkat=$meerkat && \\
bwa=$bwa && \\
samtools=$samtools && \\
blast=$blast && \\
meerkat_DB=$meerkat_DB && \\
hg19_bioDB=$hg19_bioDB && \\
	ENVs.



my $dependence = "$list/${method}_dependence.txt";
open TXT, ">$dependence" or die $!;
open I ,$i or die $!;
while(<I>){
	chomp;
	my ($sn,$normal,$tumor) = split /\s+/,$_;
	mkdir "$process/$sn";
	my $normal_basename = basename($normal);
	my $tumor_basename = basename ($tumor);


###step1 pre_process


my $pre_process = "$shell/pre_process_$sn.sh";
my $meerkat_c="$shell/meerkat_$sn.sh";
my $content="$env";
$content .="ln -fs $normal $process/$sn/ && \\\n";
$content .="ln -fs $normal.bai $process/$sn/ && \\\n";
$content .="ln -fs $tumor $process/$sn/ && \\\n";
$content .="ln -fs $tumor.bai $process/$sn/ &&\\\n";
$normal="$process/$sn/$normal_basename";
$tumor="$process/$sn/$tumor_basename";
$normal =~ s/.bam//g;
$tumor =~ s/.bam//g;
$content .="perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $normal.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools && \\\n";
$content .="perl $meerkat/pre_process.pl -t 4 -s 20 -k 1500 -q 15 -b $tumor.bam -I $hg19_bioDB/hg19.fasta -A $hg19_bioDB/hg19.fasta.fai -W $bwa -S $samtools";
AnaMethod::generateShell($pre_process,$content);
print TXT "$pre_process:$qsubMemory[2]\t$meerkat_c:$qsubMemory[2]\n";

###step2 meerkat
my $mechanism="$shell/mechanism_$sn.sh";

$content="$env ";
$content .= "mv $tumor.blacklist.gz $tumor.blacklist.real.gz &&\\\n";
$content .= "ln $normal.blacklist.gz $tumor.blacklist.gz &&\\\n";

$content .= "perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $normal.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast &&\\\n";
$content .= "perl $meerkat/meerkat.pl -s 20 -p 3 -o 1 -Q 10 -d 5 -t 8 -b $tumor.bam -F $hg19_bioDB -W $bwa -S $samtools -B $blast ";
AnaMethod::generateShell($meerkat_c,$content);
print TXT "$meerkat_c:$qsubMemory[2]\t$mechanism:$qsubMemory[2]\n";

###step 3 mechanism
my $somatic_calling="$shell/somatic_calling_$sn.sh";
$content  ="$env ";
$content .= "perl $meerkat/mechanism.pl -b $normal.bam -R $meerkat_DB/hg19_rmsk.txt && \\\n";
$content .= "perl $meerkat/mechanism.pl -b $tumor.bam -R $meerkat_DB/hg19_rmsk.txt ";
AnaMethod::generateShell($mechanism,$content);
print TXT "$mechanism:$qsubMemory[2]\t$somatic_calling:$qsubMemory[2]\n";

###step4 somatic_calling
my $somatica="$tumor.somatica.variants";
my $somaticb="$tumor.somaticb.variants";
my $somaticc="$tumor.somaticc.variants";
my $somaticd="$tumor.somaticd.variants";
my $somatice="$tumor.somatice.variants";
my $somaticf="$tumor.somaticf.variants";
my $somaticg="$tumor.somaticg.variants";

my $germline_calling="$shell/germline_calling_$sn.sh";
$content ="$env ";
$content .= "perl $meerkat/somatic_sv.pl -i $tumor.variants -o $somatica -R $meerkat_DB/hg19_rmsk.txt -F re_try/ -l 1000 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somatica -o $somaticb -R $meerkat_DB/hg19_rmsk.txt -n 1 -B $normal.bam -I $normal.isinfo  -D 5 -Q 10  -y 6 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticb -o $somaticc -R $meerkat_DB/hg19_rmsk.txt -u 1 -B $normal.bam -Q 10 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticc -o $somaticd -R $meerkat_DB/hg19_rmsk.txt -f 1 -B $normal.bam -Q 10 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticd -o $somatice -R $meerkat_DB/hg19_rmsk.txt -e 1 -B $tumor.bam -I $tumor.isinfo -D 5 -Q 10 -S $samtools&& \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somatice -o $somaticf -R $meerkat_DB/hg19_rmsk.txt -z 1 -S $samtools && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $somaticf -o $somaticg -R $meerkat_DB/hg19_rmsk.txt -d 40 -t 20 -S $samtools";
AnaMethod::generateShell($somatic_calling,$content);
print TXT "$somatic_calling:$qsubMemory[2]\t$germline_calling:$qsubMemory[2]\n";


###step 5 germline_calling

my $germa="$normal.germa.variants";
my $germb="$normal.germb.variants";
my $germc="$normal.germc.variants";
my $germd="$normal.germd.variants";
my $germe="$normal.germe.variants";

my $annotation="$shell/annotation_$sn.sh";
$content ="$env ";
$content .= "perl $meerkat/somatic_sv.pl -i $normal.variants -o $germa -R $meerkat_DB/hg19_rmsk.txt -E 0 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germa -o $germb -R $meerkat_DB/hg19_rmsk.txt -E 0 -e 1 -B $normal.bam -I $normal.isinfo  && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germb -o $germc -R $meerkat_DB/hg19_rmsk.txt -E 0 -u 1 -B $normal.bam && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germc -o $germd -R $meerkat_DB/hg19_rmsk.txt -E 0 -d 40 -t 40 && \\\n";
$content .= "perl $meerkat/somatic_sv.pl -i $germd -o $germe -R $meerkat_DB/hg19_rmsk.txt -E 0 -z 1 -p 5 -P 10 ";
AnaMethod::generateShell($germline_calling,$content);
print TXT "$germline_calling:$qsubMemory[2]\t$annotation:$qsubMemory[3]\n";

###step 6 annotation
$content ="$env ";
$content .= "perl $meerkat/fusions.pl -i $somaticg -G  $meerkat_DB/hg19_refGene.sorted.txt";
AnaMethod::generateShell($annotation,$content);

	
}
close TXT;
close I;
if(defined $pymonitor && defined $monitorOption){
	`echo "$pymonitor $monitorOption -i $dependence " >$list/${method}_qsub.sh`;
}


sub ABSOLUTE_DIR
#$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd `;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd `;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd `;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
