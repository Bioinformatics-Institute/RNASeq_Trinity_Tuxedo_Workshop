#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#################################################################################
#
#  --autopilot         automatically run the pipeline end-to-end
#
#################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
    );

if ($help_flag) {
    die $usage;
}

my $BASEDIR = $FindBin::Bin;
my $RNASEQ_DATA_DIR = "$BASEDIR/RNASEQ_data";

my %RNASEQ_DATASETS = ( 'Sp_ds' => [ "$RNASEQ_DATA_DIR/Sp_ds.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_ds.right.fq.gz" ],
                        'Sp_hs' =>  [ "$RNASEQ_DATA_DIR/Sp_hs.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_hs.right.fq.gz" ],
                        'Sp_log' => [ "$RNASEQ_DATA_DIR/Sp_log.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_log.right.fq.gz" ],
                        'Sp_plat' => [ "$RNASEQ_DATA_DIR/Sp_plat.left.fq.gz", "$RNASEQ_DATA_DIR/Sp_plat.right.fq.gz" ],
    );

my $trinity_dir = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
$ENV{PATH} .= ":$trinity_dir";  ## adding it to our PATH setting.



my $OS_type = `uname`;

## first check for tools needed.

my @tools = qw (Trinity
    bowtie
    bowtie2
    tophat2
    samtools
    igv
);

{
    my $missing_tool_flag = 0;
    foreach my $tool (@tools) {
        my $path = `which $tool`;
        unless ($path =~ /\w/) {
            print STDERR "Error, cannot find path to: $tool\n";
            $missing_tool_flag = 1;
        }
    }
    
    if ($missing_tool_flag) {
        die "\n\nTools must be in PATH setting before proceeding.\n\n";
    }

}


if (0) { 
    ## unzip the gzipped files.
    foreach my $file (<*.gz>) {
        my $unzipped_file = $file;
        $unzipped_file =~ s/\.gz//;
        unless (-s $unzipped_file) {
            my $ret = system("gunzip -c $file > $unzipped_file");
            if ($ret) {
                die "Error, could not gunzip file $file";
            }
        }
    }
}


# Run Trinity.
my @left_fqs;
my @right_fqs;
foreach my $sample_type (keys %RNASEQ_DATASETS) {
    my ($left_fq, $right_fq) = @{$RNASEQ_DATASETS{$sample_type}};
    push (@left_fqs, $left_fq);
    push (@right_fqs, $right_fq);
}

my $checkpoints_dir = $FindBin::Bin . "/__TrinDemo_checkpoints_dir";
unless (-d $checkpoints_dir) {
    mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
}


my $run_Trinity_cmd = "Trinity --seqType fq --SS_lib_type RF "
    . " --left " . join(",", @left_fqs)
    . " --right " . join(",", @right_fqs)
    . " --CPU 2 --max_memory 1G";
&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");

## TODO:  examine the trinity gmap alignments.


## align the rna-seq reads against the genome, too, for comparison 
&process_cmd("bowtie2-build GENOME_data/genome.fa genome", "$checkpoints_dir/bowtie2_build_genome.ok");

my $tophat_align_cmd = "tophat2 -I 300 -i 20 genome " 
    . join(",", @left_fqs) . " " 
    . join(",", @right_fqs);

&process_cmd($tophat_align_cmd, "$checkpoints_dir/tophat_align_genome.ok");
&process_cmd("samtools index tophat_out/accepted_hits.bam", "$checkpoints_dir/index_tophat_bam.ok");


# use IGV

my $cmd = "igv -g GENOME_data/genome.fa GENOME_data/genes.bed,tophat_out/accepted_hits.bam";

if ($AUTO_MODE) {
    $cmd .= " & ";
}
&process_cmd($cmd, "$checkpoints_dir/igv.view_all.ok");




exit(0);

####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    if (-e $checkpoint) { return; }

    
    unless ($AUTO_MODE) {
        
        my $response = "";
        while ($response !~ /^[YN]/i) {
            print STDERR "\n\n"
                . "###############################################\n"
                . "CMD: $cmd\n"
                . "###############################################\n\n"
                . "Execute (Y/N)? ";

            $response = <STDIN>;
        }

        if ($response =~ /^N/i) {
            print STDERR "\t *** Exiting on demand. ****\n\n"
                . "Goodbye. \n\n\n";
            exit(0);
        }
    }
    
    print STDERR "\tNow running:\n\t\t$cmd\n\n\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    system("touch $checkpoint");
    
    return;
}


sub show {
    my ($image) = @_;

    my $cmd;

    if ($OS_type =~ /linux/i) {
        ## use evince
        $cmd = "evince $image";
    }
    else {
        ## rely on ImageMagick:
        $cmd = "open $image";
    }
    
    if ($AUTO_MODE) {
        $cmd .= " & ";
    }
    
    &process_cmd($cmd, "$checkpoints_dir/view." . basename($image) . ".ok");

    return;
}
