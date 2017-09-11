#!/usr/bin/env perl

# Downloader script to test Fastq files from ENA and check if
# 1. they are not SSU rRNA amplicon libraries and
# 2. what is the read length

use 5.010; # For the 'say' command
use strict;
use warnings;

use LWP::Simple;    # to read data from URLs
use FindBin;
use lib $FindBin::RealBin;

# Get input from command line
my $acc = $ARGV[0]; # ENA accession number
my $CPUs = 8;

# Variables
my @fastq_urls; # List of Fastq URLs
my $nhmmer; # Nhmmer binary
my $wget;   # Wget binary
my $opsys = $^O; # Check operating system
say STDERR "Current operating system $opsys";
if ($opsys eq "darwin") {
    $nhmmer = "$FindBin::RealBin/barrnap-HGV/binaries/darwin/nhmmer";
} elsif ($opsys eq "linux") {
    $nhmmer = "$FindBin::RealBin/barrnap-HGV/binaries/linux/nhmmer";
} else {
    $nhmmer = "$FindBin::RealBin/barrnap-HGV/binaries/linux/nhmmer";
}
my $hmm = "$FindBin::RealBin/barrnap-HGV/db/ssu/ssu_ABE.hmm"; # HMM for SSU rRNA

## MAIN ########################################################################

@fastq_urls = get_fastq_urls($acc);
if (defined $fastq_urls[0]) {
    my $prop = test_SSU_amplicon($fastq_urls[0]);
    say join "\t", ($fastq_urls[0], $prop);
}

## SUBS ########################################################################

sub test_SSU_amplicon {
    # Test the first ten entries in Fastq files to check if it is an amplicon library

    # Input: URL of Fastq file from ENA
    my ($url) = @_;
    # Use wget, gzip, and head to get the first ten entries
    # Use backticks to capture STDOUT (STDERR ignored)
    my $preview = `wget -O - $url | gzip -cd | head -n40`;
    # Convert to Fasta file
    my @prevlines = split /\n/, $preview;
    open (my $fafh, ">", "tmp.$acc.fasta") or die ("Cannot open file: $!");
    for (my $i=0; $i<39; $i++) {
        if (($i+1)%4 == 1) {
            # Header line
            #my $header = $prevlines[$i] =~ s/\@/>/;
            my $header = $i;
            print $fafh ">$header";
            print $fafh "\n";
        } elsif (($i+1)%4 == 2) {
            # Sequence line
            print $fafh $prevlines[$i];
            print $fafh "\n";
        }
    }
    close ($fafh);

    # Test against HMM
    # Run nhmmer
    my @nhmmer_cmd = ($nhmmer,
                      "--cpu $CPUs",
                      "--tblout tmp.$acc.nhmmer.tblout",
                      $hmm,
                      "tmp.$acc.fasta");
    system (join " ", @nhmmer_cmd);
    # Parse output table
    my %heads;
    open (my $tblfh, "<", "tmp.$acc.nhmmer.tblout") or die ("Cannot open file: $!");
    while (my $line = <$tblfh>) {
        chomp $line;
        next if ($line =~ m/^#/); # Skip header/footer
        my @splitline = split /\s+/, $line;
        $heads{$splitline[0]}++;
    }
    close ($tblfh);
    # Report number of amplicon hits
    my $hits = scalar (keys %heads);
    return $hits;
}

sub get_fastq_urls {
    # Get the URL(s) of ENA Fastq file(s) for a given ENA entry

    # Input: Accession number of sample or run
    my ($acc) = @_;
    my @urls_arr;    # Output array containing URLs
    # Get report table using ENA REST query
    my $url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$acc&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes";
    my $tab = get $url;
    foreach my $line (split /\n/, $tab) {
        next if $line =~ m/^run_accession/; # Skip header
        my @splitline = split /\t/, $line;
        # Multiple URLs are separated by a semicolon
        my @spliturl = split /;/, $splitline[1];
        #my $run_accession = $splitline[0]; # run_accession, different from sample acc
        foreach my $fastq_url (@spliturl) {
            push @urls_arr, $fastq_url;
        }
    }
    # Return list of URLs
    return (@urls_arr);
}
