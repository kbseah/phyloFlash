#!/usr/bin/env perl

# Add records to pfrun database

=head1 NAME

B<phyloFlash_sqliteDB_add.pl> - Add completed phyloFlash runs to pFrrn database

=head1 SYNOPSIS

phyloFlash_sqliteDB_add.pl -db [SQLite db] -tgz [phyloFlash results archive] -run [run accession]

phyloFlash_sqliteDB_add.pl -man

=head1 DESCRIPTION

Add phyloFlash results to the SQLite database. This requires the tar.gz archive
from phyloFlash (produced with the -zip option). The library name will be
parsed from the phyloFlash results file. You should also supply an additional
run name, e.g. when comparing different phyloFlash results with the same file
name.

=cut

use strict;
use warnings;
use 5.010;

use FindBin;
use lib "$FindBin::RealBin/.."; # Relative path to PhyloFlash.pm
use PhyloFlash qw(err msg @msg_log);

use DBI;
use Getopt::Long;
use File::Basename qw(fileparse);
use Archive::Tar;
use Pod::Usage;

my ($dbfile,$pfarchive,$run);
my %pfrun_data;
my %pfseq_data;
my $semode;

pod2usage(-verbose=>0) if !@ARGV;

GetOptions ("db=s" => \$dbfile,       # SQLite formatted pfrun database
            "tgz=s" => \$pfarchive,   # Path to phyloFlash output TGZ archive - library name will be parsed from filename
            "run=s" => \$run,         # Unique number for the read run
            "help" => sub { pod2usage(-verbose=>1); },
            "man" => sub { pod2usage(-verbose=>2); },
            ) or pod2usage("Invalid input options");

=head1 OPTIONS

=over 8

=item -db I<FILE>

SQLite database containing phyloFlash results. This is initialized with
I<phyloFlash_sqliteDB_setupdb.pl>

=item -tgz I<FILE>

phyloFlash tar.gz archive file with results to be added to the database. The
library name will be parsed from the archive filename.

=item -run I<STRING>

Name for the phyloFlash run. This should be unique for each result within the
database. 

=back

=cut

err ("No database given") if (!defined $dbfile);
err ("No archive given") if (!defined $pfarchive);
err ("No run ID given") if (!defined $run);

# Get phyloFlash library prefix from pfarchive filename
my ($pflib,$dirs,$suffix) = fileparse($pfarchive, ".phyloFlash.tar.gz");

msg "Reading archive for phyloFlash run $pflib";
# Read phyloFlash tgz results into Archive::Tar object
my $pftgz = Archive::Tar->new;
$pftgz->read($pfarchive);

# Hash run ID
$pfrun_data{'run'} = $run;

# Read data from phyloFlash.report.csv
my $pfreport = $pftgz->get_content("$pflib.phyloFlash.report.csv");
foreach my $pfreportline (split /\n/, $pfreport) {
    my @split = split /,/, $pfreportline;
    if ($split[0] =~ m/version/) { # A dispatch table would be more elegant 
        $pfrun_data{"pfversion"} = $split[1];
    } elsif ($split[0] =~ m/library name/) {
        $pfrun_data{"libname"} = $split[1];
    } elsif ($split[0] =~ m/input reads/) {
        $pfrun_data{"input_reads"} = $split[1];
    } elsif ($split[0] =~ m/mapped SSU reads/) {
        $pfrun_data{"mapped_SSU_reads"} = $split[1];
    } elsif ($split[0] =~ m/single ended mode/) {
        if ($split[1] == 0) {
            $pfrun_data{"read pairing"} = 'paired';
        } else {
            $pfrun_data{"read pairing"} = 'single';
        }
    }
}

# Check that assembly files are in the archive
my $filechecks;
$filechecks->{'assemratio'} = $pftgz->contains_file("$pflib.assemratio.csv");
$filechecks->{'extractedSSUclass'} = $pftgz->contains_file("$pflib.phyloFlash.extractedSSUclassifications.csv");
$filechecks->{'finalfasta'} = $pftgz->contains_file("$pflib.all.final.fasta");
$filechecks->{'mapratio'} = $pftgz->contains_file("$pflib.mapratio.csv");

if ($filechecks->{'assemratio'}) {
    msg "... Archive contains $pflib.assemratio.csv";
    # Read data from assemratio.csv
    my $pfassemratio = $pftgz->get_content("$pflib.assemratio.csv");
    foreach my $pfarline (split /\n/, $pfassemratio) {
        my @split = split /,/, $pfarline;
        if ($split[0] =~ m/Unassembled/) {
            $pfrun_data{"unassembled_SSU_read_segments"} = $split[1];
        } elsif ($split[0] =~ m/Assembled/) {
            $pfrun_data{"assembled_SSU_read_segments"} = $split[1];
        }
    }
} else {
    msg "... WARNING: File $pflib.assemratio.csv not in archive: Perhaps no full-length seqs were assembled?";
}

if ($filechecks->{'mapratio'}) {
    msg "... Archive contains $pflib.mapratio.csv";
    my $pfmapratio = $pftgz->get_content("$pflib.mapratio.csv");
    foreach my $pfmapline (split /\n/, $pfmapratio) {
        my @split = split /,/, $pfmapline;
        $pfrun_data{'mapped_SSU_read_segments'} += $split[1] if $split[0] =~ m/Mapped/;
    }
} else {
    msg "... WARNING: File $pflib.mapratio.csv not in archive. Check phyloFlash version >3.0";
}


if ($filechecks->{'extractedSSUclass'}) {
    # Read data from phyloFlash.extractedSSUclassifications.csv
    my $pfclass = $pftgz->get_content("$pflib.phyloFlash.extractedSSUclassifications.csv");
    foreach my $pfclassline (split /\n/, $pfclass) {
        next if $pfclassline =~ m/^OTU,read_cov/; # skip header
        my @split = split /,/, $pfclassline;
        my $id = $split[0];
        %{$pfseq_data{$id}} = ('libname' => $pflib,
                               'read_cov' => $split[1],
                               'coverage' => $split[2],
                               'dbhit'    => $split[3],
                               'taxonomy' => $split[4],
                               'pid'      => $split[5],
                               'alnlen'   => $split[6],
                               );
        # From sequence ID get tool 
        if ($id =~ m/PFemirge/) {
            $pfseq_data{$id}{'tool'} = 'emirge';
        } elsif ($id =~ m/PFspades/) {
            $pfseq_data{$id}{'tool'} = 'spades';
        }
    }
} else {
    msg "... WARNING: File $pflib.phyloFlash.extractedSSUclassifications.csv not in archive: Perhaps no full-length seqs were assembled?";
}

if ($filechecks->{'finalfasta'}) {
    msg "... Reading in sequence data";
    # Read sequence data from all.final.fasta
    my $pffasta = $pftgz->get_content("$pflib.all.final.fasta");
    my $seq_concat;
    my $prev_seq;
    foreach my $pffastaline (split /\n/, $pffasta) {
        if ($pffastaline =~ m/^>(\S+)_[\d\.]+$/) {
            # If encounter header line
            my $seqid = $1;
            $pfseq_data{$prev_seq}{'seq'} = $seq_concat if defined $prev_seq;
            # Update sequence and clear seq_concat
            $prev_seq = $seqid;
            $seq_concat = undef;
        } else {
            # Concatenate sequence
            $seq_concat .= $pffastaline;
        }
    }
    # Sweep up last entry
    $pfseq_data{$prev_seq}{'seq'} = $seq_concat;
    
    # Hash run id
    foreach my $seqid (keys %pfseq_data) {
        $pfseq_data{$seqid}{'run'} = $run;
    }
} else {
    msg "... WARNING: File $pflib.all.final.fasta not in archive: Perhaps no full-length seqs were assembled?";
}

# Get sequence lengths and numbers of Ns
foreach my $seqid (keys %pfseq_data) {
    $pfseq_data{$seqid}{'len'} = length($pfseq_data{$seqid}{'seq'});
    $pfseq_data{$seqid}{'ns'} = () = $pfseq_data{$seqid}{'seq'} =~ m/N/g;
}


# Connect to database
msg "Connecting to database $dbfile";
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $user = "";
my $password = "";
my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
    PrintError => 0,
    RaiseError => 1,
    AutoCommit => 0,
    FetchHashKeyName => 'NAME_lc'
}) or die ("Cannot connect to database $dbfile: $!");

# Add new record to "pfrun" table in DB
msg "... Creating new pfrun record for $pflib";
my @pfrun_out;
my @params = qw(libname
                run
                pfversion
                input_reads
                mapped_SSU_reads
                mapped_SSU_read_segments
                assembled_SSU_read_segments
                unassembled_SSU_read_segments);
foreach my $param (@params) {
    push @pfrun_out, $pfrun_data{$param};
}
my $insertpfrun = $dbh->prepare("INSERT INTO pfrun (libname, run, pfversion, input_reads, mapped_SSU_reads, mapped_SSU_read_segments, assembled_SSU_read_segments, unassembled_SSU_read_segments, added) VALUES (?, ?, ?, ?, ?, ?, ?, ?, date('now'))");
$insertpfrun->execute(@pfrun_out);
$dbh->commit() or die $dbh->errstr;

if ($filechecks->{'finalfasta'}) {
    my $insertpfseq = $dbh->prepare('INSERT INTO pfseq (pfseq_id, libname, run, tool, seq, len, ns, read_cov, coverage, dbhit, taxonomy, pid, alnlen) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
    # For each sequence, add new record to "pfseq" table in DB
    my $counter = 0;
    foreach my $pfseq_id (keys %pfseq_data) {
        msg "... Creating new pfseq record $pfseq_id";
        $counter++;
        my @pfseq_out;
        push @pfseq_out, $pfseq_id;
        my @params = qw(libname run tool seq len ns read_cov coverage dbhit taxonomy pid alnlen);
        foreach my $param (@params) {
            push @pfseq_out, $pfseq_data{$pfseq_id}{$param};
        }
        $insertpfseq->execute(@pfseq_out);
        $dbh->commit() or $dbh->errstr if $counter % 500 == 0; # Commit every 500 records
    }
} else {
    msg "... WARNING: No sequences added to database as none found in phyloFlash archive file";
}

# Disconnect DB
$dbh->commit() or $dbh->errstr;
$dbh->disconnect;


## Reminders ##

#pfrun
#* pfrun_id (key) - from commandline (libname)
#* run (xref to run table) - from commandline
#* pfversion (pF version) - phyloFlash.report.csv
#* input_reads  - phyloFlash.report.csv
#* mapped_SSU_reads - phyloFlash.report.csv
#* assembled_SSU_reads - assemratio.csv
#* unassembled_SSU_reads - assemratio.csv
#
#pfseq
#* pfseq_id (key)
#* pfrun_id (xref to pfrun)
#* run (xref to run)
#* tool (spades or emirge) - parse from record
#* seq (16S rRNA sequence) - all.final.fasta
#* read_cov - phyloFlash.extractedSSUclassifications.csv
#* coverage
#* dbhit
#* taxonomy
#* pid
#* alnlen

