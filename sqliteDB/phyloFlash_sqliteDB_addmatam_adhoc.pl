#!/usr/bin/env perl

# Add MATAM assembled sequences and Vsearch closest hits to SILVA to pF database

=head1 NAME

B<phyloFlash_sqliteDB_addmatam_adhoc.pl> - Add MATAM assembled sequences and Vsearch closest hits to SILVA to pF database

=head1 SYNOPSIS

phyloFlash_sqliteDB_addmatam_adhoc.pl -db [SQLite db] -libname [phyoFlash library name] -seq [matam fasta] -vsearch [vsearch table]

phyloFlash_sqliteDB_addmatam_adhoc.pl -man

=head1 DESCRIPTION

Add matam assembly results

=cut

use strict;
use warnings;
use 5.010;

use DBI;
use Getopt::Long;
use File::Basename qw(fileparse);
use Archive::Tar;
use Pod::Usage;

my ($dbfile,$libname,$matam_fasta,$matam_vsearch_tsv,$run,$overwrite);
my %pfseq_data;

pod2usage(2) if (!@ARGV);

GetOptions ('db=s' => \$dbfile,
            'seq=s' => \$matam_fasta,
            'vsearch=s' => \$matam_vsearch_tsv,
            'libname=s' => \$libname,
            'run=s' => \$run,
            'help' => sub { pod2usage(1) },
            'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
            ) or pod2usage(2);

=head1 INPUT ARGUMENTS

=over 8

=item -db FILE

Path to pFrrn SQLite database file

=item -libname STRING

Name of library

=item -seq FILE

Path to final_assembly.fa file produced by MATAM

=item -vsearch FILE

Path to TSV results from running vsearch usearch_global on MATAM assembled
sequences vs. SILVA database (for closest db hit of assembled seqs)

=item -help

Short help message

=item -man

Full manual page

=back

=cut

die ("No DB file specified") if (!defined $dbfile);
die ("No Fasta file specified") if (!defined $matam_fasta);
die ("No Vsearch results given") if (!defined $matam_vsearch_tsv);

# Connect to database
say "Connecting to database $dbfile";
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $user = "";
my $password = "";
my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
    PrintError => 0,
    RaiseError => 1,
    AutoCommit => 0,
    FetchHashKeyName => 'NAME_lc'
}) or die ("ERROR: Cannot connect to database $dbfile: $!");

# Read sequence data from MATAM fasta file
say "Reading in sequence data from MATAM assembly for library $libname";
open(my $fastafh, "<", $matam_fasta) or die ("Cannot open file $matam_fasta: $!");
my $seq_concat;
my $prev_seq;
while (my $fastaline = <$fastafh>) {
    chomp $fastaline;
    if ($fastaline =~ m/^>(\d+) count=([\d\.]+)$/) {
        # If encounter header line
        my $seqid = "$libname.matam_$1";
        $pfseq_data{$prev_seq}{'seq'} = $seq_concat if (defined $prev_seq);
        # Update sequence and clear seq_concat
        $prev_seq = $seqid;
        $seq_concat = undef;
    } else {
        # Concatenate sequence
        $seq_concat .= $fastaline;
    }
}
# Sweep up last entry
$pfseq_data{$prev_seq}{'seq'} = $seq_concat if defined $prev_seq; 

# Get sequence lengths and number of Ns
foreach my $key (keys %pfseq_data) {
    $pfseq_data{$key}{'len'} = length($pfseq_data{$key}{'seq'});
    $pfseq_data{$key}{'ns'} = () = $pfseq_data{$key}{'seq'} =~ m/N/g;
}


# Read data from vsearch output table
open(my $csvfh, "<", $matam_vsearch_tsv) or die ("Cannot open file $matam_vsearch_tsv: $!");
foreach my $classline (<$csvfh>) {
    next if $classline =~ m/^OTU,read_cov/; # skip header
    my @split = split /\t/, $classline;
    
    # Parse ID and Count from MATAM sequence ID
    my $id = $split[0];
    my ($seqid, $count);
    if ($id =~ m/(\d+) count=([\d\.]+)/) {
        $seqid = "$libname.matam_$1";
        $count = $2;
    } else {
        say STDERR "Entry $id does not appear to be a MATAM sequence ID, skipping...";
        next;
    }
    
    # Parse dbhit and taxonomy string from SILVA header
    my $silvahit = $split[1];
    my ($dbhit, $taxonstring);
    if ($silvahit =~ m/(\w+\.\d+\.\d+) (.+)/) {
        $dbhit = $1;
        $taxonstring = $2;
    }
    
    # Hash results
    $pfseq_data{$seqid}{'libname'}  = $libname;
    $pfseq_data{$seqid}{'run'}      = "$libname\_matam";
    $pfseq_data{$seqid}{'read_cov'} = $count;
    $pfseq_data{$seqid}{'coverage'} = $count;
    $pfseq_data{$seqid}{'dbhit'}    = $dbhit;
    $pfseq_data{$seqid}{'taxonomy'} = $taxonstring;
    $pfseq_data{$seqid}{'pid'}      = $split[2];
    $pfseq_data{$seqid}{'alnlen'}   = $split[3];
    $pfseq_data{$seqid}{'tool'}     = 'matam';
}


# For each sequence, add new record to "pfseq" table in DB
my $insert_pfseq = $dbh->prepare('INSERT INTO pfseq (pfseq_id, libname, run, tool, seq, read_cov, coverage, dbhit, taxonomy, pid, alnlen, ns, len) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)');
foreach my $pfseq_id (keys %pfseq_data) {
    say "Creating new pfseq record $pfseq_id";
    my @pfseq_out;
    push @pfseq_out, $pfseq_id;
    my @params = qw(libname run tool seq read_cov coverage dbhit taxonomy pid alnlen ns len);
    foreach my $param (@params) {
        push @pfseq_out, $pfseq_data{$pfseq_id}{$param};
    }
    $insert_pfseq->execute(@pfseq_out);
}
$dbh->commit() or $dbh->errstr;
say "Disconnecting from database";
$dbh->disconnect;

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 - Brandon Seah <kbseah@mpi-bremen.de>

LICENCE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
