#!/usr/bin/env perl

=head1 NAME

B<pFrrn_makeblastdb.pl> - Make Blast database from pFrrn database

=head1 SYNOPSIS

pFrrn_makeblastdb.pl -db [SQLite db] -blastdb [Blast db prefix] -title="[Title]"

pFrrn_makeblastdb.pl -man

=head1 DESCRIPTION

Create Blastn database of full-length SSU sequences produced by phyloFlash, from
a pFrrn SQLite database file.

=cut 

# Export pF SSU sequences from database and build blastdb

use strict;
use warnings;
use 5.010;

use DBI;
use Getopt::Long;
use File::Basename qw(fileparse);
use File::Temp qw(tempfile);
use Pod::Usage;

my $makeblastdb = "makeblastdb"; # Binary for makeblastdb
my $ENA_acc;
my ($dbfile,    # SQLite database
    $blastdb    # Blastdb path
    );
my $blastdbtitle;
my $minlen = 800;
my $maxNs = 10;

pod2usage(2) if (!@ARGV); # Help message 

GetOptions ("db=s" => \$dbfile,
            "blastdb=s" => \$blastdb,
            "ENA" => \$ENA_acc,
            "title=s" => \$blastdbtitle,
            "len=i" => \$minlen,
            "N=i" => \$maxNs,
            'help' => sub { pod2usage(1) },
            'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
            );

=head1 INPUT ARGUMENTS

=over 8

=item -db FILE

Path to pFrrn SQLite database file

=item -blastdb PREFIX

Filename prefix for Blast database files

=item -ENA

Parse accession numbers as ENA accession numbers (Default: No)

=item -title STRING

Title for Blast database

=item -len INTEGER

Minimum length of assembled/reconstructed SSU sequences to include in Blast db
(Default: 800)

=item -N INTEGER

Maximum number of Ns in assembled/reconstructed SSU sequences to allow, to
include sequences in Blast db (Default: 10)

=item -help

Short help message

=item -man

Full manual page

=back

=cut 


say STDERR "Creating blastdb $dbfile from pFrrn Sqlite db $dbfile";
say STDERR "from sequences with min length $minlen and max $maxNs Ns";

# If no blastdb title given, make one from $dbfile
my ($dbfilebase, $dbfiledirs, $dbfilesuffix) = fileparse($dbfile, (".db",".sqlite",".sql"));
if (!defined $blastdbtitle) {
    $blastdbtitle = $dbfilebase;
}

# Connect to database
say "Connecting to database $dbfile";
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $user = "";
my $password = "";
my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
    PrintError => 0,
    RaiseError => 1,
    AutoCommit => 1,
    FetchHashKeyName => 'NAME_lc'
}) or die ("Cannot connect to database $dbfile: $!");

# Export all sequences meeting criteria in Fasta format to a temporary file
my ($fastafh, $fastafile) = tempfile(SUFFIX=>'.fasta',
                                     UNLINK=>1);
my $sth = $dbh->prepare('SELECT pfseq_id, run, seq FROM pfseq WHERE ns <= ? AND len >= ?');
$sth->execute ($maxNs, $minlen);
# If ENA accessions are included then use them in header
if (defined $ENA_acc) {
    while (my @rec = $sth->fetchrow_array) {
        say $fastafh ">gnl|pFrrnENA|$rec[1]_$rec[0]";
        say $fastafh $rec[2];
    }
} else {
    while (my @rec = $sth->fetchrow_array) {
        say $fastafh ">gnl|pFrrn|$rec[0]";
        say $fastafh $rec[2];
    }
}
close $fastafh;

# Makeblastdb
my @cmd = ($makeblastdb,
           '-dbtype nucl',
           '-parse_seqids',
           "-in $fastafile",
           "-out $blastdb",
           "-title=\"$blastdbtitle\"");
system (join " ", @cmd);

# Disconnect from DB
$dbh->disconnect;

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2017 - Brandon Seah <kbseah@mpi-bremen.de>

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