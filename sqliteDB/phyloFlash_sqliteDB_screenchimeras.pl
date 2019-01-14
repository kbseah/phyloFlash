#!/usr/bin/env perl

# Sort SPAdes, EMIRGE, or MATAM sequences from a single pF run by coverage
# into a FASTA file, and run Vsearch Uchime_ref

=head1 NAME

B<phyloFlash_sqliteDB_screen_chimeras.pl> - Screen sequences in pFrrn database for chimeras

=head1 SYNOPSIS

phyloFlash_sqliteDB_screen_chimeras.pl -db [SQLite db] -num [num seq to screen] -refdb [Vsearch db] -CPU [num CPU]

phyloFlash_sqliteDB_screen_chimeras.pl -man

=head1 DESCRIPTION

Screen reconstructed sequences in the pfseq table of a pFrrn database for
sequence chimerism, by comparing vs a database with the uchime_ref algorithm,
as implemented in Vsearch 2.4.0+. Output will be parsed and saved to database.

=cut

use strict;
use warnings;
use 5.010;

use DBI;
use Getopt::Long;
use File::Temp qw(tempfile);
use Pod::Usage;
use FindBin;

my $dbfile;
my $tool = 'all';
my $num = 100;
my $CPU = 8;
my $vsearch_db = "$FindBin::RealBin/../132/SILVA_SSU.noLSU.masked.trimmed.fasta";
my $vsearch_path = 'vsearch';
my %seqhash;
my @seqarr;

pod2usage(2) if (!@ARGV);

GetOptions ('db=s' => \$dbfile,
            'num=i' => \$num,
            'tool=s' => \$tool,
            'refdb=s' => \$vsearch_db,
            'CPU=i' => \$CPU,
            'help' => sub { pod2usage(1) },
            'man' => sub { pod2usage(-exitval=>0, -verbose=>2) },
           ) or pod2usage(2);

=head1 INPUT ARGUMENTS

=over 8

=item -db FILE

Path to pFrrn SQLite database file. Default: none

=item -num INTEGER

Max number of sequences to process. Script checks for records in the pfseq
table which have empty uchime_ref_score field. Default: 100

=item -tool STRING

Screen only sequences produced by this tool. Valid values: "all", "spades",
"emirge", "matam". Default: "all"

=item -refdb FILE

Path to reference database (Fasta format) used by Vsearch uchime_ref to screen
for chimeras. Default: ../132/SILVA_SSU.noLSU.masked.trimmed.fasta

=item -CPU INTEGER

Number of processors used by Vsearch. Default: 2

=item -help

Short help message

=item -man

Full manual page

=back

=cut

# Die if mandatory options not given
die ('No database specified') if !defined $dbfile;
die ('No number of sequences specified') if !defined $num;
die ('No tool specified') if !defined $tool;

# Connect to database
say STDERR "Connecting to database $dbfile";
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $user = "";
my $password = "";
my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
    PrintError => 0,
    RaiseError => 1,
    AutoCommit => 0,
    FetchHashKeyName => 'NAME_lc'
}) or die ("ERROR: Cannot connect to database $dbfile: $!");

my $sth = $dbh->prepare('SELECT pfseq_id, libname, run, read_cov, tool, seq FROM pfseq WHERE uchime_ref_score IS NULL LIMIT ?');
$sth->execute($num);
while (my $row_href = $sth->fetchrow_hashref) {
    push @seqarr, $row_href if $tool eq 'all' || $tool eq $row_href->{'tool'};
}

# Create temp Fasta file for sequences 
my ($fastafh,$fastafile) = tempfile(SUFFIX=>".fasta",
                                    UNLINK=>1
                                    );
say STDERR "Temporary fasta file $fastafile";
# Write sequences to temp file
for (my $i=0; $i <= $#seqarr; $i++) {
    say $fastafh ">$i";
    say $fastafh $seqarr[$i]->{'seq'};
}

# Create temp file for vsearch uchime_ref output
my ($uchimeoutfh, $uchimeoutfile) = tempfile (SUFFIX=>".uchimeout",
                                              UNLINK=>1
                                              );
say STDERR "UCHIME output written to $uchimeoutfile";

# Construct uchime ref command
my @vsearch_args = ("--uchime_ref $fastafile",
                    "--db $vsearch_db",
                    "--uchimeout $uchimeoutfile",
                    "--threads $CPU");
my $vsearch_cmd = join " ", ($vsearch_path, @vsearch_args);
my $vsearch_return = system ($vsearch_cmd);

say STDERR "Vsearch completed with return value $vsearch_return";

# Read results into DB
open (my $uchime_fhin, "<", $uchimeoutfile) or die ("$!");
my $counter = 0;
while (my $line = <$uchime_fhin>) {
    $counter ++;
    chomp $line;
    my @split = split /\t/, $line;
    my @curr_key = ($seqarr[$split[1]]->{'pfseq_id'},
               $seqarr[$split[1]]->{'libname'},
               $seqarr[$split[1]]->{'run'}
              );
    
    #say STDERR "Updating record with Uchime results for sequence $split[1]";
    $dbh->do('UPDATE pfseq SET
                uchime_ref_score = ?,
                uchime_ref_parent_A = ?,
                uchime_ref_parent_B = ?,
                uchime_ref_idQM = ?,
                uchime_ref_idQA = ?,
                uchime_ref_idQB = ?,
                uchime_ref_idAB = ?,
                uchime_ref_idQT = ?,
                uchime_ref_LY = ?,
                uchime_ref_LN = ?,
                uchime_ref_LA = ?,
                uchime_ref_RY = ?,
                uchime_ref_RN = ?,
                uchime_ref_RA = ?,
                uchime_ref_div = ?,
                uchime_ref_YN = ?,
                uchime_ref_db = ?
              WHERE pfseq_id = ? AND libname = ? AND run = ?',
            undef,
            $split[0],
            @split[2, 3, 5 .. 17],
            $vsearch_db,
            @curr_key);
    if ($counter%100 == 0) {
        $dbh->commit or die $dbh->errstr; # Commit every 100 sequences
        say STDERR "Added $counter sequences to database";
    }
}
$dbh->commit or die $dbh->errstr; # Commit any leftover transactions
say STDERR "Added $counter sequences to database";
close ($uchime_fhin);

# Disconnect from db
say STDERR "Disconnecting from database";
$sth->finish;
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
