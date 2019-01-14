#!/usr/bin/env perl

# Set up SQLite database to store phyloFlash results

=head1 NAME

B<phyloFlash_sqliteDB_setup.pl> - Initialize SQLite database for phyloFlash results

=head1 SYNOPSIS

phyloFlash_sqliteDB_setup.pl -db [database filename]

phyloFlash_sqliteDB_setup.pl -help

=head1 DESCRIPTION

Initialize SQLite database for phyloFlash results. 

After initializing the database, use B<phyloFlash_sqliteDB_add.pl> to populate the
database with data from actual phyloFlash results.

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::RealBin/.."; # Relative path to PhyloFlash.pm
use PhyloFlash qw(err msg @msg_log);
#use XML::LibXML;
use DBI;
use Getopt::Long;
use Pod::Usage;

my $dbfile = 'test.sqlite';  # Name for DB file

pod2usage(-verbose=>0) if !@ARGV;

GetOptions ("dbfile=s" => \$dbfile,
            "help" => sub { pod2usage(-verbose=>1); },
            "man" => sub { pod2usage(-verbose=>2); },
            ) or pod2usage("Invalid options");

=head1 OPTIONS

=over 8

=item --db I<FILE>

Name of SQLite database file to be created.

Default: I<test.sqlite>

=item --help

Short help message

=item --man

Full manual page

=back

=cut

# Check that filenames are supplied
die ("No SQLite db filename specified") if (!defined $dbfile);

# Initialize SQLite database
msg "Initializing SQLite database file";
my $dsn = "dbi:SQLite:dbname=$dbfile";
my $user = "";
my $password = "";
my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
    PrintError => 0,
    RaiseError => 1,
    AutoCommit => 1,
    FetchHashKeyName => 'NAME_lc'
});

# Create table "pfrun"
msg "Creating tables for SQLite database";
my $pfrun_table_sql = << 'END_SQL';
CREATE TABLE pfrun (
    libname                         VARCHAR(25),
    run                             VARCHAR(20),
    pfversion                       VARCHAR(50),
    input_reads                     INTEGER,
    mapped_SSU_reads                INTEGER,
    mapped_SSU_read_segments        INTEGER,
    assembled_SSU_read_segments     INTEGER,
    unassembled_SSU_read_segments   INTEGER,
    added                           DATE,
    PRIMARY KEY (libname, run)
)
END_SQL

$dbh->do($pfrun_table_sql);

# Create table "pfseq"
my $pfseq_table_sql = << 'END_SQL';
CREATE TABLE pfseq (
    pfseq_id                  VARCHAR(50),
    libname                   VARCHAR(25),
    run                       VARCHAR(20),
    tool                      VARCHAR(10),
    seq                       VARCHAR(2500),
    len                       INTEGER,
    ns                        INTEGER,
    read_cov                  INTEGER,
    coverage                  REAL,
    dbhit                     VARCHAR(20),
    taxonomy                  VARCHAR(350),
    pid                       REAL,
    alnlen                    INTEGER,
    uchime_ref_score          NUMERIC,
    uchime_ref_parent_A       TEXT,
    uchime_ref_parent_B       TEXT,
    uchime_ref_idQM           NUMERIC,
    uchime_ref_idQA           NUMERIC,
    uchime_ref_idQB           NUMERIC,
    uchime_ref_idAB           NUMERIC,
    uchime_ref_idQT           NUMERIC,
    uchime_ref_LY             INTEGER,
    uchime_ref_LN             INTEGER,
    uchime_ref_LA             INTEGER,
    uchime_ref_RY             INTEGER,
    uchime_ref_RN             INTEGER,
    uchime_ref_RA             INTEGER,
    uchime_ref_div            NUMERIC,
    uchime_ref_YN             TEXT,
    uchime_ref_db             TEXT,
    PRIMARY KEY (pfseq_id, libname, run)
)
END_SQL

$dbh->do($pfseq_table_sql);

# Disconnect DB file
msg "Disconnecting database, setup complete";
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