#!/usr/bin/perl
use strict;
use warnings;
use DBI;
my $chrSql = $ARGV[0]; 
my $driver   = "SQLite";
my $dsn = "DBI:$driver:dbname=$chrSql";
my $userid = "";
my $password = "";
my $db_chr = DBI->connect($dsn, $userid, $password, { RaiseError => 1, AutoCommit => 1 })
                      or die $DBI::errstr;

my $sth = $db_chr->table_info('', 'main', '%', 'TABLE');
my @table_chrs;
while(my $r = $sth->fetchrow_hashref){
	push @table_chrs, $r->{TABLE_NAME}; 
}
foreach my$table_chr(@table_chrs){
	my $sth_db_chr = $db_chr->prepare("select CHR,COORD1,NT1,NT2,SCORE from $table_chr where SCORE!='' AND SCORE >= 0");
	my @rows;
	$sth_db_chr->execute();
	while (@rows = $sth_db_chr->fetchrow_array()){
			$rows[0]=~ s/chr//;
			print"$rows[0]\t$rows[1]\t$rows[2]\t$rows[3]\t$rows[4]\n";
	}
}
$db_chr->disconnect(); 
