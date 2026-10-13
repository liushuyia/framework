#!/usr/bin/perl

my $tot_num = 1;

open(infile,"type.dat");

my $num = 1;
while(<infile>){
    my $input = $_;
    chomp($input);

    open(infile1,"list.dat");

    while(<infile1>){
	my $input1 = $_;
	chomp($input1);
	@temp = split(/\s+/);
	print "$temp[4]    $temp[3]     $num   ./processed_checkout_rootfiles/$temp[0]\/$input\.root          ./hist_rootfiles/XsFlux/cov_$num\.root       0         $tot_num       $temp[2]      $temp[1]\n";
	$tot_num ++;
    }
    close(infile1);
    $num ++;
}
print "-1        end        -1        end             end         -1            -1           -1-        -1\n";
