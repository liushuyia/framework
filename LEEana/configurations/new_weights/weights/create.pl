#!/usr/bin/perl
open(infile,"list.dat");

while(<infile>){
    my $temp1 = $_;
    chomp($temp1);
    @temp = split(/\s+/,$temp1);
    my $dir = $temp[0];
    open(infile1,"type.dat");
    while(<infile1>){
	my $type = $_;
	chomp($type);
	print "./processed_checkout_rootfiles/checkout_$dir\.root   /data0/wgu/1019_checkout_eventweight_sep24/$dir\/$type\.root       ./processed_checkout_rootfiles/$dir\/$type\.root  $type\n";
    }
    close(infile1);
}
print "end    end     end     end\n";
