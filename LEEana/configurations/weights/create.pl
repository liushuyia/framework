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
	#print "./processed_checkout_rootfiles/checkout_$dir\.root   /data0/xqian/MicroBooNE/numi/$dir\/$type\.root       ./processed_checkout_rootfiles/$dir\/$type\.root  $type\n";
	#print "./processed_checkout_rootfiles/$dir\.root   /data0/xqian/MicroBooNE/new_checkout_files/$dir\/$type\.root       ./processed_checkout_rootfiles/$dir\/$type\.root  $type\n";
	#print "./processed_checkout_rootfiles/checkout_$dir\.root   /data0/wgu/checkout_links/$dir\/$type\.root       ./processed_checkout_rootfiles/$dir\/$type\.root  $type\n";
	print "./processed_checkout_rootfiles/checkout_$dir\.root   /data0/xqian/MicroBooNE/nc_files/$dir\/$type\.root       ./processed_checkout_rootfiles/$dir\/$type\.root  $type\n";
    }
    close(infile1);
}
print "end    end     end     end\n";
