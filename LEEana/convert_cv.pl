#!/usr/bin/perl

open(infile,"filelist_cv");
my $num = 0;
while(<infile>){
    $filename = $_;
    chomp($filename);
    @temp = split(/\s+/,$filename);

    if ($temp[0] eq "end") {
	last;
    }

    if ($num % 12 == 11){
	system("bin/convert_cv_spec $temp[2]$temp[1] $temp[3]$temp[1] ");
    }else{
	system("bin/convert_cv_spec $temp[2]$temp[1] $temp[3]$temp[1] &");
    }
    
    $num++;
}

