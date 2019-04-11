#! /usr/bin/perl
#perl get_gene_pairs.pl ~/result/adjp_cut.txt ~/result/pcc_cut.txt ~/result/gene_pairs.txt
open I1, "<$ARGV[0]";
open I2, "<$ARGV[1]";
open OUT, ">$ARGV[2]";
my %list;
while(<I1>){
    chomp;
    @data0 = split/\t/,$_;
    if(/^coding(.*)/){}
    else{
		$list{$data0[0].$data0[1]} = $data0[2];
	}
}
while(<I2>){
    chomp;
	@data = split/\t/,$_;
	if(/^coding(.*)/){
		print OUT "mRNA\tlncRNA\tpcc\tadjp\n";
	}
	if(exists $list{$data[0].$data[1]}){
		$adjp = $list{$data[0].$data[1]};
		print OUT "$_\t$adjp\n";
	}
}		
close I1;
close I2;
close OUT;

