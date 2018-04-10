#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:26:19 2018

@author: mac
"""

with open ('1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf', 'r') as vcf:
    with open('SNPS.txt', 'w') as SNPS:
        counter = 0
        SNPS.write('"num"\t"global"\t"coords"\t"Chr_snp"\t"from"\t"to"')
        for line in vcf:
            if (line[0:1] == '##') or (line[0:2] == "POS"):
                continue
            else:
                counter += 1
                line_list = line.split('\t')
                if line[5] == ".":
                    continue
                else:
                    SNPS.write(counter+'\t'+line_list[1]+'\t'+line_list[1]+'\t'+line_list[0]+'\t'+line_list[3]+line_list[4])
            
        
        
