#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:26:19 2018

@author: mac
"""
with open ('ANNOTATIONS/arab_SNP.recode.vcf', 'r') as vcf:
    with open('SNPS.txt', 'w') as SNPS:
        counter = 0
        SNPS.write('"num"\t"global"\t"coords"\t"Chr_snp"\t"from"\t"to"')
        for line in vcf:
            if (line[0] == '#'):
                continue
            else:
                counter += 1
                line_list = line.split('\t')
                if line[5] == ".":
                    continue
                else:
                    #print(str(counter)+'\t'+str(line_list[1])+'\t'+str(line_list[1])+'\t'+str(line_list[0])+'\t'+str(line_list[3])+str(line_list[4]))
                    SNPS.write(str(counter)+'\t'+str(line_list[1])+'\t'+str(line_list[1])+'\t'+str(line_list[0])+'\t'+str(line_list[3])+str(line_list[4]))

        
        
