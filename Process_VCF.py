#!/usr/bin/env python
# coding: utf-8

#------------------------------
# Import the needed libraries 
#------------------------------
import pandas as pd
import numpy as np
import csv, os, sys, re
import logging
import argparse

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)






'''

this file contains a VCF class to help parse VCF files 

can be used by sourcing this file in python:
    exec(open("path_to_script/Process_VCF.py").read())

then run:
    vcf = VCF("path_to_VCF_file")
    vcf.configure_vcf()
    vcf.parseINFO()
    vcf.parseEncoding()
    vcf.fullVCF()
If you have edited a vcf and want to print the new vcf in the correct vcf format,
then run the following.
    vcf.writeVCF(output = <outout path file name>)


'''
def unique(list_obj):
    unique_list = []
    for i in list_obj:
        if i not in unique_list:
            unique_list.append(i)
    return unique_list


class VCF:
    '''
    Define the VCF object 
    '''
    def __init__(self, vcf_path):
        '''
        Create properties for the new object and assign values to them
            Parse the arguments given to the script 
            should include: The VCF file of interest 
                            HapMap file destination 
                            Output file location 
        '''
        
        #***********************************
        # create and assign the attributes 
        #***********************************
        self.vcf_path = vcf_path 
        


        #------------------------------------------
        # Read in the VCF file as a pandas dataframe and store it in the object
        #------------------------------------------
        logger.info("Reading VCF File")
        self.vcf = pd.read_csv(self.vcf_path,
                                    sep='\t', low_memory=False, comment='#', header =None, index_col = False,
                                    names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"])



    #********************************************
    # Additional methods for the Adaboost Object 
    #********************************************
    
    def configure_vcf(self):
        '''
        Function to read in and configure the VCF files
            Configure the VCF and separate the lines into a dictionary by chromosome 
            Return the VCF header and dictionary
        '''
        logger.info("Processing the VCF")
        vcf_header = []
        vcf_dic = {}

        # Check if the file was gunzipped and treat it appropriately
        if self.vcf_path.endswith('.gz'):
            vcf = gzip.open( self.vcf_path, mode='rt' )
        else:
            vcf = open( self.vcf_path, "r" )

        for line in csv.reader(vcf, delimiter = "\t"):
            # check to see if header, is so append the line to header 
            if line[0][0] == "#":
                vcf_header.append(line)
                continue
            # else append the line to the dictionary based on the chromosome
            chromosome = line[0]
            position = line[1]
            new_key = str(chromosome) + ":" + str(position)
            vcf_dic[new_key]= line
        vcf.close()

        # Add the vcf header and body to the new object 
        self.vcf_header = vcf_header
        self.vcf_body = vcf_dic
        return( self )


    def parseINFO(self):
        '''
        Parse the vcf INFO Column 
            Load in the data VCF then convert INFO column into pandas dataframe 
        '''

        logger.info("Processing the input data.")
        
        ## subset the data to get the get 'Chr', 'Pos','REF','ALT'
        df_vcf_subset = self.vcf[['CHROM', 'POS','REF','ALT']]


        #---------------------------------------------------------------------
        # Load in the data VCF and convert info column into pandas dataframe 
        #---------------------------------------------------------------------
        ## Read in the data header of the vcf file to get info column ID's
        ## Separate the info column of the vcf for each variant 
        ## create a pandas dataframe for the information column 


        # lists to hold the ID's and the values 
        info_id = []
        all_info = []
        # Get the ID's from the header
        for i in self.vcf_header:
            if i[0][0:11] == "##INFO=<ID=":
                info_id.append(str(i[0].split(",")[0][11:])) 
        # print(info_id)

        # Iterate through each variant
        for i in self.vcf_body:
            info_num = [None]*len(info_id)
            ## split the info section
            info = self.vcf_body[i][7].split(";")
            if "" in info:
                info.remove("")
            for value in info:
                ## pull out the ID and value 'ID=Value'
                temp = value.split("=")
                ## If the ID has no value (given by IndexError), make binary indicator, 1=pressent 
                try: 
                    info_num[info_id.index(temp[0])] = temp[1]
                except IndexError:
                    info_num[info_id.index(temp[0])] = 1
            all_info.append(info_num)
        df_info = pd.DataFrame(data = all_info)
        df_info.columns = info_id
        # print(df_info.head())
        self.info = df_info

        return self


    def parseEncoding(self):
        '''
        Parse the ENCODING column in the VCF
        '''
        logger.info("Parsing the encoding.")
        
        # Get all possible names 
        FORMAT_row_lists = list(self.vcf.FORMAT.str.split(":"))
        total = [j for i in FORMAT_row_lists for j in i]
        column_names = unique(total)

        # Parse the encodings 
        encoding = self.vcf.ENCODING.str.split(":")
        # go over each row and combine the FORMAT with the ENCODING
        all_rows = []
        for i in range(len(FORMAT_row_lists)):
            # make the dictionary and append to list 
            a = dict(zip(FORMAT_row_lists[i], encoding[i]))
            all_rows.append(a)
        # convert the list dictionaries into a dataframe 
        encoding_df = pd.DataFrame(all_rows)
        
        self.encoding = encoding_df
        
        return self


    def fullVCF(self):
        '''
        Adding the parsed info columnn
        '''
        logger.info("Adding the parsed info and encoding sections.")
        
        # check
        # make sure each variant has its associated info 
        # checkDF(self.vcf, self.info)

        # combine the two data frames
        self.vcf_full = pd.concat([self.vcf, self.info], axis=1)

        # check
        # make sure each variant has its associated ENCODING 
        # checkDF(self.vcf_full, encoding)

        # combine the two data frames
        self.vcf_full = pd.concat([self.vcf_full, self.encoding], axis=1)
        
        return self

    

    def writeVCF(self, output):
        '''
        Option to write the vcf to a new outout vcf file.
        '''
        logger.info("Writing new VCF file.")
        
        #~~~~~~~~~~~~~
        # Header
        #~~~~~~~~~~~~~
        # Open the new file and write the header 
        header = self.vcf_header
        with open(output,'w') as csv_file:
            ## Format the header into a string
            header[-1] = ["\t".join(header[-1])]
            for i in header:
                i.insert(len(i),"\n")
            header_string = "".join(sum(header,[]))
            ## Write the header 
            csv_file.write(header_string)

        #~~~~~~~~~~~~~
        # Body
        #~~~~~~~~~~~~~
        # append the variants to the output file 
        self.vcf.to_csv(output, sep='\t', index=False, mode = "a", header=None, quoting=csv.QUOTE_NONE)




###################################
# option to edit the contig label 

# def configEdit(vcf, )