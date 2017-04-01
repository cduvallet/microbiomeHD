"""
Code written by Thomas Gurry, thomasgurry@gmail.com

OVERVIEW:

Python module for creating, updating, and interpreting individual dataset summary files and objects.  These summary files are tab-delimited, standard format files that are intended to standardize and facilitate interactions between the user, the data, and different computing requests associated with the data.  They contain they information about the dataset in addition to the location of dataset-specific files.

Data from the summary file can be loaded into a Python object, modified, and the changes written back to file.  In this manner, information about the dataset can be efficiently accessed and manipulated.


A summary file object is created in Python as follows:

obj = SummaryParser('summary_file.txt')


An example summary file is found between the dotted lines:

=============================================
DATASET_ID	obio_testdata

#16S_start
DATASET_ID	obio_testdata
RAW_FASTQ_FILE	obio.raw.fastq
PRIMERS_FILE	obio.primers.lst
BARCODES_MAP	obio.barcodes.lst
BARCODES_MODE   1
PROCESSED	True
ASCII_ENCODING  ASCII_64
AMPLICON	V1-V2
KEYWORDS	fmt
#16S_end

=============================================

"""

import numpy as np
import os
import os.path
import sys

class SummaryParser():
    def __init__(self, summary_file):
        self.summary_file = summary_file
        self.datasetID = None

        # Initialize 16S attributes
        self.attribute_value_16S = {'PROCESSED': "N/A"}

        # Initialize ITS attributes
        self.attribute_value_ITS = {'PROCESSED': "N/A"}

        # Check summary file for integrity
        self.SummaryFileChecker()

    def SummaryFileChecker(self):
        # Checks summary file for minimal entries and appropriate format
        with open(self.summary_file, 'r') as summary_fid:
            all_lines = summary_fid.readlines()
        # Check for contents
        if len(all_lines) == 0:
            print "Summary file appears to be empty.  Check its contents before proceeding."
            raise NameError("Summary file appears to be empty.  Check its contents before proceeding.")

        try:
            # Check for 16S lines
            checksum = 0
            for line in all_lines:
                if line.split('\t')[0].rstrip('\n\r') == "#16S_start":
                    checksum += 1
                elif line.split('\t')[0].rstrip('\n\r') == "#16S_end":
                    checksum += 1
            if checksum != 2:
                raise NameError("No 16S lines found in summary file.")
            # Check for tab delimitation and empty space characters
            for line in all_lines:
                if len(line.rstrip('\n\r')) > 0 and line[0] != "#":   # If not one of the empty lines or lines defining sections, check for tab delimitation
                    if len(line.split(' ')) > 1:
                        print "Empty space characters detected in summary file.  Please make tab-delimited."
                        raise NameError("Empty space characters detected in summary file.  Please make tab-delimited.")
                    if len(line.split('\t')) == 1:
                        print "No tab characters found in summary file attribute line '" + line + "'.  Please make tab-delimited."
                        raise NameError("No tab characters found in summary file attribute line '" + line + "'.  Please make tab-delimited.")
            # Check for minimal contents for 16S processing.
            checksum = 0
            for line in all_lines:
                if line.split('\t')[0] == "DATASET_ID" and len(line.split('\t')[1]) > 0:
                    checksum += 1
                if line.split('\t')[0] == "BARCODES_MAP" and len(line.split('\t')[1]) > 0:
                    checksum += 1
                if line.split('\t')[0] == "PRIMERS_FILE" and len(line.split('\t')[1]) > 0:
                    checksum += 1
                if (line.split('\t')[0] == "RAW_FASTQ_FILE" or line.split('\t')[0] == "RAW_FASTA_FILE" or line.split('\t')[0] == "RAW_FASTQ_FILES" or line.split('\t')[0] == "RAW_FASTA_FILES") and len(line.split('\t')[1]) > 0:
                    checksum += 1
            if checksum != 4:
                print "Failed to find DATASET_ID, BARCODES_MAP, PRIMERS_FILE and/or RAW_FASTQ_FILE/RAW_FASTA_FILE/RAW_FASTQ_FILES lines.  Check integrity of summary file."
                raise NameError("Failed to find DATASET_ID, BARCODES_MAP, PRIMERS_FILE and/or RAW_FASTQ_FILE/RAW_FASTA_FILE/RAW_FASTQ_FILES lines.  Check integrity of summary file.")
        except:
            try:
                # Check for ITS lines
                checksum = 0
                for line in all_lines:
                    if line.split('\t')[0].rstrip('\n\r') == "#ITS_start":
                        checksum += 1
                    elif line.split('\t')[0].rstrip('\n\r') == "#ITS_end":
                        checksum += 1
                if checksum != 2:
                    raise NameError("No ITS lines found in summary file.")
                # Check for tab delimitation and empty space characters
                for line in all_lines:
                    if len(line.rstrip('\n\r')) > 0 and line[0] != "#":   # If not one of the empty lines or lines defining sections, check for tab delimitation
                        if len(line.split(' ')) > 1:
                            raise NameError("Empty space characters detected in summary file.  Please make tab-delimited.")
                        if len(line.split('\t')) == 1:
                            raise NameError("No tab characters found in summary file attribute line '" + line + "'.  Please make tab-delimited.")
                # Check for minimal contents for ITS processing.
                checksum = 0
                for line in all_lines:
                    if line.split('\t')[0] == "DATASET_ID" and len(line.split('\t')[1]) > 0:
                        checksum += 1
                    if line.split('\t')[0] == "BARCODES_MAP" and len(line.split('\t')[1]) > 0:
                        checksum += 1
                    if line.split('\t')[0] == "PRIMERS_FILE" and len(line.split('\t')[1]) > 0:
                        checksum += 1
                    if (line.split('\t')[0] == "RAW_FASTQ_FILE" or line.split('\t')[0] == "RAW_FASTA_FILE" or line.split('\t')[0] == "RAW_FASTQ_FILES" or line.split('\t')[0] == "RAW_FASTA_FILES") and len(line.split('\t')[1]) > 0:
                        checksum += 1
                if checksum != 4:
                    raise NameError("Failed to find DATASET_ID, BARCODES_MAP, PRIMERS_FILE and/or RAW_FASTQ_FILE/RAW_FASTA_FILE/RAW_FASTQ_FILES lines.  Check integrity of summary file.")
            except:
                raise NameError("No 16S or ITS sections found.")
                # WGS should go here.


    def Extract16SLines(self):
        # Description:  Extracts the line numbers pertaining to 16S
        with open(self.summary_file,'r') as summary_fid:
            all_lines = summary_fid.readlines()
            for i in range(len(all_lines)):
                line = all_lines[i].split()
                if(len(line)>0):
                    if(line[0] == "#16S_start"):
                        startline = i+1
                    if(line[0] == "#16S_end"):
                        endline = i
                        break
        return [startline, endline]


    def ExtractITSLines(self):
        # Description:  Extracts the line numbers pertaining to ITS
        with open(self.summary_file,'r') as summary_fid:
            all_lines = summary_fid.readlines()
            for i in range(len(all_lines)):
                line = all_lines[i].split()
                if(len(line)>0):
                    if(line[0] == "#ITS_start"):
                        startline = i+1
                    if(line[0] == "#ITS_end"):
                        endline = i
                        break
        return [startline, endline]


    def ReadSummaryFile(self):
        # Description:  Loads the summary file specified in the object.
        with open(self.summary_file,'r') as summary_fid:
            summary_file_lines = summary_fid.readlines()
            # Read dataset ID
            self.datasetID = summary_file_lines[0].split('\t')[1].rstrip('\n\r')

            # Read 16S attributes, if any.
            try:
                [startline, endline] = self.Extract16SLines()
                for i in range(startline, endline):
                    line = summary_file_lines[i].split('\t')
                    attribute = line[0]
                    value = line[1].rstrip('\n\r')
                    self.attribute_value_16S[attribute] = value
            except:
                pass

            # Read ITS attributes, if any.
            try:
                [startline, endline] = self.ExtractITSLines()
                for i in range(startline, endline):
                    line = summary_file_lines[i].split('\t')
                    attribute = line[0]
                    value = line[1].rstrip('\n\r')
                    self.attribute_value_ITS[attribute] = value
            except:
                pass


    def WriteSummaryFile(self):
        # Description :  Writes a tab-delimited summary file using the values currently in the summary object.
        with open(self.summary_file,'w') as summary_fid:

            # Populate summary file
            summary_fid.write('DATASET_ID' + '\t' + self.datasetID + '\n')
            summary_fid.write('\n')

            # 16S portion
            summary_fid.write('#16S_start' + '\n')
            for [key, value] in self.attribute_value_16S.items():
                summary_fid.write(key + '\t' + value + '\n')
            summary_fid.write('#16S_end' + '\n')


            # ITS portion
            summary_fid.write('#ITS_start' + '\n')
            for [key, value] in self.attribute_value_ITS.items():
                summary_fid.write(key + '\t' + value + '\n')
            summary_fid.write('#ITS_end' + '\n')
