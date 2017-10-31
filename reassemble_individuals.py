#!/usr/bin/python
"""
This takes the output from Msatallele (the MLG 1 col per locus file)
and attempts to "re-construct" each individual, taking into account where a
locus has been genotyped more than once in an individual.

It assumes that an individual is spread across several rows. This is likely if
you have repeat data or have used Fragman for your peak scoring

The first thing it will do is flag up any inconsistencies in your data
Look at these more closely as fix your input table

Once all inconsistencies are removed it will give you 2 outputs.
One which contains just one row per individual
The second which contains count of how many times a locus was repeatedly genotyped
and you got a consistent result

This should help with calculating error rates and the like

written Sep 2017
Graeme Fox
g.fox@mmu.ac.uk
"""

import argparse, re
parser = argparse.ArgumentParser(description='file_converter_for_msatallele_output')
parser.add_argument('-i','--input1', help='Msatallele output - the MLG 1 col per locus', \
                    required=True)
parser.add_argument('-o','--output1', help='Genotype Data Outputt', \
                    required=True)
parser.add_argument('-c','--counts', help='Count Data Outputt', \
                    required=True)
args = parser.parse_args()

dicts = {}

with open (args.output1, 'w') as output:
    with open (args.counts, 'w') as count_output:
    #    output.write("line\tPanel\tMarker\tSize_1\tSize_2\n")
        with open(args.input1, 'r') as input:
            first_line = input.readline()
            ### get list of markers
            markers = []
            for x in first_line.split("\t")[1:]:
                markers.append(x.rstrip("\n"))

            unique_samples = set()
            # get number of markers
            marker_count = len(first_line.split("\t")[1:])
            print("Found " + str(marker_count) + " markers.")
            for line in input:
                ### get the sample name
                # do some find and replace for inconsistencies in sample name
                line = line.replace("--", "_")
                line = line.replace("LizPoint", "LizardPoint")
                line = line.replace("Lizard_Point", "LizardPoint")
                line = line.replace("Liz_Point", "LizardPoint")
                line = line.replace("Lizard_point", "LizardPoint")
                line = line.replace("MIS", "MidIrishSea_")
                line = line.replace("DON", "Donegal_")
                line = line.replace("FOF", "FirthofForth_")
                line = line.replace("IOS", "IslesofScilly_")
                line = line.replace("County_Clare", "CountyClare")
                line = line.replace("F005", "005")
                line = line.replace("Wexford_144", "Wexford_1044")
                line = line.replace("Wexford_1043", "Waterford_1043")
                line = line.replace("Wexford_1046", "Waterford_1046")
                line = line.replace("Wexford_1047", "Waterford_1047")
                line = line.replace("Wexford_1048", "Waterford_1048")
                line = line.replace("Wexford_1049", "Waterford_1049")

                #chop every pre-decing the sample name off the front
                num = re.sub(r'data_[0-9_]*(mplex[0-9]*_)+', "", line)
                num2 = re.sub(r'_[A-H0-9]*_[0-9]*.fsa', "", num)
                num3 = re.sub(r'_mplex[0-9]*', "", num2)
                sample_number = int(num3.split("\t")[0].split("_")[1])
                sample_name = (num3.split("\t")[0].split("_")[0] + "_" + "%04d" % sample_number)
                unique_samples.add(sample_name)

                # if this sample does not currently have a dict, create it
                if not sample_name in dicts:
                    dicts[sample_name] = {}
                for i, marker in zip(range(0,marker_count), markers):
                    # if that sample/marker combination contains an allele
                    if str(line.split("\t")[i+1].rstrip("\n")) != "":
                        # it this marker already has an entry for this sample
                        if marker in dicts[sample_name]:
                            # but it curretnly holds an NA
                            if dicts[sample_name][marker.rstrip("\n")] == "NA":
                                # replace it with the proper data you've found
                                dicts[sample_name][marker.rstrip("\n")] = str(line.split("\t")[i+1].rstrip("\n"))
                                dicts[sample_name][marker + "_count_data"] = 1
                            ## if it already exists and is already holding some data
                            else:
                                #check if it is the same value as the one you have just found
                                if dicts[sample_name][marker.rstrip("\n")] == str(line.split("\t")[i+1].rstrip("\n")):
                                    dicts[sample_name][marker + "_count_data"] = dicts[sample_name][marker + "_count_data"] + 1
                                    # and if it is, do nothing
                                    #print("success. You are consistent.")
                                    continue
                                else:
                                    # if it is a repeat genotype, but there is inconcistency between the two results
                                    # flag it as being something to double check
                                    print("Multiple genotypes found for sample: " \
                                            + sample_name + " at marker: " \
                                            + marker)
                                    continue
                        # if this sample/marker does not currently have an entry
                        if not marker in dicts[sample_name]:
                        # if we have a genotype to fill in
                        # create a marker entry in the sample dict and put in the genotype
                            dicts[sample_name][marker.rstrip("\n")] = str(line.split("\t")[i+1].rstrip("\n"))
                            # as this is a new entry in the dictionary, we can also create its count_data key and assign it as zero
                            dicts[sample_name][marker + "_count_data"] = 1
                    # if there is no genotype available, create the marker and fill with NA
                    else:
                        if not marker in dicts[sample_name]:
                            dicts[sample_name][marker.rstrip("\n")] = "NA"
                            dicts[sample_name][marker.rstrip("\n") + "_count_data"] = "0"
        header = "Sample"
        structure_header = ""
        for x in markers:
            header = header + "\t" + x
            structure_header = structure_header + "\t " + x
        output.write(header + "\n")
        count_output.write(header + "\n")

        for x in unique_samples:
            genotypeoutput = str(x) + "\t"
            countoutput = str(x)
            for y in markers:
                genotypeoutput = genotypeoutput + str(dicts[x][y] + "\t")
                countoutput = countoutput + "\t" + str(dicts[x][y + "_count_data"])

            output.write(genotypeoutput + "\n")
            count_output.write(countoutput + "\n")


        ### give it all the lat, long info of the different sites
        lat_long = {}
        location_codes = {}
        lat_long["FirthofForth"] = {}
        lat_long["FirthofForth"]["latitude"] = str(56.130284)
        lat_long["FirthofForth"]["longtitude"] = str(-2.763474)
        lat_long["FirthofForth"]["code"] = 1
        lat_long["Wexford"] = {}
        lat_long["Wexford"]["latitude"] = str(52.316916)
        lat_long["Wexford"]["longtitude"] = str(-6.325791)
        lat_long["Wexford"]["code"] = 2
        lat_long["Donegal"] = {}
        lat_long["Donegal"]["latitude"] = str(54.555871)
        lat_long["Donegal"]["longtitude"] = str(-8.377615)
        lat_long["Donegal"]["code"] = 3
        lat_long["Portreath"] = {}
        lat_long["Portreath"]["latitude"] = str(50.286217)
        lat_long["Portreath"]["longtitude"] = str(-5.308941)
        lat_long["Portreath"]["code"] = 4
        lat_long["CountyClare"] = {}
        lat_long["CountyClare"]["latitude"] = str(52.895294)
        lat_long["CountyClare"]["longtitude"] = str(-9.532182)
        lat_long["CountyClare"]["code"] = 5
        lat_long["LizardPoint"] = {}
        lat_long["LizardPoint"]["latitude"]= str(49.956634)
        lat_long["LizardPoint"]["longtitude"]= str(-5.205915)
        lat_long["LizardPoint"]["code"] = 6
        lat_long["Waterford"] = {}
        lat_long["Waterford"]["latitude"] = str(52.169763)
        lat_long["Waterford"]["longtitude"] = str(-6.936047)
        lat_long["Waterford"]["code"] = 2
        lat_long["MidIrishSea"] = {}
        lat_long["MidIrishSea"]["latitude"] = str(53.332399)
        lat_long["MidIrishSea"]["longtitude"] = str(-5.302501)
        lat_long["MidIrishSea"]["code"] = 7
        lat_long["IslesofScilly"] = {}
        lat_long["IslesofScilly"]["latitude"] = str(49.935254)
        lat_long["IslesofScilly"]["longtitude"] = str(-6.321319)
        lat_long["IslesofScilly"]["code"] = 8



        ## write genepop format
        with open(args.output1.rstrip(".txt") + "_genepop.txt", 'w') as genepop:
            with open(args.output1.rstrip(".txt") + "_popgenreport.txt", 'w') as popgenreport:
                with open(args.output1.rstrip(".txt") + "_structure.txt", 'w') as structure_output:
                    structure_output.write(structure_header.strip("\t") + "\n")
                    genepop.write("A descriptive title\n")
                    header = header.replace("\t", ",")
                    header = header.replace("Sample,", "")
                    genepop.write(header + "\n")
                    popgenreport.write("ind,pop,lat,long," + header + "\n")
                    genepop.write("POP\n")
                    for x in unique_samples:
                        pop = x.split("_")[0]
                        genotypeoutput = str(x) + ",\t"

                        popgenoutput = str(x) + "," + pop + "," + lat_long[pop]["latitude"] \
                                    + "," + lat_long[pop]["longtitude"] + ","

                        structure_output_line = str(x) + "\t" + str(lat_long[pop]["code"]) + "\t"

                        countoutput = str(x)

                        for y in markers:
                            genotypeoutput = genotypeoutput + str(dicts[x][y].replace("NA","000000") + "\t")
                            allele = str(dicts[x][y].replace("NA","000000") + "\t")
                            popgenoutput = popgenoutput + allele[0:3] + "/" + allele[3:6] + ","
                            popgenoutput = popgenoutput.replace("000/000", "NA")
                            structure_output_line = structure_output_line + allele[0:3] + "\t" + allele[3:6] + "\t"
                        popgenreport.write(popgenoutput.rstrip(",") + "\n")
                        genepop.write(genotypeoutput.rstrip(",") + "\n")
                        structure_output.write(structure_output_line + "\n")
