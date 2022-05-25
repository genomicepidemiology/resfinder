#!/usr/bin/env python3

import json
import sys
import os
import pandas as pd
from tabulate import tabulate
from argparse import ArgumentParser

# JSON class
class WriteResultFile:
    """Recreation of all ResFinder result files from ResFinder json file."""
    def __init__(self, json_file, db_path, outdir, sample, mp):
        self.json_file = json_file
        self.db_point_path = db_path
        self.outdir = outdir
        self.sample = sample
        self.mapping_method = mp
        # Load json data and check format
        self.data = WriteResultFile.open_and_read_json(self)
        WriteResultFile.check_json_std_format(self)
        # Load global variable
        self.species = self.data.get("provided_species")
        self.gene_lst= WriteResultFile.point_gene_list(self)

    def parseJSON(self):
        """
        Iterating through seq_regions and seq_variations classes in json file to extract information for result files.
        All Resfinder result files will be written in this function.
        """

        # Extract databases to list
        self.db_list = WriteResultFile.extract_databases(self)
        # Extract phenotypes to dict
        self.phenotype_dict,self.pheno_dict,self._pheno_species_dict = WriteResultFile.read_phenotypes(self)
        
        # Initiate result files
        self.ResFinder_results_tab = open(self.outdir + "/ResFinder_results_tab.txt", "w")
        self.ResFinder_results_tab.write("Resistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\t" +
        "Position in reference\tContig or Depth\tPosition in contig\tPhenotype\tAccession no.\tNotes\n")
        self.ResFinder_disinf_results_tab = open(self.outdir + "/ResFinder_disinf_results_tab.txt", "w")
        self.ResFinder_disinf_results_tab.write("Resistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\t" +
        "Position in reference\tContig or Depth\tPosition in contig\tPhenotype\tAccession no.\tNotes\n")
        self.PointFinder_results_tab = open(self.outdir + "/PointFinder_results.txt", "w")
        self.PointFinder_results_tab.write("Mutation\tIdentity\tNucleotide change\tAmino acid change\tResistance\tPMID\n")
        self.ResFinder_phenotable = open(self.outdir + "/pheno_table.txt", "w")
        self.ResFinder_phenotable.write("# ResFinder phenotype results.\n")
        self.ResFinder_phenotable = open(self.outdir + "/pheno_table.txt", "w")
        self.ResFinder_phenotable.write("# ResFinder phenotype results.\n")
        self.ResFinder_phenotable_species = open(self.outdir + "/pheno_table_" + self.species +".txt", "w")
        self.ResFinder_phenotable_species.write("# ResFinder phenotype results for " + self.species + ".\n")
        
        
        # Iterate through the json seq regions class
        self.res_matrix = []
        self.point_gene_hits = set()
        for j in self.data["seq_regions"].items():

            self.data_row = []
            if j[1].get("ref_database")[0] == self.db_list[0]:     # pointfinder hit
                self.point_gene_hits.add(j[1].get("name"))
            
            elif j[1].get("ref_database")[0] == self.db_list[1] or j[1].get("ref_database")[0] == self.db_list[2]:

                # Add pheno_table data from Resfinder results
                if float(j[1].get("identity")) == 100.0:
                    if j[1].get("alignment_length") == j[1].get("ref_seq_lenght"):
                        for phenotype in j[1].get("phenotypes"):
                            if phenotype in self.pheno_species_dict:
                                self.pheno_species_dict[phenotype][2].append(3)
                                self.pheno_species_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                            self.pheno_dict[phenotype][2].append(3)
                            self.pheno_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                    elif j[1].get("alignment_length") < j[1].get("ref_seq_lenght"):
                        for phenotype in j[1].get("phenotypes"):
                            if phenotype in self.pheno_species_dict:
                                self.pheno_species_dict[phenotype][2].append(1)
                                self.pheno_species_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                            self.pheno_dict[phenotype][2].append(1)
                            self.pheno_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                elif float(j[1].get("identity")) < 100.0:
                    if j[1].get("alignment_length") == j[1].get("ref_seq_lenght"):
                        for phenotype in j[1].get("phenotypes"):
                            if phenotype in self.pheno_species_dict:
                                self.pheno_species_dict[phenotype][2].append(2)
                                self.pheno_species_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                            self.pheno_dict[phenotype][2].append(2)
                            self.pheno_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                    elif j[1].get("alignment_length") < j[1].get("ref_seq_lenght"):
                        for phenotype in j[1].get("phenotypes"):
                            if phenotype in self.pheno_species_dict:
                                self.pheno_species_dict[phenotype][2].append(1)
                                self.pheno_species_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                            self.pheno_dict[phenotype][2].append(1)
                            self.pheno_dict[phenotype][3].append("{} ({})".format(j[1].get("name"),j[1].get("ref_id")))
                
                # Extract informations for tab and txt
                self.data_row.append(j[1].get("name"))
                self.data_row.append(j[1].get("identity"))
                self.data_row.append(str(j[1].get("alignment_length")) + ".."+ str(j[1].get("ref_seq_lenght")))
                self.data_row.append(j[1].get("coverage"))
                self.data_row.append(str(j[1].get("ref_start_pos")) + ".."+ str(j[1].get("ref_end_pos")))
                if self.mapping_method == "kma":
                    self.data_row.append(j[1].get("depth"))
                elif self.mapping_method == "blast":
                    self.data_row.append(j[1].get("query_id"))
                self.data_row.append(str(j[1].get("query_start_pos")) + ".."+ str(j[1].get("query_end_pos")))
                phenotypes = set()
                for pheno in j[1].get("phenotypes"):
                    if pheno in self.phenotype_dict:
                        phenotypes.add(str(self.phenotype_dict.get(j[1].get("phenotypes")[0]))[2:-2].capitalize())
                for i in phenotypes:
                    self.data_row.append(i + " resistance")
                self.data_row.append(j[1].get("ref_acc"))
                self.data_row.append(j[1].get("notes")[0])
                
                if j[1].get("ref_database")[0] == self.db_list[1]:  # resfinder hit
                    self.res_matrix += [self.data_row]
                    # Write Resfinder tab result file
                    WriteResultFile.write_tab(self, self.ResFinder_results_tab)
                elif j[1].get("ref_database")[0] == self.db_list[2]: # disinfectant hit
                    # Write to Disinfectant tab result file
                    WriteResultFile.write_tab(self, self.ResFinder_disinf_results_tab )
                
            else:
                print("Error: unknown reference database has excluded seq region:" + j[0])
        self.ResFinder_results_tab.close()

        # Write pheno_table files
        WriteResultFile.write_pheno_table(self, self.pheno_dict, self.ResFinder_phenotable)
        WriteResultFile.write_pheno_table(self, self.pheno_species_dict, self.ResFinder_phenotable_species)

        # Write Resfinder result txt
        WriteResultFile.write_txt(self)

        # Iterate through the json sequence variations class
        self.point_matrix = []
        self.gene_set = {}
        for j in self.data["seq_variations"].items():

            # Add genes with pointfinder hits
            self.gene_set.add(j[0].split("_",1)[0])

            # Store data if any phenotype
            self.data_row = []
            if j[1].get("phenotypes"):
                self.data_row.append(j[1].get("genes")[0] + j[1].get("seq_var"))
                # RNA mutations
                if j[0].split("_",1)[0] in ("16S","23S"):
                    self.data_row.append(j[1].get("ref_codon").upper() + "->"+ j[1].get("var_codon").upper())
                    self.data_row.append("RNA mutations")
                else:
                    # DNA mutations
                    self.data_row.append(j[1].get("codon_change").upper())
                    self.data_row.append(j[1].get("ref_aa").upper() + "->"+ j[1].get("var_aa").upper())
                self.data_row.append(j[1].get("phenotype").capitalize())
                self.data_row.append(j[1].get("pmids"))
                self.point_matrix += [self.data_row]
                
                # Write PointFinder tab result file
                WriteResultFile.write_tab(self, self.PointFinder_results_tab)
        self.PointFinder_results_tab.close()
        
        # Write PointFinder_table.txt result file
        WriteResultFile.write_point_table(self)

        
    def open_and_read_json(self):
        """Open and load data from json file"""
        self.f = open(self.json_file, 'r')
        data = json.load(self.f)
        self.f.close()
        return data

    def check_json_std_format(self):
        """Checks if all classes in the std ResFinder json format is contained in the json infile."""
        json_file_keys = set(self.data.keys())
        json_std_format_keys = {"type", "databases", "seq_regions", "seq_variations",
        "phenotypes", "software_name", "software_version", "software_commit", "run_date",
        "key", "provided_species" , "result_summary"}
        if not json_file_keys.union(json_std_format_keys):
            sys.exit("Error: json file does not follow the std format.")

    def extract_databases(self):
        """Store reference databases from databases class. Returns a list of database names"""
        self.db_list = []
        for db in self.data.get("databases").keys():
            if db.startswith("P"):      # PointFinder DB
                self.db_list.insert(0,db)
            elif db.startswith("R"):    # ResFinder DB
                self.db_list.insert(1,db)
            elif db.startswith("D"):    # DisinFinder DB
                self.db_list.insert(2,db)
            else:
                return print("Error: Unknown database in json file: " + db)
        return self.db_list

    def read_phenotypes(self):
        """
        Store phenotypes predicted in phneotype class in dicts
        Returns 3 dicts:
            phenotype_dict for Pointfinder and Resfinder hits
            pheno_dict for pheno table overview
            self.pheno_species_dict for species specific pheno table overview
        """
        self.phenotype_dict = {}
        self.pheno_dict = {}
        self.pheno_species_dict = {}
        for j in self.data["phenotypes"].items():

            if j[1].get("amr_resistant"):
                # store data for pheno_table
                if j[1].get("amr_species_relevant"):
                    self.pheno_species_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"Resistant",[],[]]
                    self.pheno_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"Resistant",[],[]]
                else:
                    self.pheno_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"Resistant",[],[]]
                
                if j[1].get("ref_database") == self.db_list[0]:
                    # store phenotype for pointfinder results
                    for k in range(len(j[1].get("seq_variations"))):
                        self.phenotype_dict[j[1].get("seq_variations")[k]] = j[1].get("amr_resistance")
                else:
                    # store phenotype
                    self.phenotype_dict[j[1].get("amr_resistance")] = j[1].get("amr_classes")
            else:
                # store data for pheno_table
                if j[1].get("amr_species_relevant"):
                    self.pheno_species_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"No resistance",[0],[]]
                    self.pheno_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"No resistance",[0],[]]
                else:
                    self.pheno_dict[j[1].get("amr_resistance")] = [j[1].get("amr_classes")[0],"No resistance",[0],[]]

        return self.phenotype_dict, self.pheno_dict, self.pheno_species_dict
    
    def group_matrix(self):
        """Group Pointfinder hits by genes"""
        self.res = {i: [] for i in self.gene_lst}
        for row in self.point_matrix:
            self.res[row[0].split()[0]].append(row)
        return self.res
     
    def write_point_table(self):
        self.PointFinder_table = open(self.outdir + "/PointFinder_table.txt", "w")
        self.PointFinder_table.write("Chromosomal point mutations - Results\nSpecies: {}\nGenes: {}" + 
        "\nMapping methode: {}\n\n\nKnown Mutations\n\n".format(self.species,", ".join(sorted(self.gene_lst)),self.mapping_method))
        res = WriteResultFile.group_matrix(self)
        for gene in res.keys():
            self.PointFinder_table.write(gene + "\n")
            if self.res.get(gene):
                self.PointFinder_table.write("Mutation\tIdentity\tNucleotide change\tAmino acid change\tResistance\tPMID\n")
                WriteResultFile.write_tab(self, self.PointFinder_table)
                self.PointFinder_table.write("\n")
            elif gene in self.point_gene_hits and gene not in self.gene_set:
                self.PointFinder_table.write("No known mutations found in {}\n\n".format(gene))
            else:
                self.PointFinder_table.write("No hit found\n\n")
        self.PointFinder_table.close()

    def write_pheno_table(self, dict, outfile):
        """Write pheno table files from dict."""
        outfile.write("#\n# " + self.sample + "\n" +
        "#The phenotype 'No resistance' should be interpreted with\n" +
        "#caution, as it only means that nothing in the used\n" +
        "# database indicate resistance, but resistance could exist\n" +
        "# from 'unknown' or not yet implemented sources.\n#\n" +
        "# The 'Match' column stores one of the integers 0, 1, 2, 3.\n" +
        "#      0: No match found\n" +
        "#      1: Match < 100% ID AND match length < ref length or Match = 100% ID AND match length < ref length\n" +
        "#      2: Match < 100% ID AND match length = ref length\n" +
        "#      3: Match = 100% ID AND match length = ref length\n" +
        "# If several hits causing the same resistance are found,\n" +
        "# the highest number will be stored in the 'Match' column.\n\n" +
        "# Antimicrobial\tClass\tWGS-predicted phenotype\tMatch\tGenetic background\n")
        for phenotype in dict.keys():
            values = dict.get(phenotype)
            outfile.write("{}\t{}\t{}\t{}\t{}\n".format(phenotype,values[0],values[1],max(values[2]),(", ".join(sorted(values[3])))))
        outfile.close()

    def write_tab(self, outfile):
        """Write a list into a tabseperated file."""
        outfile.write("\t".join(map(str,self.data_row)) + "\n")

    def write_txt(self):
        """Write ResFinder_results.txt from a matrix A=[[]] grouped by phneotypes."""
        df = pd.DataFrame(self.res_matrix, columns=["Resistance gene","Identity",
        "Alignment Length/Gene Length","Coverage","Position in reference",
        "Contig or Depth","Position in contig","Phenotype","Accession no.","Notes"])
        # group by phenotype
        grouped_df = df.sort_values(by=['Phenotype']).groupby(['Phenotype'])

        txt_str = ""
        for key, item in grouped_df:
            title = key
            headers = list(grouped_df.get_group(title).columns)
            rows = list(grouped_df.get_group(title).values)
            table = WriteResultFile.text_table(title, headers, rows)
            txt_str += table

        # print entire table
        with open(self.outdir + "/ResFinder_results.txt", "w") as f:
            f.write(txt_str)

    def text_table(title, headers, rows, table_format='psql'):
        ''' Create text table

        USAGE:
            >>> from tabulate import tabulate
            >>> title = 'My Title'
            >>> headers = ['A','B']
            >>> rows = [[1,2],[3,4]]
            >>> print(text_table(title, headers, rows))
            +-----------+
            | My Title  |
            +-----+-----+
            |    A |    B |
            +=====+=====+
            |    1 |    2 |
            |    3 |    4 |
            +-----+-----+
        '''
        # Create table
        table = tabulate(rows, headers, tablefmt=table_format)
        # Prepare title injection
        width = len(table.split('\n')[0])
        tlen = len(title)
        if tlen + 4 > width:
            # Truncate oversized titles
            tlen = width - 4
            title = title[:tlen]
        spaces = width - 2 - tlen
        left_spacer = ' ' * int(spaces / 2)
        right_spacer = ' ' * (spaces - len(left_spacer))
        # Update table with title
        table = '\n'.join(['+%s+' % ('-' * (width - 2)),
                          '|%s%s%s|' % (left_spacer,
                          title, right_spacer),
                          table, '\n'])
        return table
    
    def point_gene_list(self):
        """Extract genes from resistance overview in Pointfinder 
        for relevant species. A list with genes is returned"""
        gene_set = set()
        file = os.path.join(self.db_point_path, self.species, "resistens-overview.txt")
        with open(file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                gene = line.split()[0].split("_", maxsplit=1)[0]
                gene_set.add(gene)
        self.gene_lst = list(gene_set)
        return self.gene_lst


##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-i", "--infile",
                    help="JSON file from ResFinder",
                    required=True)
parser.add_argument("-o", "--outdir",
                    help="Path to output directory",
                    default=".")
parser.add_argument("-db_point", "--db_path_point",
                    help=("Path to the databases for PointFinder."),
                    required=True, default=None)
parser.add_argument("-mp", "--mapping_method",
                    help="Mapping method used for ResFinder run, either KMA or BLAST",
                    required=True)
parser.add_argument("-s", "--sample_name",
                    help="Name of sample file ran with ResFinder",
                    required=True)

args = parser.parse_args()

# Load and check input file(s)
if args.infile is None:
    sys.exit("Input Error: No input file provided!\n")
infile = args.infile
if not os.path.exists(infile):
    sys.exit("Input Error: Input file not found at expected location: %s"%(infile))

# Check if valid output directory is provided
if args.outdir:
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
else:
    outdir = os.getcwd()

# Check if valid database is provided
if args.db_path_point:
    db_path = os.path.abspath(args.db_path_point)
    if not os.path.exists(args.db_path_point):
       sys.exit("Input Error: The specified database directory does not "
                "exist:" + db_path + "\n") 
# Check existence of the resistens-overview file in a species
res_overview_file = "%s/campylobacter/resistens-overview.txt" % (db_path)
if not os.path.isfile(res_overview_file):
    sys.exit("Resistance-overview.txt not found at expected location: %s"%(res_overview_file))


# Defining varibales
method_path = args.mapping_method
sample = args.sample_name


# Call JSON class to recreate ResFinder result files
results = WriteResultFile(infile, db_path, outdir, sample, method_path)
results.check_json_std_format()
results.parseJSON()
