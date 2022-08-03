#!/usr/bin/env python3
import sys
from argparse import ArgumentParser
import pickle
import json

from cgelib.output.result import Result
from cgelib.alignment.read_alignment import KMAAlignment

from resfinder.cge.config import Config
from resfinder.cge.resfinder import ResFinder
from resfinder.cge.pointfinder import PointFinder
from resfinder.cge.output.std_results import ResFinderResultHandler
from resfinder.cge.output.std_results import PointFinderResultHandler
from resfinder.cge.kma_manager import KMAManager

#  Modules used to create the extended ResFinder output (phenotype output)
from resfinder.cge.phenotype2genotype.isolate import Isolate
from resfinder.cge.phenotype2genotype.res_profile import PhenoDB
from resfinder.cge.phenotype2genotype.res_sumtable import ResSumTable
from resfinder.cge.phenotype2genotype.res_sumtable import PanelNameError

from resfinder import __version__

# TODO list:
# TODO: JSON output summary is not species dependent


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def handle_results(finder, res, std_result, conf, db_name, method):

    if (method == ResFinder.TYPE_BLAST
        and db_name == 'ResFinder'):

        ResFinderResultHandler.standardize_results_old(std_result,
                                                       res.results,
                                                       db_name)
    elif (method == PointFinder.TYPE_BLAST
            and db_name == 'PointFinder'):
        results_pnt = finder.find_best_seqs(res, conf.pf_gene_cov)

        PointFinderResultHandler.standardize_results_old(std_result,
                                                         results_pnt,
                                                         db_name,
                                                         method)
    else:
        filenames = [x.split('/')[-1] for x in res.alignment_files]
        out_path = \
            [x.rsplit('/', 1)[0] for x in res.alignment_files][0]

        kma_alignment = KMAAlignment(
            output_path=out_path,
            filenames=filenames,
            result_file=["Result", "Alignment"])

        if (method == PointFinder.TYPE_KMA
            and db_name == 'PointFinder'):
            for kmahit in kma_alignment.parse_hits():

                if (float(kmahit['template_coverage']) < conf.pf_gene_cov
                        or float(
                            kmahit['template_identity']) < conf.pf_gene_id):
                    continue

                if conf.specific_gene:
                    if kmahit['templateID'].split('-', 1)[0] != conf.specific_gene:
                        continue

                PointFinderResultHandler.standardize_results_new(std_result,
                                                                 kmahit,
                                                                 db_name,
                                                                 method,
                                                                 finder)
        else: #method == Resfinder.TYPE_KMA and db_name 'ResFinder:
            for kmahit in kma_alignment.parse_hits():
                if (float(kmahit['template_coverage']) < conf.rf_gene_cov
                        or float(kmahit['template_identity']) < conf.rf_gene_id):
                    continue

                ResFinderResultHandler.standardize_results_new(std_result,
                                                               kmahit,
                                                               db_name)

def main():

    # version = read_version("VERSION")
    version = __version__

    ##########################################################################
    # PARSE COMMAND LINE OPTIONS
    ##########################################################################

    parser = ArgumentParser(allow_abbrev=False)

    # General options
    parser.add_argument("-ifa", "--inputfasta",
                        help="Input fasta file.",
                        default=None)
    parser.add_argument("-ifq", "--inputfastq",
                        help="Input fastq file(s). Assumed to be single-end fastq \
                              if only one file is provided, and assumed to be \
                              paired-end data if two files are provided.",
                        nargs="+",
                        default=None)
    parser.add_argument("--nanopore",
                        action="store_true",
                        dest="nanopore",
                        help="If nanopore data is used",
                        default=False)
    parser.add_argument("-o", "--outputPath",
                        help=("Output directory. If it doesn't exist, it will be "
                              "created."),
                        required=True,
                        default=None)
    parser.add_argument("-j", "--out_json",
                        help=("Specify JSON filename and output directory. If the "
                              "directory doesn't exist, it will be created."),
                        default=None)
    parser.add_argument("-b", "--blastPath",
                        help="Path to blastn",
                        default=None)
    parser.add_argument("-k", "--kmaPath",
                        help="Path to KMA",
                        default=None)
    parser.add_argument("-s", "--species",
                        help="Species in the sample",
                        default=None)
    parser.add_argument("--ignore_missing_species",
                        action="store_true",
                        help="If set, species is provided and --point flag is set, "
                             "will not throw an error if no database is found for "
                             "the provided species. If species is not found. Point "
                             "mutations will silently be ignored.",
                        default=False)


    # Acquired resistance options
    parser.add_argument("-db_res", "--db_path_res",
                        help=("Path to the databases for ResFinder. Defaults to "
                              "'db_resfinder' in the ResFinder application "
                              "directory."),
                        default=None)
    parser.add_argument("-db_res_kma", "--db_path_res_kma",
                        help=("Path to the ResFinder databases indexed with KMA. "
                              "Defaults to the value of the --db_res flag."),
                        default=None)
    parser.add_argument("-acq", "--acquired",
                        action="store_true",
                        help="Run resfinder for acquired resistance genes",
                        default=None)
    parser.add_argument("-ao", "--acq_overlap",
                        help="Genes are allowed to overlap this number of\
                              nucleotides. Default: {}.".format(
                            Config.DEFAULT_VALS["acq_overlap"]),
                        type=int,
                        default=None)
    parser.add_argument("-l", "--min_cov",
                        help=("Minimum (breadth-of) coverage of ResFinder within "
                              "the range 0-1."),
                        type=float,
                        default=None)
    parser.add_argument("-t", "--threshold",
                        help=("Threshold for identity of ResFinder within the "
                              "range 0-1."),
                        type=float,
                        default=None)
    # Disinfectant resistance options
    parser.add_argument("-d", "--disinfectant",
                        action="store_true",
                        help="Run resfinder for disinfectant resistance genes",
                        default=False)
    parser.add_argument("-db_disinf", "--db_path_disinf",
                        help=("Path to the databases for DisinFinder."),
                        default=None)
    parser.add_argument("-db_disinf_kma", "--db_path_disinf_kma",
                        help=("Path to the DisinFinder databases indexed with "
                              "KMA. Defaults to the value of the --db_res "
                              "flag."),
                        default=None)

    # Point resistance option
    parser.add_argument("-c", "--point",
                        action="store_true",
                        help="Run pointfinder for chromosomal mutations",
                        default=None)
    parser.add_argument("-db_point", "--db_path_point",
                        help="Path to the databases for PointFinder",
                        default=None)
    parser.add_argument("-db_point_kma", "--db_path_point_kma",
                        help="Path to the PointFinder databases indexed with KMA. \
                              Defaults to the 'kma_indexing' directory inside the \
                              given database directory.",
                        default=None)
    parser.add_argument("-g", "--specific_gene",
                        nargs='+',
                        help="Specify genes existing in the database to \
                              search for - if none is specified all genes are \
                              included in the search.",
                        default=None)
    parser.add_argument("-u", "--unknown_mut",
                        action="store_true",
                        help="Show all mutations found even if in unknown to the\
                              resistance database",
                        default=None)
    parser.add_argument("-l_p", "--min_cov_point",
                        help=("Minimum (breadth-of) coverage of Pointfinder within "
                              "the range 0-1. If None is selected, the minimum "
                              "coverage of ResFinder will be used."),
                        type=float,
                        default=None)
    parser.add_argument("-t_p", "--threshold_point",
                        help=("Threshold for identity of Pointfinder within the "
                              "range 0-1. If None is selected, the minimum "
                              "coverage of ResFinder will be used."),
                        type=float,
                        default=None)
    parser.add_argument("--ignore_indels",
                        action="store_true",
                        help="Ignore frameshift-causing indels in Pointfinder.",
                        default=None)
    parser.add_argument("--ignore_stop_codons",
                        action="store_true",
                        help="Ignore premature stop codons in Pointfinder.",
                        default=None)
    parser.add_argument("-v", "--version", action="version",
                        version=__version__,
                        help="Show program's version number and exit")

    # Temporary option only available temporary
    parser.add_argument("--pickle",
                        action="store_true",
                        help="Create a pickle dump of the Isolate object. \
                              Currently needed in the CGE webserver. Dependency \
                              and this option is being removed.",
                        default=False)

    args = parser.parse_args()

    # Parse and check all arguments and expected files.
    conf = Config(args)

    # Initialise result dict
    std_result = Result.init_software_result(
        name="ResFinder",
        gitdir=f"{conf.resfinder_root}/../../")

    init_result_data = {
        "provided_species": conf.species,
        "software_version": __version__,
        "key": f"ResFinder-{__version__}",
    }
    std_result.add(**init_result_data)

    if(conf.acquired):
        std_result.init_database("ResFinder", conf.db_path_res)
    if(conf.point):
        std_result.init_database("PointFinder", conf.db_path_point_root)
    if(conf.disinf):
        std_result.init_database("DisinFinder", conf.db_path_disinf)

    kma_resultfiles = ["Result", "Fragments"]
    ##########################################################################
    # ResFinder
    ##########################################################################

    if(conf.acquired is True):

        acquired_finder = ResFinder(db_conf_file=conf.db_config_file,
                                    db_path=conf.db_path_res,
                                    pheno_file=conf.phenotype_file,
                                    notes=conf.db_notes_file,
                                    db_path_kma=conf.db_path_res_kma)

        if(conf.inputfasta):
            method = ResFinder.TYPE_BLAST
            resfinder_result = acquired_finder.blast(
                inputfile=conf.inputfasta,
                out_path=conf.outPath_res_blast,
                min_cov=conf.rf_gene_cov,
                threshold=conf.rf_gene_id,
                blast=conf.blast,
                allowed_overlap=conf.rf_overlap)
        else:
            method = ResFinder.TYPE_KMA
            if(conf.nanopore):
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_res_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True,
                                  custom_args="-ont -md 5")
            else:
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_res_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True)

            kma_manager = KMAManager(min_cov=conf.rf_gene_cov,
                                     threshold=conf.rf_gene_id,
                                     databases=acquired_finder.databases,
                                     db_path_kma=conf.db_path_res_kma,
                                     sample_name=conf.sample_name)

            resfinder_result = kma_manager.run_KMAAligner(conf, kma_resultfiles,
                                                          kma_params)

        handle_results(finder=acquired_finder, res=resfinder_result,
                       std_result=std_result, conf=conf, db_name="ResFinder",
                       method=method)

    ##########################################################################
    # DisinFinder
    ##########################################################################
    if(conf.disinf is True):

        # Actually running DisinFinder (for disinfectant resistance)
        disinf_finder = ResFinder(db_conf_file=conf.db_config_disinf_file,
                                  db_path=conf.db_path_disinf,
                                  pheno_file=conf.disinf_file,
                                  notes=conf.db_notes_disinf_file,
                                  db_path_kma=conf.db_path_disinf_kma)

        if(conf.inputfasta):
            method = ResFinder.TYPE_BLAST
            disin_result = disinf_finder.blast(inputfile=conf.inputfasta,
                                                out_path=conf.outPath_disinf_blast,
                                                min_cov=conf.dis_gene_cov,
                                                threshold=conf.dis_gene_id,
                                                blast=conf.blast,
                                                allowed_overlap=conf.dis_overlap)
        else:
            method = ResFinder.TYPE_KMA
            if(conf.nanopore):
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_disinf_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True,
                                  custom_args="-ont -md 5")
            else:
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_disinf_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True)

            kma_manager = KMAManager(min_cov=conf.dis_gene_cov,
                                     threshold=conf.dis_gene_id,
                                     databases=disinf_finder.databases,
                                     db_path_kma=conf.db_path_disinf_kma,
                                     sample_name=conf.sample_name)

            disin_result = kma_manager.run_KMAAligner(conf, kma_resultfiles,
                                                     kma_params)

        handle_results(finder=disinf_finder, res=disin_result,
                       std_result=std_result, conf=conf, db_name="DisinFinder",
                       method=method)

    ##########################################################################
    # PointFinder
    ##########################################################################

    if(conf.point):

        finder = PointFinder(db_path=conf.db_path_point,
                             species=conf.species_dir,
                             gene_list=conf.specific_gene,
                             ignore_indels=conf.ignore_indels,
                             ignore_stop_codons=conf.ignore_stop_codons)

        if(conf.inputfasta):
            method = PointFinder.TYPE_BLAST
            blast_run = finder.blast(inputfile=conf.inputfasta,
                                     out_path=conf.outPath_point_blast,
                                     min_cov=0.01,  # Sorts on coverage later
                                     threshold=conf.pf_gene_id,
                                     blast=conf.blast,
                                     cut_off=False)
            results = blast_run.results

        else:
            method = PointFinder.TYPE_KMA
            if(conf.nanopore):
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_point_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True,
                                  custom_args="-ont -md 5")
            else:
                kma_params = dict(inputfile_1=conf.inputfastq_1,
                                  inputfile_2=conf.inputfastq_2,
                                  output=conf.outPath_point_kma,
                                  cge=True,
                                  apm="p",
                                  kma_1t1=True)

            kma_manager = KMAManager(db_path_kma=conf.db_path_point,
                                     databases=[conf.species_dir],
                                     min_cov=0.01,
                                     threshold=conf.pf_gene_id,
                                     sample_name=conf.sample_name)

            results = kma_manager.run_KMAAligner(conf, kma_resultfiles,
                                                 kma_params)

        handle_results(finder=finder, res=results, std_result=std_result,
                       conf=conf, db_name="PointFinder", method=method)

    ##########################################################################
    # Phenotype to genotype
    ##########################################################################
    # Load genotype to phenotype database
    res_pheno_db = PhenoDB(
        abclassdef_file=conf.abclassdef_file, acquired_file=conf.phenotype_file,
        point_file=conf.point_file, disinf_file=conf.disinf_file,
        disclassdef_file=conf.disclassdef_file)
    # Isolate object store results
    isolate = Isolate(name=conf.sample_name, species=conf.species,
                      amr_panel_file=conf.db_panels_file)

    if(conf.acquired or conf.disinf):
        isolate.load_finder_results(std_table=std_result,
                                    phenodb=res_pheno_db,
                                    type="seq_regions")
    if(conf.point):
        isolate.load_finder_results(std_table=std_result,
                                    phenodb=res_pheno_db,
                                    type="seq_variations")

    isolate.calc_res_profile(res_pheno_db)
    ResFinderResultHandler.load_res_profile(std_result, isolate,
                                            conf.amr_abbreviations)

    if(conf.out_json):
        std_result_file = conf.out_json
    else:
        std_result_file = "{}/{}.json".format(
            conf.outputPath, conf.sample_name.replace("_R1", "").split(".")[0])
    with open(std_result_file, 'w') as fh:
        fh.write(std_result.json_dumps())

    # Create and write the downloadable tab file
    pheno_profile_str = isolate.profile_to_str_table(header=True)

    # TODO: REMOVE THE NEED FOR THE PICKLED FILE
    if(conf.pickle):
        isolate_pickle = open("{}/isolate.p".format(conf.outputPath), "wb")
        pickle.dump(isolate, isolate_pickle, protocol=2)

    pheno_table_file = "{}/pheno_table.txt".format(conf.outputPath)
    with open(pheno_table_file, 'w') as fh:
        fh.write(pheno_profile_str)

    if(conf.species is not None):
        # Apply AMR panel
        input_amr_panels = "{}/phenotype_panels.txt".format(conf.db_path_res)
        res_sum_table = ResSumTable(pheno_profile_str)
        res_sum_table.load_amr_panels(input_amr_panels)

        try:
            panel_profile_str = res_sum_table.get_amr_panel_str(
                panel_name_raw=conf.species, header=True)
        # If specified species does not have an associated panel, just ignore it
        # and exit.
        except PanelNameError:
            eprint("Warning: No panel was detected for the species: {}"
                   .format(conf.species))
            sys.exit()

        amr_panel_filename = conf.species.replace(" ", "_")

        panel_tabel_file = ("{tablename}_{species}.txt"
                            .format(tablename=pheno_table_file[:-4],
                                    species=amr_panel_filename))

        with open(panel_tabel_file, "w") as fh:
            fh.write(panel_profile_str)

    return 0

if __name__ == '__main__':
    sys.exit(main())
