from os.path import splitext
import logging
import argparse


def Input_func():

    import sys

    parser = argparse.ArgumentParser(prog="SIEVE-Score_v1.1",
                                     description="Calculate SIEVE-Score",
                                     fromfile_prefix_chars='@')

    parser.add_argument("-i", "--input", action="append", required=True,
                        help="input Glide files")
    parser.add_argument("-o", "--output", nargs="?",
                        default="SIEVE-Score.csv", help="output csv file")
    parser.add_argument("-l", "--log", nargs="?", default="SIEVE-Score.log",
                        help="log file")
    parser.add_argument("--hits", help="hits file")

    parser.add_argument("-t", "--title", default="",
                        help="job name")
    parser.add_argument("-n", "--num_cluster", type=int, default=10,
                        help="number of clusters")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="debug mode")
    parser.add_argument("-a", "--annotate", action="store_true",
                        help="annotate score/name onto result plot for debug")
    parser.add_argument("-z", "--zeroneg", action="store_true",
                        help="If set, compounds which are not in the hits " +
                             "file are treated as negative, " +
                             "for evaluation set.")
    parser.add_argument("-p", "--plus", type=float, default=1.0,
                        help="plus score of SIEVE-Score_1.0")
    parser.add_argument("-m", "--minus", type=float, default=-1.0,
                        help="minus score of SIEVE-Score_1.0")
    parser.add_argument("--opt_coef", nargs=2, default=False,
                        help="optimize score coefficient.\n" +
                             "require 2 arguments of file or int: " +
                             "actives, decoys")
    parser.add_argument("--propose", type=int, default=1000,
                        help="max number of output compounds")
    parser.add_argument("--noclustering", dest="clustering",
                        action="store_false",
                        help="If set, clustering will not used.")
    parser.add_argument("--score_type", default="normal", 
                        choices=["normal", "exp", "Gauss"],
                        help="score function type")
    parser.add_argument("--dist", default="euclidean",
                        choices=["euclidean", "mahalanobis"],
                        help="distance metric of score function #not implemented")
    parser.add_argument("--score_cutoff", dest="cutoff", type=float,
                        default=1.0, help="cutoff_value")
    parser.add_argument("--score_dim", dest="dim", type=float, default=1.0,
                        help="change parameter of score function")
    parser.add_argument("--score_scale", dest="scale", type=float, default=1.0,
                        help="change parameter of score function")
    parser.add_argument("--score_sigma", dest="sigma", type=float, default=1.0,
                        help="change parameter of score function")
    parser.add_argument("--nprocs", type=int, default=1,
                        help="number of cpus to use, default=1")

    args = parser.parse_args(sys.argv[1:])
    args = Check_options(args)
    return args


def set_log_info(args):
    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(format='%(asctime)s: %(levelname)s:%(message)s',
                        level=log_level, filename=args.log, filemode="w")

    logger = logging.getLogger(__name__)
    return logger


def Check_options(args):

    logger = set_log_info(args)

    if args.opt_coef is not False:
        # opt_coef = [actives, decoys]
        actives = num_molecule(args.opt_coef[0])
        decoys = num_molecule(args.opt_coef[1])
        args.p = decoys / actives * float()

    return args


def num_molecule(x):
    if x.isdigit():
        return int(x)
    else:
        try:
            from schrodinger import structure
        except ImportError:
            logger.exception("error in p_optimize. " +
                             "if you want to count molecules of file, " +
                             "please run in schrod env.",
                             exc_info=True)
            quit()

        st = structure.StructureReader(x)
        return sum(1 for _ in st)
