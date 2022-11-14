from os.path import splitext
import logging
import argparse


def Input_func():
    import sys

    parser = argparse.ArgumentParser(prog="SIEVE-Score_v1.1",
                                     description="Calculate SIEVE-Score",
                                     fromfile_prefix_chars='@')

    parser.add_argument("-i", "--input", required=True,
                        help="input Glide pv file or .interation file")
    parser.add_argument("-o", "--output", nargs="?",
                        default="SIEVE-Score.csv", help="output csv file")
    parser.add_argument("-l", "--log", nargs="?", default="SIEVE-Score.log",
                        help="log file")
    parser.add_argument("--hits", help="hits file")
    
    parser.add_argument("-m", "--mode", default="cv",
                        choices=["cv", "comparesvm", "paramsearch",
                                 "datasize", "interaction", "screen"],
                        help=
                        """Mode. Available patterns: 
                        cv (default): performance evaluation.
                        comparesvm: compare SVM and RF.
                        paramsearch: perform parameter search.
                        datasize: perform dataset size analysis.
                        screen: do screening using independent training and test.
                        interacion: interaction energy vector extraction only.
                        """)
    parser.add_argument("-t", "--title", default="",
                        help="job name")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="debug mode")
    parser.add_argument("--annotate", action="store_true",
                        help="annotate score/name onto result plot for debug")
    parser.add_argument("-z", "--zeroneg", action="store_true",
                        help="If set, compounds which are not in the hits " +
                             "file are treated as negative, " +
                             "for evaluation set.")
    parser.add_argument("--active", default=None,
                        help="Number of active compounds.\n" + 
                             "Only for calculate EF. Type: file or int")
    parser.add_argument("--decoy", default=None,
                        help="Number of decoy compounds.\n" + 
                             "Only for calculate EF. Type: file or int")
    parser.add_argument("--propose", type=int, default=100000,
                        help="max number of output compounds")
    parser.add_argument("--nprocs", type=int, default=1,
                        help="number of cpus to use, default=1")
    parser.add_argument("--inter_only", action="store_true",
                        help="If set, only output interaction file.")
    parser.add_argument("--use_docking_score", action="store_true",
                        help="If set, add docking score for feature.")
    parser.add_argument("--model", type=str, default="RF",
                        help="Model type, RF or SVM.")
    parser.add_argument("--reverse", action="store_true",
                        help="Reverse the final order.")
    parser.add_argument("--n_splits", type=int, default=5,
                        help="datasize only. param for StratifiedShuffleSplit.")
    parser.add_argument("--random_state", type=int, default=None,
                        help="datasize only. param for StratifiedShuffleSplit.")
    parser.add_argument("--train_size", default="0.9",
                        help="datasize only. param for StratifiedShuffleSplit.")
    parser.add_argument("--n_iter", type=int, default=1,
                        help="eval only. number of iteration.")
    parser.add_argument("--testdata", default=None,
                        help="screen only. test compound file.")
    parser.add_argument("--testhits", default=None,
                        help="screen only. test compound hits file if available.")
    parser.add_argument("--ignore", default=None,
                        help="ignore compounds file to train/test.")
    args = parser.parse_args(sys.argv[1:])
    logger = set_log_info(args)

    if args.active is not None:
        args.active = num_molecule(args.active)
    if args.decoy is not None:
        args.decoy = num_molecule(args.decoy)

    try:
        args.train_size = int(args.train_size)
    except ValueError:
        args.train_size = float(args.train_size)
    
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


def num_molecule(x):
    if x.isdigit():
        return int(x)
    else:
        try:
            from schrodinger import structure
        except ImportError:
            logger.exception("if you want to count number of molecules in file, " +
                             "please run in schrodinger env.",
                             exc_info=True)
            quit()

        st = structure.StructureReader(x)
        return sum(1 for _ in st)
