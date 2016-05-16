from os.path import splitext
import logging

def Input_func(options,settings):
    for i,x in enumerate(settings):
        if x.startswith("-"):
            option_name = x[1:]
            if option_name == "i":
                options[option_name].append(settings[i+1].strip())
            else:
                options[option_name] = settings[i+1].strip()

    return Check_options(options)

def set_log_info(options):
    if options.has_key("debug") and options["debug"]:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(format='%(asctime)s: %(levelname)s:%(message)s',
                        level=log_level, filename=options["log"],filemode="w")

    logger = logging.getLogger(__name__)
    return logger


def Check_options(options):
    tf = ['show', 'score', 'zeroneg', 'score_correction', 'p_opt', 'debug']
    for x in tf:
        if x in options.keys():
            if options[x] in [False, 'False', 'false', 'No', 'no', '0']:
                options[x] = False
            else:
                options[x] = True

    ints = ['cl', 'skip', 'propose']
    for x in ints:
        if x in options.keys():
            options[x] = int(options[x])

    floats = ['p', 'm', 'threshold', 'opt_coef','cutoff','score_dim',
              'scale', 'sigma']
    for x in floats:
        if x in options.keys():
            options[x] = float(options[x])

    if options['log'] is None:
        if options['o'] is None:
            print("parameters are not set.")
            quit()
        else:
            options['log'] = splitext(options['o'])[0]+'.log'

    logger = set_log_info(options)

    if options['i'] == []:
        logger.error('No input assigned.')
        quit()

    if options['o'] is None:
        logger.error('No output assigned.')
        quit()

    if 'title' not in options.keys():
        options['title'] = options['i']

    if options['p_opt']:
        if (options["actives"] is None or
            options["decoys"]  is None):
            logger.error("error in p_optimize, actives/decoys is not set")
            quit()
        try:
            from schrodinger import structure
            st = structure.StructureReader(options["actives"])
            actives = sum(1.0 for _ in st)
            st = structure.StructureReader(options["decoys"])
            decoys  = sum(1.0 for _ in st)
            options['p'] = decoys/actives * float(options['opt_coef'])

        except:
            logger.exception("error in p_optimize. please run in schrod env.",
                             exc_info=True)
            quit()

    return options
