def Input_func(options,settings):
    for i,x in enumerate(settings):
        x=x.strip()
        if x.startswith("-"):
            option_name = x[1:]
            if option_name == "i":
                options[option_name].append(settings[i+1].strip())
            else:
                options[option_name] = settings[i+1].strip()

    return Check_options(options)

def Check_options(options):
    if options['i'] == []:
        print('No input assigned. quit.')
        quit()

    if options['o'] == None:
        print('No output assigned. quit.')
        quit()

    if options['title'] == None:
        options['title'] = options['i']

    tf = ['show', 'score', 'zeroneg', 'score_correction', 'p_opt']
    for x in tf:
        if options[x] in [False, 'False', 'false', 'No', 'no', '0']:
            options[x] = False
        else:
            options[x] = True

    ints = ['cl', 'skip', 'propose', 'cutoff']
    for x in ints:
        if x in options.keys():
            options[x] = int(options[x])

    floats = ['p', 'm', 'threshold']
    for x in floats:
        if x in options.keys():
            options[x] = float(options[x])

    if options['p_opt']:
        try:
            actives = sum(1.0 for line in open(options["actives"]))
            decoys = sum(1.0 for line in open(options["decoys"]))
            options['p'] = decoys/actives
        except IOError:
            pass

    return options
