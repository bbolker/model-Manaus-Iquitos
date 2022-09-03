setClass(
         'pomp',
         representation(
                        data = 'matrix',
                        t0 = 'numeric',
                        particles = 'function',
                        rprocess = 'function',
                        dmeasure = 'function',
                        rmeasure = 'function',
                        userdata = 'list'
                        )
         )

setClass(
         'mif',
         representation(
                        'pomp',
                        ivps = 'character',
                        pars = 'character',
                        stvs = 'character',
                        Nmif = 'integer',
                        alg.pars = 'list',
                        coef = 'numeric',
                        random.walk.sd = 'numeric',
                        pred.mean = 'matrix',
                        pred.var = 'matrix',
                        filter.mean = 'matrix',
                        conv.rec = 'matrix',
                        cond.loglik = 'numeric',
                        eff.sample.size = 'numeric',
                        loglik = 'numeric'
                        )
         )

