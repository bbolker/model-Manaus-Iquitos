pomp <- function (data, t0,
                  particles, rprocess, dmeasure, rmeasure,
                  ...) {
  new(
      'pomp',
      data = data,
      t0 = t0,
      particles = particles,
      rprocess = rprocess,
      dmeasure = dmeasure,
      rmeasure = rmeasure,
      userdata = list(...)
      )
}

