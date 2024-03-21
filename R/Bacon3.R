Bacon3 <- function (core = "MSB2K", thick = 5, coredir = "", prob = 0.95, 
                    d.min = NA, d.max = NA, add.bottom = TRUE, d.by = 1, seed = NA, 
                    depths.file = FALSE, depths = c(), depth.unit = "cm", age.unit = "yr", 
                    unit = depth.unit, acc.shape = 1.5, acc.mean = 20, mem.strength = 10, 
                    mem.mean = 0.5, boundary = NA, hiatus.depths = NA, hiatus.max = 10000, 
                    add = c(), after = 1e-04/thick, cc = 1, cc1 = "IntCal20", 
                    cc2 = "Marine20", cc3 = "SHCal20", cc4 = "ConstCal", cc.dir = c(), 
                    postbomb = 0, delta.R = 0, delta.STD = 0, t.a = 3, t.b = 4, 
                    normal = FALSE, suggest = TRUE, accept.suggestions = FALSE, 
                    reswarn = c(10, 200), remember = TRUE, ask = TRUE, run = TRUE, 
                    defaults = "defaultBacon_settings.txt", sep = ",", dec = ".", 
                    runname = "", slump = c(), remove = FALSE, BCAD = FALSE, 
                    ssize = 4000, th0 = c(), burnin = min(500, ssize), youngest.age = c(), 
                    oldest.age = c(), MinAge = c(), MaxAge = c(), cutoff = 0.01, 
                    plot.pdf = TRUE, dark = 1, date.res = 100, age.res = 200, 
                    yr.res = age.res, close.connections = TRUE, older.than = c(), 
                    younger.than = c(), save.elbowages = FALSE, verbose = TRUE,
                    suppress.plots = TRUE, bacon.change.thick = FALSE,
                    ...) 
{
  # Temporarily attach all internal functions from package rbacon
  attach(loadNamespace("rbacon"), name = "rbacon_all", warn.conflicts = FALSE)
  attach(loadNamespace("rintcal"), name = "rintcal_all", warn.conflicts = FALSE)
  
  
  coredir <- assign_coredir(coredir, core, ask, isPlum = FALSE)
  if (core == "MSB2K" || core == "RLGH3") {
    if (!dir.exists(file.path(coredir, core))) {
      dir.create(file.path(coredir, core), showWarnings = FALSE, 
                 recursive = TRUE)
      fileCopy <- system.file(paste0("extdata/Cores/", 
                                     core), package = "rbacon")
      file.copy(fileCopy, coredir, recursive = TRUE, overwrite = FALSE)
    }
  }
  if (length(cc.dir) == 0) 
    cc.dir <- system.file("extdata", package = "rintcal")
  cc.dir <- validateDirectoryName(cc.dir)
  defaults <- system.file("extdata", defaults, package = packageName())
  dets <- read.dets(core, coredir, sep = sep, dec = dec, cc = cc)
  if (ncol(dets) > 4 && length(cc) > 0) {
    cc.csv <- unique(dets[, 5])
    if (verbose) {
      if (length(cc.csv) == 1) {
        if (cc.csv != cc) 
          message(" Using calibration curve specified within the .csv file,", 
                  cc.csv, "\n")
      }
      else if (min(cc.csv) == 0) 
        message(" Using a mix of cal BP and calibrated C-14 dates\n")
      else message(" Using several C-14 calibration curves\n")
    }
  }
  if (suggest) {
    sugg <- sapply(c(1, 2, 5), function(x) x * 10^(-1:2))
    ballpacc <- lm(dets[, 2] * 1.1 ~ dets[, 4])$coefficients[2]
    ballpacc <- abs(sugg - ballpacc)
    ballpacc <- ballpacc[ballpacc > 0]
    sugg <- sugg[order(ballpacc)[1]]
    if (!sugg %in% acc.mean) 
      if (accept.suggestions) {
        acc.mean <- sugg
        message("Adapting acc.mean to ", sugg, " ", 
                age.unit, "/", depth.unit)
      }
    else {
      ans <- readline(message(" Ballpark estimates suggest changing the prior for acc.mean to ", 
                              sugg, " ", age.unit, "/", depth.unit, ". OK? (y/N) "))
      if (tolower(substr(ans, 1, 1)) == "y") 
        if (bacon.change.thick) {
          message(" Setting thick to ", sugg, "\n")
          thick <- sugg
          acc.mean <- sugg
        }
        
      else message(" No problem, using the provided prior")
    }
  }
  if (thick < d.by) 
    warning("Please set d.by to a value smaller than that of thick", 
            .call = TRUE)
  if (mem.mean < 0 || mem.mean > 1) 
    stop("The prior for the mean of the memory should be between 0 and 1", 
         call. = FALSE)
  if (length(mem.mean) > 1) 
    stop("Can only use one value for mem.mean across a core", 
         call. = FALSE)
  if (length(mem.strength) > 1) 
    stop("Can only use one value for mem.strength across a core", 
         call. = FALSE)
  if (!is.na(boundary[1])) 
    boundary <- sort(unique(boundary))
  if (!is.na(hiatus.depths[1])) {
    hiatus.depths <- sort(unique(hiatus.depths))
    if (length(acc.mean) == 1) 
      acc.mean <- rep(acc.mean, length(hiatus.depths) + 
                        1)
  }
  info <- Bacon.settings(core = core, coredir = coredir, dets = dets, 
                         thick = thick, remember = remember, d.min = d.min, d.max = d.max, 
                         d.by = d.by, depths.file = depths.file, slump = slump, 
                         acc.mean = acc.mean, acc.shape = acc.shape, mem.mean = mem.mean, 
                         mem.strength = mem.strength, boundary = boundary, hiatus.depths = hiatus.depths, 
                         hiatus.max = hiatus.max, BCAD = BCAD, cc = cc, postbomb = postbomb, 
                         cc1 = cc1, cc2 = cc2, cc3 = cc3, cc4 = cc4, depth.unit = depth.unit, 
                         normal = normal, t.a = t.a, t.b = t.b, delta.R = delta.R, 
                         delta.STD = delta.STD, prob = prob, defaults = defaults, 
                         runname = runname, ssize = ssize, dark = dark, youngest.age = youngest.age, 
                         oldest.age = oldest.age, cutoff = cutoff, age.res = age.res, 
                         after = after, age.unit = age.unit)
  assign_to_global("info", info)
  info$coredir <- coredir
  if (is.na(seed)) 
    seed <- sample(1:1e+06, 1)
  set.seed(seed)
  info$seed <- seed
  info$isplum <- FALSE
  if (length(MinAge) > 0) 
    warning("please do not use MinAge, instead use the option youngest.age", 
            call. = FALSE)
  if (length(MaxAge) > 0) 
    warning("please do not use MaxAge, instead use the option oldest.age", 
            call. = FALSE)
  if (info$t.b - info$t.a != 1) 
    warning("t.b - t.a should always be 1, check the manual", 
            call. = FALSE)
  if (min(acc.shape) < 1) 
    warning("\nWarning, using values <1 for acc.shape might cause unexpected results\n", 
            call. = TRUE)
  negativeages <- FALSE
  if (info$postbomb == 0) 
    if (ncol(info$dets) == 4) {
      if (info$cc > 0) 
        if (min(info$dets[, 2]) < 0) 
          negativeages <- TRUE
    }
  else if (max(info$dets[, 5]) > 0) 
    if (min(info$dets[which(info$dets[, 5] > 0), 2]) < 
        0) 
      negativeages <- TRUE
  if (negativeages) 
    stop("you have negative C14 ages so should select a postbomb curve", 
         call. = FALSE)
  info$calib <- bacon.calib(dets, info, date.res, cc.dir = cc.dir, 
                            cutoff = cutoff)
  info$rng <- c()
  for (i in 1:length(info$calib$probs)) {
    tmp <- info$calib$probs[[i]]
    info$rng <- range(info$rng, tmp[, 1])
  }
  if (length(th0) == 0) 
    info$th0 <- round(rnorm(2, max(youngest.age, dets[1, 
                                                      2]), dets[1, 3]))
  info$th0[info$th0 < info$youngest.age] <- info$youngest.age
  if (length(depths) == 0) 
    depths <- seq(info$d.min, info$d.max, by = d.by)
  if (depths.file) {
    dfile <- paste0(info$coredir, info$core, "/", info$core, 
                    "_depths.txt")
    if (!file.exists(dfile)) 
      stop("I cannot find the file ", paste0(info$coredir, 
                                             info$core, "/", info$core, "_depths.txt"), call. = FALSE)
    depths <- fastread(dfile, header = FALSE)[, 1]
    if (!is.numeric(depths[1])) 
      stop("File should contain numbers only, no headers", 
           call. = FALSE)
  }
  info$depths <- depths
  if (min(depths) < info$d.min) 
    info$d.min <- min(depths)
  if (max(depths) > info$d.max) 
    info$d.max <- max(depths)
  info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), 
                     by = thick)
  if (add.bottom) 
    info$elbows <- c(info$elbows, max(info$elbows) + thick)
  info$K <- length(info$elbows)
  info$cK <- info$d.min + (info$thick * info$K)
  if (!is.na(info$hiatus.depths[1]) || !is.na(info$boundary[1])) {
    ifelse(is.na(info$boundary[1]), hd <- info$hiatus.depths, 
           hd <- info$boundary)
    if (min(hd) < info$d.min) 
      stop("cannot have hiatus above the core's top depth. Adapt hiatus.depths or d.min.", 
           call. = FALSE)
    if (max(hd) + info$thick > info$d.max) 
      stop("the age-depth model should have at least one section below the one containing the deepest hiatus. Adapt thick or d.max?", 
           call. = FALSE)
    if (length(hd) > 1) {
      above <- c()
      for (i in hd) above <- c(above, max(which(info$elbows <= 
                                                  i)))
      if (any(diff(above) < 2)) 
        stop("we need at least 2 section elbows between hiatuses. Choose fewer hiatuses, different depths, more sections (decrease thick) or a different d.min.\n ", 
             call. = FALSE)
    }
  }
  ans <- "n"
  if (suggest) 
    if (length(reswarn) == 2) 
      if (info$K < min(reswarn)) {
        sugg <- pretty(thick * (info$K/min(reswarn)), 
                       10)
        sugg <- min(sugg[sugg > 0])
        if (accept.suggestions) 
          ans <- "y"
        else ans <- readline(message(" Warning, the current value for thick, ", 
                                     thick, ", will result in very few age-model sections (", 
                                     info$K, ", not very flexible). Suggested maximum value for thick: ", 
                                     sugg, " OK? (y/n) "))
      }
  else if (info$K > max(reswarn)) {
    sugg <- max(pretty(thick * (info$K/max(reswarn))))
    if (accept.suggestions) 
      ans <- "y"
    else ans <- readline(message(" Warning, the current value for thick, ", 
                                 thick, ", will result in very many age-model sections (", 
                                 info$K, ", possibly hard to run). Suggested minimum value for thick: ", 
                                 sugg, " OK? (y/n) "))
  }
  if (tolower(substr(ans, 1, 1)) == "y") {
    message(" Setting thick to ", sugg, "\n")
    thick <- sugg
    info$thick <- thick
    info$elbows <- seq(floor(info$d.min), ceiling(info$d.max), 
                       by = thick)
    if (length(info$slump) > 0) 
      info$elbows <- seq(floor(info$d.min), toslump(ceiling(info$d.max), 
                                                    info$slump, remove = remove), by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
  }
  if (length(slump) > 0) {
    if (length(slump)%%2 == 1) 
      stop("slumps need both upper and lower depths. Please check the manual", 
           call. = FALSE)
    slump <- matrix(sort(slump), ncol = 2, byrow = TRUE)
    info$slump <- slump
    slumpdmax <- toslump(ceiling(info$d.max), slump, remove = remove)
    info$elbows <- seq(floor(info$d.min), slumpdmax, by = thick)
    info$K <- length(info$elbows)
    info$cK <- info$d.min + (info$thick * info$K)
    info$slumpfree <- toslump(depths, slump, remove = remove)
    info$slumphiatus <- toslump(info$hiatus.depths, slump, 
                                remove = remove)
    if (!is.na(info$boundary[1])) {
      info$slumpboundary <- toslump(info$boundary, slump, 
                                    remove = remove)
      info$slumphiatus <- info$slumpboundary
    }
    slumpdets <- info$dets
    slumpdets[, 4] <- toslump(slumpdets[, 4], slump, remove = remove)
    info$slumpdets <- slumpdets[!is.na(slumpdets[, 4]), 
    ]
  }
  info$prefix <- paste0(coredir, core, "/", core, runname, 
                        "_", info$K)
  info$coredir <- coredir
  info$bacon.file <- paste0(info$prefix, ".bacon")
  if (!file.exists(outfile <- paste0(info$prefix, ".out"))) 
    file.create(outfile)
  if (file.mtime(outfile) < file.mtime(paste0(info$coredir, 
                                              core, "/", core, ".csv"))) 
    message("Warning! The file with the dates seems newer than the run you are loading. If any dates have been added/changed/removed?, then please run Bacon.cleanup()")
  if (BCAD) 
    info$BCAD <- TRUE
  if (!is.na(boundary[1])) {
    if (length(slump) > 0) 
      boundary <- info$slumpboundary
    info$hiatus.depths <- boundary
    if (length(add) == 0) 
      add <- max(1, 1.5 * max(info$acc.mean))
    info$hiatus.max <- add
  }
  assign_to_global("info", info)
  prepare <- function() {
    if (suppress.plots == FALSE){
    pn <- c(1, 2, 3, 3)
    if (!is.na(info$hiatus.depths[1])) 
      if (is.na(info$boundary[1])) 
        pn <- c(1, 2, 3, 4, 4, 4)
    layout(matrix(pn, nrow = 2, byrow = TRUE), heights = c(0.3, 
                                                           0.7))
    oldpar <- par(mar = c(3, 3, 1, 1), mgp = c(1.5, 0.7, 
                                               0), bty = "l", yaxs = "i")
    on.exit(par(oldpar))
    PlotAccPrior(info$acc.shape, info$acc.mean, depth.unit = depth.unit, 
                 age.unit = age.unit)
    PlotMemPrior(info$mem.strength, info$mem.mean, thick)
    if (!is.na(info$hiatus.depths)[1]) 
      if (is.na(info$boundary)[1]) 
        PlotHiatusPrior(info$hiatus.max, info$hiatus.depths)
    calib.plot(info, BCAD = BCAD)
    legend("topleft", core, bty = "n", cex = 1.5)
    }
  }
  cook <- function() {
    bacon.its(ssize, burnin, info)
    txt <- paste0(info$prefix, ".bacon")
    bacon(txt, as.character(outfile), ssize + burnin, cc.dir)
    scissors(burnin, info)
    if (suppress.plots == FALSE){
    agedepth(info, BCAD = BCAD, depths.file = depths.file, 
             depths = depths, verbose = TRUE, age.unit = age.unit, 
             depth.unit = depth.unit, ...)
    }
    if (plot.pdf) 
      if (dev.interactive()) 
        dev.copy2pdf(file = paste0(info$prefix, ".pdf"))
    else {
      pdf(file = paste0(info$prefix, ".pdf"))
      agedepth(info, BCAD = BCAD, depths.file = depths.file, 
               depths = depths, verbose = FALSE, age.unit = age.unit, 
               depth.unit = depth.unit, ...)
      dev.off()
    }
  }
  write.Bacon.file(info, younger.than = younger.than, older.than = older.than)
  if (!run) 
    prepare()
  else if (!ask) 
    cook()
  else {
    prepare()
    if (accept.suggestions) 
      ans <- "y"
    else ans <- readline(message(" Run ", core, " with ", 
                                 info$K, " sections? (Y/n) "))
    ans <- tolower(substr(ans, 1, 1))[1]
    if (ans == "y" || ans == "") 
      cook()
    else message("  OK. Please adapt settings")
  }
  if (save.elbowages) {
    saved <- sapply(info$elbows, Bacon.Age.d)
    write.table(saved, paste0(info$prefix, "_", info$d.min, 
                              "_", thick, "_elbowages.txt"), row.names = FALSE, 
                col.names = FALSE)
  }
  # detach internal rbacon functions
  detach("rbacon_all")
  detach("rintcal_all")
  }


