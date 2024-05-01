library(seqtrie)
library(stringdist)
library(Rcpp)
library(dplyr)
library(qs)
library(ggplot2)

NITER <- 3

tic <- function() { .time <<- Sys.time() }
toc <- function() { as.numeric(Sys.time() - .time, units = "secs") }

encode_source <- function(file, width = 160) {
  n <- file.info(file)$size
  x <- readChar(con = file, nchars=n, useBytes = TRUE)
  x <- qserialize(x, preset = "custom", algorithm = "zstd", compress_level = 22)
  x <- base91_encode(x)
  starts <- seq(1,nchar(x), by=width)
  x <- sapply(starts, function(i) {
    substr(x, i, i+width-1)
  })
  x <- gsub('\\"', "\'", x)
  dput(x)
}

decode_source <- function(x) {
  x <- paste0(x, collapse = "")
  x <- gsub("\\'", '\\"', x)
  x <- base91_decode(x)
  qdeserialize(x)
}

# Compare with original implementation -- requires C++17
ogsource <- c(
"un]'BAAA@QRtHACAAAAAAAjn6BAABdk1kux.im:AwoDf|vw(wxy(B'B'1u@cml1Qlc1,lq.{^M^W015+H<yKma_%c3=C_*MC/<fXA'k%WLq#$$FAY'/ic,@n,wG&e5K?kmD.ax^T_1bN;Dvt&$kw=+xwINU+[>+[", 
"n&'q/(a,Z}8|b6`Vc_HnG8}n2@8gO;qKg[+LC#sWWFj~ojJ/`$5:)9v>Q^yOUnnpT(HZ4SG=weaYl2:4wzxel*j%G?*zR`[ZNZuv7RoHDx^d`tdM{l7mg82*xG0#u<xf~heiFHYMtZ<E,(6ja/#l<$Qg|iOaZ847", 
"8aDsFeTP|B|Lx(X@8QJP:.[^p5l<8&(rvveS8M{Jt49D*whd6RT!8HxyOwh15f{VErZaFyQ5Vf>n_PP*W*SW`G,Sd(4_Yr(nUF8w8^Qc{4H3u!+ia|df8=7xWIYku4xI3oAa*83Tu7/1L%9@`xQ#hir#^M,Hs{<_", 
"Me!L4*dPtzP/mD^+Eo([?Gd0'5T6o}xP$=nf(rz_H7PZjOBW|w0*!Y$6D:[L_Vx9%bsCZz<NwVxS]gIl+Ob{1:mhyJORTNip0&!C)Dl}~<?pkeLgI%]XK+sTh%n#,l>:.Ec>W35]yU/#lg1+nnFP8H7(O=.JBz!@", 
";eSC&}%r!r]Mcg4KM:o+E,2p*otfWnbKiy.[mLUp%e?e2N9mWbi1(@&lNU1qIJu=3pSQ1b+^+ir|&A||8+;<$KL`LOvG=(E9hlR^78AWq+VkJrzI(.dQ=jJcAb;`)3'.,REsKCwU|QMRRM8gwgDiS#ALR/W/Qnb5", 
"adW[k}?61TM?N4Bm8rb'|v/u<l32(R?*55'nTCn:/88a@r@X%1Z'u4aFvtbYkQ$D>J[c$+N]{6GSbS0U8MO_?byTm;/{weP8Fa.Z(d=*8ax+3@O?dJp/~xpIU(f&$ECL$&lLSwaq1t.EmrHM:aTX[DICu*^D:Ts+", 
"GEsLmu(:iG1Ez/=Nr<Z|QShD`Jai4g%n~pmB5)Z@9Qhhbc@^gIfHSZX|Du*Y$PZKEdzA/FAAKO95AApZTXn2RdpD4xXEI7<c=5+R^XCtCABtgAGANcAAY4)IT@,;;'>#OOv*1KBB[ng(U_:Y3YDVh]M{B9jpmc/Z", 
"<qPOb)%.n[Zxj/Q8(dx9DS|}hReg8zA&4:Czu[n$YyGxEGU91iG:g[tcgs.EFn!H?t{<W&u0g9Qz{HUSe2Wak60YJEhj46?oehN0A,xIM*K1*m.DO`H?,u2q+g3}Df?ga`_ajCiDal5ICNF14GJHsv{a.y2ZU@{L", 
";*KV;V'q:@=v?'$e3MbI%kd]UX*vQkcqCF'T&i`I7Mhtfx9gXE5=^Ynej~Qas_X_(&3bw8AovP&U<bXo),WQ&(Ts*'7rOJg!fVxD)BKS,f){uQ[*nF_|>rXyz<zwM;/!{0WTtUzh`3gby3!5+~xjxT?~@RvZcRXp", 
"Ua2NClD/6_E?=^[DzN3O%zq9n@3&Xnt_k+2D7`.=qHgErF~NGiz~:)gQhJJz|3@6W=PZ:#mfCEp>7;xykS~IVM'(vlyKNp_<HbjkmMv:/|z]O'lmf*'GK+.L`g*yA*h}Aveavu|rjf3&(.M+Q1GHnlqx|R*l!ZV,", 
"^d0]Vn,DxMV~/~LiA`.YChJCZ;<m?|sU3(/|n;rZImy:(*}+Y#=3gOwL4('BOE^B~V@6m|$y[N]^ykqYqw,@mJp@M__QRa{XFb%kE3=ho.bseCqQs1LaqY(oWZ{qG}_jVYB{dcy7>JIdQ,83[~@T=t.xh>s^)I3S", 
"dki(|C|x|E3ENOpuMCkQQ`0%1uSYMBP9+|Du+_z)@P)~E05zR~atrERVGz7>G*tc4.1#{5ODYvw_l6vXt6!8xec+]qfhhKH6pZqkg<`y%?cYwG)si>XY`REyahjn*N:*l_?jkIKfO%T2ay4T>k=b,sJKS7ae~9ch", 
"0|>m/o=Mxc,`6AEz*|!|6N|DR}x87Y~9ag[3'>XUt;4$Df?!O|SWPR.`mut_&Ktb*F*5g%G!Y~>10E&gS1S}uh>=66aN#@z4#=Pl|1>VOw^+H)J!kC{N_B.^I8uZl)vS)}3ESdrM^dmPi:!t1KF>GKtIWs]'k:{1", 
"oYEW$gwzt1#)IW2DqE2L+T{K@acd!%[0G.G!VF.pouqN0&J@o%4G0EJa}(sW5V3%>Tq@GuL9R!)uOvS2HWvaf'K`_G):7}(wh=]k%k%d,V#}e)pCU)kSQZ.>u}/we4f@&FB'C;z[`S?%'C%*|HTB/9L4`5,{ccQn", 
"sr9eH'cjnDxC>z@hrT]BblzvVdXn?l!T)i%__BKiBXSD%N;Uh[5D,vjj|I</ItSiH`jNI$6e{.8b1<)EVYvUT6YZ<zoCus{2~O|/X!=;rUyf32jjgn[[A{;,Q;b/C`8dYTcq+Ue^/ze0<Oj6a?IVwccBz6Qax#/a", 
"u>tCP8c!.qEBcC(ES.yRWNH>W`*l^z4a:flaa6;GzM@|R!soPT[B<0f?Wq:$(xDwCPz!>:Qic5{ca665am:(`[L<1v,qo+E~'(')Cvs;=ze_WkLDbL&PFns;eb]q375:QR8y.NAVE+Nb]OKnuQ1X9HDq.c#~U8_`", 
"L&M:NFQTXR*)1b~X,rmfxj?~C=uewOk@fVNUS=@W$L>`fy:.q5QP|5s'B}upI$,jmh[h};CXs.XWzTn*9q>apPu*@]'vq01E+,'6]yCaK#AT^J58Kn7[`O'}OP;K@.wTB7#(90H2P:OH>}*o}b'B<ur7:FF0q43X", 
"6%pS^e+fg?frqYI+T8kj8^7[@)DV}.',Gh<;oioh<.5s&F%uIP?f=EN`:4&GIB]c6`%<4Dei%S+MI@k,vK|xGhE27::kKM?HU+'gK7Jy)_%~xV<`Q!,Y$xBZxnO$E<#Lsvz#tcDSjit:9yhBC&6idCOfjQ~fREp4", 
".T4?o9^$&h9l$A)~jxwTrFYE*YV?]F0?9&ZSeSF=uJd1wswalB~2=_E[^HlDkZAcIjwx17oci_Ud1O>4~)yavI*`Vv'HD")
sourceCpp(code = decode_source(ogsource))

run_og <- function(query, target, max_distance, show_progress = F) {
  results <- og_levenshteinSearch(query, target, max_distance = max_distance)
  results$query <- query[results$query+1]
  results$target <- target[results$target+1]
  results %>% arrange(query, target)
}

# run_dnatree <- function(query, target, max_distance=NULL, max_fraction=NULL, mode = "levenshtein", show_progress = FALSE, nthreads = 8) {
#   x <- treedist::DNATree$new()
#   x$insert(target)
#   x$search(query, max_distance = max_distance, max_fraction = max_fraction, mode = mode, show_progress=show_progress, nthreads=nthreads) %>%
#     arrange(query, target)
# }

run_radixtree <- function(query, target, max_distance=NULL, max_fraction=NULL, mode = "levenshtein", show_progress = FALSE, nthreads = 8) {
  x <- seqtrie::RadixTree$new()
  x$insert(target)
  x$search(query, max_distance = max_distance, max_fraction = max_fraction, mode = mode, show_progress=show_progress, nthreads=nthreads) %>% 
    arrange(query, target)
}

run_radixforest <- function(query, target, max_distance=NULL, max_fraction=NULL, mode = "levenshtein", show_progress = FALSE, nthreads = 8) {
  x <- seqtrie::RadixForest$new()
  x$insert(target)
  x$search(query, max_distance = max_distance, max_fraction = max_fraction, mode = mode, show_progress=show_progress, nthreads=nthreads) %>% 
    arrange(query, target)
}

# run_prefixtree <- function(query, target, max_distance=NULL, max_fraction=NULL, mode = "levenshtein", show_progress = FALSE, nthreads = 8) {
#   x <- treedist::PrefixTree$new()
#   x$insert(target)
#   x$search(query, max_distance = max_distance, max_fraction = max_fraction, mode = mode, show_progress=show_progress, nthreads=nthreads) %>%
#     arrange(query, target)
# }

run_stringdist <- function(query, target, max_distance=NULL, max_fraction=NULL, nthreads = 8, show_progress = F) {
  results <- stringdist::stringdistmatrix(query, target, method = "lv", nthread=nthreads)
  results <- data.frame(query = rep(query, times=length(target)), 
                        target = rep(target, each=length(query)), 
                        distance = as.vector(results))
  results <- dplyr::filter(results, is.finite(distance))
  if(!is.null(max_distance)) {
    results <- filter(results, distance <= max_distance)
  } else {
    results <- filter(results, distance <= max_fraction * nchar(query))
  }
  results$distance <- as.integer(results$distance)
  dplyr::arrange(results, query, target)
}

# methods <- list(run_dnatree, run_radixtree, run_radixforest, run_prefixtree, run_stringdist, run_og)
# names(methods) <- c("DNATree", "RadixTree", "RadixForest", "PrefixTree", "stringdist", "OG")

methods <- list(run_radixtree, run_radixforest, run_og)
names(methods) <- c("RadixTree", "RadixForest", "OG")

# data("covid_cdr3")
cc3_subset <- sample(covid_cdr3, size = 1000)

# sd_results <- run_stringdist(cc3_subset, cc3_subset, 2)
# og_results <- run_og(cc3_subset, cc3_subset, 2)
# dt_results <- run_dnatree(cc3_subset, cc3_subset, max_distance = 2)
# rt_results <- run_radixtree(cc3_subset, cc3_subset, max_distance = 2)
# rf_results <- run_radixforest(cc3_subset, cc3_subset, max_distance = 2)
# pt_results <- run_prefixtree(cc3_subset, cc3_subset, max_distance = 2)

# stopifnot(identical(sd_results, og_results))
# stopifnot(identical(sd_results, dt_results))
# stopifnot(identical(sd_results, rt_results))
# stopifnot(identical(sd_results, rf_results))
# stopifnot(identical(sd_results, pt_results))

################################################################################

grid <- expand.grid(nseqs = c(100,300,1000,3000,10000), maxdist = c(2,3), iter = 1:NITER, method = names(methods)) %>% sample_n(nrow(.))
# grid <- filter(grid, nseqs <= 1000 | method %in% c("DNATree", "RadixTree", "RadixForest", "PrefixTree"))
grid$time <- rep(0, nrow(grid))
for(i in 1:nrow(grid)) {
  print(grid[i,])
  set.seed(grid$iter[i])
  x <- sample(covid_cdr3, size = grid$nseqs[i])
  tic()
  methods[[grid$method[i]]](x, x, max_distance = grid$maxdist[i], show_progres = TRUE)
  grid$time[i] <- toc()
  if(grid$method[i] != "OG") grid$time[i] <- grid$time[i] * 8
}
maxdist_results <- grid

grid <- expand.grid(nseqs = c(100,300,1000,3000,10000,30000), maxfrac = c(0.035,0.15), iter = 1:NITER, method = c("RadixTree", "RadixForest")) %>% sample_n(nrow(.))
grid$time <- rep(0, nrow(grid))
for(i in 1:nrow(grid)) {
  print(grid[i,])
  set.seed(grid$iter[i])
  x <- sample(covid_cdr3, size = grid$nseqs[i])
  tic()
  methods[[grid$method[i]]](x, x, max_fraction = grid$maxfrac[i], show_progres = TRUE)
  grid$time[i] <- toc()
  if(grid$method[i] != "OG") grid$time[i] <- grid$time[i] * 8
}
maxfrac_results <- grid

maxdist_results %>% group_by(nseqs, method, maxdist) %>% summarize(time = mean(time)) %>% as.data.frame %>% print
maxfrac_results %>% group_by(nseqs, method, maxfrac) %>% summarize(time = mean(time)) %>% as.data.frame %>% print

ggplot(maxfrac_results, aes(x = nseqs, y = time, color = method)) + geom_point() + geom_smooth(fill = NA) +
  scale_x_log10() +
  facet_wrap(~maxfrac, scales = "free") + 
  theme_bw(base_size = 16)

ggplot(maxdist_results, aes(x = nseqs, y = time, color = method)) + geom_point() + geom_smooth(fill = NA) +
  scale_x_log10() +
  facet_wrap(~maxdist, scales = "free") + 
  theme_bw(base_size = 16)
# ggsave(g, file = "benchmark_plot.png", width = 6, height = 4)
