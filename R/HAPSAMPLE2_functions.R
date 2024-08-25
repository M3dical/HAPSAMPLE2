vcfR_parallel = function(indir, infilename, cols, n_variants, minibatch_size, type = "gt") {
  n_minibatches = trunc(n_variants/minibatch_size) + 1*(n_variants %% minibatch_size > 0)
  u = pmin((0:(n_minibatches-1))*minibatch_size+1, n_variants)
  v = pmin((1:n_minibatches)*minibatch_size,  n_variants)
  vcf_in = foreach::foreach(b = 1:n_minibatches, .packages = c("vcfR"), .combine = "rbind") %dopar% {
    try({
      vcf_mb = vcfR::read.vcfR(paste0(indir, "/", infilename), skip = u[b] - 1, nrows = v[b] - (u[b] - 1), cols = cols)
      if(type == "fix") {
        vcf_mb@fix
      }
      if(type == "gt") {
        vcf_mb@gt
      }
    })
  }
  return(vcf_in)
}

gt_to_hap = function(x) {
  n_variants = length(x)
  out1 = rep(0L, n_variants)
  out2 = rep(0L, n_variants)
  out1[x == "1|0" | x == "1|1" | x == "1/0" | x == "1/1"] = 1L
  out2[x == "0|1" | x == "1|1" | x == "0/1" | x == "1/1"] = 1L
  out = cbind(out1, out2)
}

vcf_to_rds = function(vcfObj_gt, outdir, chrom_name, sample_info) {
  n_ind = ncol(vcfObj_gt)
  n_biallelic = nrow(vcfObj_gt)
  ind_name = colnames(vcfObj_gt)
  pop_vec = sample_info[match(ind_name, sample_info[, "sample"]), "pop"]
  pop_u = sort(unique(sample_info[, "pop"]))
  n_pop = length(pop_u)
  count_by_pop = rep(0, n_pop)
  names(count_by_pop) = pop_u
  ma_counts = matrix(0, n_biallelic, n_pop)
  colnames(ma_counts) = pop_u
  ma_count_by_hap = rep(NA, n_ind*2)
  if(sum(is.na(pop_vec)) > 0){
    stop("Individual(s) not found in Sample Information File")
  }
  for (i in 1:n_ind) {
    hap_mat = gt_to_hap(vcfObj_gt[,i])
    ind_name_i = ind_name[i]
    ind_pop_i = pop_vec[i]
    hap_vec_1 = which(hap_mat[,1] == 1)
    hap_vec_2 = which(hap_mat[,2] == 1)
    outfilename1 = paste0(paste(chrom_name, ind_name_i, 1, sep = "_"), ".rds")
    outfilename2 = paste0(paste(chrom_name, ind_name_i, 2, sep = "_"), ".rds")
    outfile1 = paste0(outdir, "/", outfilename1)
    outfile2 = paste0(outdir, "/", outfilename2)
    saveRDS(hap_vec_1, file = outfile1, compress = F)
    saveRDS(hap_vec_2, file = outfile2, compress = F)
    ma_counts_i = rowSums(hap_mat)
    ma_counts[, ind_pop_i] = ma_counts[, ind_pop_i] + ma_counts_i
    count_by_pop[ind_pop_i] = count_by_pop[ind_pop_i] + 2
    ma_count_by_hap[i*2-1] = length(hap_vec_1)
    ma_count_by_hap[i*2] = length(hap_vec_2)
  }
  attr(ma_counts, "count_by_pop") = count_by_pop
  attr(ma_counts, "ma_count_by_hap") =  ma_count_by_hap
  return(ma_counts)
}

vcfFile_to_rds = function(indir, outdir, infilename, max_ram, max_inds, sample_filename, print_ind_count, n_cores = 4, n_fixed_cols) {
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl, cores = n_cores)
  sample_info = utils::read.table(paste0(indir,sample_filename), sep = "\t", header = T, stringsAsFactors = F)
  fix_cols = c(1:(n_fixed_cols - 1))
  vcf_fix = vcfR::read.vcfR(paste0(indir,infilename), cols = fix_cols)
  bi_allelic = nchar(vcf_fix@fix[, "REF"]) == 1 & nchar(vcf_fix@fix[, "ALT"]) == 1
  n_biallelic = sum(bi_allelic)
  ind_names = sample_info[, "sample"]
  pop_u = sort(unique(sample_info[, "pop"]))
  n_pop = length(pop_u)
  count_by_pop = rep(0, n_pop)
  names(count_by_pop) = pop_u
  ma_counts = matrix(0, n_biallelic, n_pop)
  ma_freqs = matrix(NA, n_biallelic, n_pop)
  colnames(ma_counts) = pop_u
  colnames(ma_freqs) = pop_u
  ma_count_by_hap = rep(NA, max_inds*2)
  n_variants = nrow(vcf_fix@fix)
  batch_size = trunc(max_ram/.1 * 10^6*10^1/n_variants)
  n_batch = trunc(max_inds/batch_size) + 1*(max_inds %% batch_size > 0)
  chrom_num = as.numeric(vcf_fix@fix[1, "CHROM"])
  chrom_name = paste0("chr", chrom_num)
  loci_names = vcf_fix@fix[bi_allelic, "ID"]
  loci_pos = vcf_fix@fix[bi_allelic, "POS"]
  meta_data = c(batch_size, n_batch, max_ram, max_inds, sample_filename)
  names(meta_data) = c("batch_size", "n_batch", "max_ram", "max_inds", "sample_filename")
  meta_data_file = paste0(chrom_name, "_meta.rds")
  saveRDS(meta_data, file = paste0(outdir, meta_data_file))
  loci_info = data.frame(loci_names, as.integer(loci_pos))
  colnames(loci_info) = c("name", "pos")
  loci_info_file = paste0(chrom_name, "_loci_info.rds")
  saveRDS(loci_info, file = paste0(outdir, loci_info_file))
  ma_freqs_file = paste0(chrom_name, "_freqs.rds")
  ma_count_by_hap_file = paste0(chrom_name, "_ma_count.rds")
  s = 0
  t = 0
  for (b in 1:n_batch) {
    s = t + 1
    t = min(t + batch_size*2, max_inds*2)
    cols = c(n_fixed_cols + (1:batch_size)+(b-1)*batch_size)
    cols = cols[cols <= max_inds + n_fixed_cols]
    minibatch_size = trunc(n_variants/n_cores)+1
    vcf_b = vcfR_parallel(indir = indir, infilename = infilename, cols = cols, n_variants = n_variants, minibatch_size = minibatch_size, type = "gt")
    vcfObj_gt = vcf_b[bi_allelic,]
    ma_counts_b = vcf_to_rds(vcfObj_gt = vcfObj_gt, outdir = outdir, chrom_name = chrom_name, sample_info = sample_info)
    ma_counts_file_b = paste0(chrom_name, "_", "ma_counts_B", b, ".rds")
    saveRDS(ma_counts_b, file = paste0(outdir, ma_counts_file_b))
    ma_counts[] = ma_counts + ma_counts_b
    count_by_pop_b = attr(ma_counts_b, "count_by_pop")
    ma_count_by_hap_b = attr(ma_counts_b, "ma_count_by_hap")
    count_by_pop = count_by_pop + count_by_pop_b
    ma_count_by_hap[s:t] = ma_count_by_hap_b
    ma_count_by_hap_file_b = paste0(chrom_name, "_", "ma_count_by_hap_B", b, ".rds")
    saveRDS(ma_count_by_hap_b, file = paste0(outdir, ma_count_by_hap_file_b))
    if(print_ind_count == TRUE) {
      print(paste0(min(batch_size*b, max_inds), " individual .rds files written for Chr", chrom_num))
    }
  }
  for (k in 1:n_pop) {
    ma_freqs[,k] = ma_counts[,k]/count_by_pop[k]
  }
  saveRDS(ma_freqs, file = paste0(outdir, ma_freqs_file))
  saveRDS(ma_count_by_hap, file = paste0(outdir, ma_count_by_hap_file))
  parallel::stopCluster(cl)
}

retAlleleFreqs = function(dir, infiles, n_loci) {
  n_files = length(infiles)
  counts = foreach::foreach(i = 1:n_files, .combine = cbind) %dopar% {
    zeroes = matrix(rep(0L, n_loci), nrow = n_loci, ncol = 1)
    inhap = readRDS(paste0(dir,infiles[i]))
    zeroes[inhap, ] = 1L
    return(zeroes)
  }
  freqs = rowMeans(counts)
  return(freqs)
}

writeAlleleFreqs_W = function(dir, chr_vec, sample_filename, n_cores = 4) {
  n_chr = length(chr_vec)
  sample_info = utils::read.table(paste0(dir, sample_filename), sep = "\t", header = T, stringsAsFactors = FALSE)
  n_ind = nrow(sample_info)
  hap_tags = paste0(rep(sample_info[, "sample"], 2), "_", sort(rep(c(1, 2), each = n_ind)))
  hap_tags = sort(hap_tags)
  pops = sort(unique(sample_info[, "pop"]))
  pop_counts = table(sample_info[, "pop"])
  n_pops = length(pops)
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl, cores=n_cores)
  for (h in 1:n_chr) {
    chr_h = chr_vec[h]
    chrom_tag = paste0("chr", chr_h)
    loci_info_h = readRDS(paste0(dir, chrom_tag, "_loci_info.rds"))
    n_loci_h = nrow(loci_info_h)
    freqs_h = matrix(0, nrow = n_loci_h, ncol = n_pops)
    for (k in 1:n_pops) {
      pop_k = pops[k]
      sample_name = substr(hap_tags, 1, (nchar(hap_tags)-2))
      inds_hk = hap_tags[sample_name %in% sample_info[sample_info[, "pop"] %in% pop_k, "sample"]]
      infiles_hk = paste0(chrom_tag, "_", inds_hk, ".rds")
      freqs_hk = retAlleleFreqs(dir = dir, infiles = infiles_hk, n_loci = n_loci_h)
      freqs_h[, k] = freqs_hk
      outfile_hk = paste0(chrom_tag, "_", pop_k, "_freqs.rds")
      saveRDS(freqs_hk, paste0(dir, outfile_hk))
    }
    outfile_h = paste0(chrom_tag, "_", "ALL", "_freqs.rds")
    all_freqs_h = rowMeans(freqs_h)
    saveRDS(all_freqs_h, paste0(dir, outfile_h))
  }
  parallel::stopCluster(cl)
}

sampleChromo = function(chr, indir, outdir, sample_info, retain_every = 5, retain_snp = NULL, maf_threshold = 0.02) {
  chrom_tag = paste0("chr", chr)
  n_ind = nrow(sample_info)
  n_hap = n_ind*2
  n_snp_chr = 0
  loci_info_file = paste0(chrom_tag, "_loci_info.rds")
  loci_info = as.data.frame(readRDS(paste0(indir,loci_info_file)), stringsAsFactors = FALSE)
  loci_names = loci_info[, "name"]
  n_loci_orig = nrow(loci_info)
  orig_range = 1:n_loci_orig
  retain_snp_chr = NULL
  pops = unique(sample_info[, "pop"])
  n_pops = length(pops)
  pop_freq_files = paste0(chrom_tag, "_", pops, "_freqs.rds")
  all_freqs_file = paste0(chrom_tag, "_ALL", "_freqs.rds")
  freqs = readRDS(paste0(indir,all_freqs_file))
  meets_thresh = as.numeric(freqs >= maf_threshold & (1-freqs) >= maf_threshold)
  order_map = data.frame(orig = orig_range, meets_thresh = meets_thresh)
  if (is.null(retain_snp) == FALSE) {
    retain_snp = unique(retain_snp)
    retain_snp_chr = retain_snp[retain_snp %in% loci_names]
    n_snp_chr = length(retain_snp_chr)
    if (n_snp_chr > 0) {
      retain_snp_order = match(retain_snp_chr, loci_names)
      order_map_retain = order_map[retain_snp_order,]
    }
  }
  order_map = order_map[order_map[, "meets_thresh"] == 1,]
  pt_range = 1:nrow(order_map)
  order_map = order_map[pt_range %% retain_every == 0,]
  if (n_snp_chr > 0) {
    order_map = rbind(order_map, order_map_retain[!(order_map_retain[, "orig"] %in% order_map[, "orig"]),])
  }
  order_map = order_map[order(order_map[, "orig"]),]
  order_map[, "post_re"] = 1:nrow(order_map)
  order_map[, "new"] = order_map[, "post_re"]
  hap_tags = paste0(rep(sample_info[, "sample"], 2), "_", sort(rep(c(1, 2), each = n_ind)))
  infiles = paste0(indir, chrom_tag, "_", hap_tags, ".rds")
  outfiles = paste0(outdir, chrom_tag, "_", hap_tags, ".rds")
  loci_info_out = loci_info[order_map[, "orig"],]
  saveRDS(loci_info_out, paste0(outdir,loci_info_file))
  retain_loci = order_map[, "orig"]
  for (k in 1:n_pops) {
    in_freq_k = readRDS(paste0(indir, pop_freq_files[k]))
    out_freq_k = in_freq_k[retain_loci]
    saveRDS(out_freq_k, paste0(outdir, pop_freq_files[k]))
  }
  out_freqs = freqs[retain_loci]
  saveRDS(out_freqs, paste0(outdir, all_freqs_file))
  catchOut = foreach::foreach(i = 1:n_hap) %dopar% {
    try({
      inhap_i = readRDS(paste0(infiles[i]))
      retain_i = inhap_i[inhap_i %in% retain_loci]
      outhap_i = order_map[match(retain_i, order_map[, "orig"]), "new"]
      saveRDS(outhap_i, outfiles[i])
    })
  }
}

sampleChromo_W = function(chr_vec, indir, outdir, sample_filename, retain_every = 5, retain_snp = NULL, maf_threshold = 0.02, n_cores = 4) {
  if (maf_threshold >= .5) { stop("maf_threshold should be < 0.5") }
  n_chr = length(chr_vec)
  sample_info = utils::read.table(paste0(indir,sample_filename), sep = "\t", header = T, stringsAsFactors = FALSE)
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  for (i in 1:n_chr) {
    chr_i = chr_vec[i]
    catchOut = sampleChromo(chr = chr_i, indir = indir, outdir = outdir, sample_info = sample_info, retain_every = retain_every, retain_snp = retain_snp, maf_threshold = maf_threshold)
  }
  parallel::stopCluster(cl)
}

makeLociInfo = function(chr, dir) {
  chrom_tag = paste0("chr", chr)
  chr_gen_map <- utils::read.table(paste0(dir,"genetic_map_",chrom_tag,"_combined_b37.txt"),sep="",header=T,stringsAsFactors=F)
  chr_gen_map <- chr_gen_map[,c(1,3)]
  colnames(chr_gen_map) <- c("pos","gen_map_cm")
  loci_info_file = paste0(chrom_tag, "_loci_info.rds")
  loci_info = as.data.frame(readRDS(paste0(dir,loci_info_file)), stringsAsFactors = F)
  n_missing = sum(is.na(loci_info[, "name"]))
  names_for_missing = paste0(chrom_tag, "_missing_", 1:n_missing)
  loci_info[is.na(loci_info[, "name"]), "name"] = names_for_missing
  loci_names = loci_info[, "name"]
  n_biallelic = nrow(loci_info)
  loci_info[, "gen_map_cm"] = NA
  loci_info[, "gen_map_cm"] = chr_gen_map[match(loci_info[, "pos"], as.character(chr_gen_map[, "pos"])), "gen_map_cm"]
  gen_map_lm = stats::lm(gen_map_cm ~ pos, chr_gen_map)
  gen_map_int = summary(gen_map_lm)$coefficients[1, "Estimate"]
  gen_map_coef = summary(gen_map_lm)$coefficients[2, "Estimate"]
  gen_map_min = min(loci_info[, "gen_map_cm"], na.rm = T)
  gen_map_max = max(loci_info[, "gen_map_cm"], na.rm = T)
  min_pos = as.numeric(loci_info[which(loci_info[, "gen_map_cm"] == gen_map_min)[1], "pos"])
  max_pos = as.numeric(loci_info[which(loci_info[, "gen_map_cm"] == gen_map_max)[1], "pos"])
  if (is.na(loci_info[1, "gen_map_cm"])) {
    loci_info[1, "gen_map_cm"] = gen_map_min + gen_map_coef*(as.numeric(loci_info[1, "pos"])-min_pos)
  }
  if (is.na(loci_info[n_biallelic, "gen_map_cm"])) {
    loci_info[n_biallelic, "gen_map_cm"] = gen_map_max + gen_map_coef*(as.numeric(loci_info[n_biallelic, "pos"])-max_pos)
  }
  loci_info[, "gen_map_interp"] = stats::approx(x = as.numeric(loci_info[, "pos"]), y = loci_info[, "gen_map_cm"], xout = as.numeric(loci_info[, "pos"]))$y
  chr_cross = diff(loci_info[, "gen_map_interp"], lag = 1)
  chr_recomb = chr_cross/100
  out = list(loci_info = loci_info, chr_recomb = chr_recomb)
  return(out)
}

makeLociInfo_W = function(chr_vec, dir) {
  n_chr = length(chr_vec)
  loci_info = vector(mode = "list", length = n_chr)
  chr_recomb = vector(mode = "list", length = n_chr)
  for (i in 1:n_chr) {
    out_i = makeLociInfo(chr = chr_vec[i], dir = dir)
    loci_info[[i]] = out_i[["loci_info"]][, c(1:2)]
    chr_recomb[[i]] = out_i[["chr_recomb"]]
  }
  names(loci_info) = chr_vec
  names(chr_vec) = chr_recomb
  out = list(loci_info = loci_info, chr_recomb = chr_recomb)
}

calcChrAlleleFreqs = function(dir, snp_order, chr, sample_info, n_cores = 4) {
  sample_names = sample_info$sample
  chrom_tag = paste0("chr", chr)
  chr_sample_names = paste0(chrom_tag, "_", sample_names)
  chr_sample_hap = expand.grid(chr_sample_names, c(1,2), stringsAsFactors = FALSE)
  chr_sample_files = paste0(chr_sample_hap[,1], "_", chr_sample_hap[,2], ".rds", sep = "")
  all_files = list.files(path = dir)
  chr_sample_files_avail = intersect(all_files, chr_sample_files)
  n_hap = length(chr_sample_files_avail)
  ind_avail = unlist(strsplit(chr_sample_files_avail, split = "_", fixed = T))
  ind_avail = ind_avail[-grep(".rds", ind_avail, fixed = T)]
  ind_avail = ind_avail[-grep("chr", ind_avail, fixed = T)]
  ind_avail_u = unique(ind_avail)
  n_ind = length(unique(ind_avail_u))
  sample_ma_info = as.data.frame(matrix(NA, nrow = n_hap, ncol = length(snp_order) + 4))
  colnames(sample_ma_info) = c("pop", "ind", "file", "ma_count", names(snp_order))
  sample_ma_info[, "ind"] = ind_avail
  sample_ma_info[, "file"] = chr_sample_files_avail
  sample_ma_info[, "pop"] = sample_info[match(sample_ma_info[, "ind"], sample_info[, "sample"]), "pop"]
  n_pop = length(unique(sample_ma_info[, "pop"]))
  for (i in 1:n_hap){
    hap_i = readRDS(paste0(dir,sample_ma_info[i, "file"]))
    sample_ma_info[i, c("ma_count", names(snp_order))] =  c(length(hap_i), snp_order %in% hap_i)
  }
  print_msg = paste0("Chr", chr, " min. allele freqs for ",
                     n_ind, " individual samples, from ",
                     n_pop, " populations, calculated")
  print(print_msg)
  out = sample_ma_info
  return(out)
}

calcChrAlleleFreqs_W = function(snp, loci_info, chr_vec, dir, sample_filename, n_cores=4, elig_pops) {
  sample_info = utils::read.table(paste0(dir, sample_filename), sep = "\t", header = T, stringsAsFactors = FALSE)
  all_pops = sort(unique(sample_info[, "pop"]))
  if (sum(elig_pops %in% all_pops) < length(elig_pops)){
    unfound_pops = elig_pops[!(elig_pops %in% all_pops)]
    out_msg = paste0("elig_pops contains populations not found in the sample info file")
    stop(out_msg)
  }
  sample_info = sample_info[sample_info[, "pop"] %in% elig_pops,]
  n_snp = length(snp)
  n_chr = length(chr_vec)
  snp_map = as.data.frame(matrix(NA, nrow = n_snp, ncol = 3))
  colnames(snp_map) = c("chr", "snp", "order")
  snp_map[, "snp"] = snp
  sample_ma_list = vector(mode = "list", n_chr)
  for (i in 1:n_chr) {
    loci_info_i = loci_info[[i]]
    snp_map[snp_map[, "snp"] %in% loci_info_i[, "name"], "chr"] = chr_vec[i]
  }
  for (i in 1:n_chr) {
    loci_info_i = loci_info[[i]]
    snp_map[snp_map[, "chr"] == chr_vec[i], "order"] = match(snp_map[snp_map[, "chr"] == chr_vec[i], "snp"], loci_info_i[, "name"])
  }
  for (i in 1:n_chr) {
    chr_i = chr_vec[i]
    snp_order_i = snp_map[snp_map[, "chr"] == chr_i, "order"]
    names(snp_order_i) = snp_map[snp_map[, "chr"] == chr_i, "snp"]
    sample_ma_i = calcChrAlleleFreqs(dir, snp_order = snp_order_i, chr = chr_vec[i], sample_info = sample_info, n_cores = n_cores)
    sample_ma_list[[i]] = sample_ma_i
  }
  snp_map = snp_map[order(snp_map[, "chr"], snp_map[, "order"]),]
  names(sample_ma_list) = chr_vec
  out = list(snp_map = snp_map, sample_ma_list = sample_ma_list)
  return(out)
}

calcAlleleFreq_by_pop = function(sample_ma_list, snp_map) {
  snp = snp_map[, "snp"]
  chr_vec = sort(unique(snp_map[, "chr"]))
  n_chr = length(chr_vec)
  n_snp = nrow(snp_map)
  for (i in 1:n_chr){
    chr_i = as.character(chr_vec[i])
    if(i == 1){
      pop_u = unique(sample_ma_list[[chr_i]][, "pop"])
      n_pop = length(pop_u)
      ma_freq_by_pop = matrix(NA, nrow = n_pop, ncol = n_snp)
      rownames(ma_freq_by_pop) = sort(pop_u)
      colnames(ma_freq_by_pop) = snp
    }
    snp_i = snp[snp_map[, "chr"] == chr_vec[i]]
    sample_ma_i = cbind(sample_ma_list[[chr_i]][, "pop", drop = F], sample_ma_list[[chr_i]][, -c(1:4), drop = F])
    ma_i = stats::aggregate(. ~ pop, sample_ma_i, FUN = mean)
    ma_freq_by_pop[, snp_i] = as.matrix(ma_i[, snp_i])
  }
  out = t(ma_freq_by_pop)
}

makeGenFreqAnc = function(ma_freq) {
  n_loci = nrow(ma_freq)
  n_pop = ncol(ma_freq)
  gen_freq = vector(mode = "list", n_loci)
  for (i in 1:n_loci) {
    gen_freq[[i]] = matrix(0, nrow = 3, ncol = n_pop)
    gen_freq[[i]][1,] = (1 - ma_freq[i,])^2
    gen_freq[[i]][3,] = (ma_freq[i,])^2
    gen_freq[[i]][2,] = 1 - gen_freq[[i]][3,] - gen_freq[[i]][1,]
  }
  gf_by_anc = matrix(0, nrow = n_loci*3, n_pop)
  for (k in 1:n_pop) {
    for (i in 1:length(gen_freq)) {
      gf_by_anc[(i*3-2):(i*3),k] = gen_freq[[i]][,k]
    }
  }
  return(gf_by_anc)
}

mostFrequentGenos = function(pg_mat, max_combs) {
  n_snp = nrow(pg_mat)
  combs_cutoff = min(max_combs, 3^n_snp)
  g_mat = matrix(c(0,1,2), nrow = n_snp, ncol = 3, byrow = T)
  pg_maxes = matrixStats::rowMaxs(pg_mat)
  pg_maxes_ord = order(-pg_maxes)
  pg_mat_ord = pg_mat
  g_mat_ord = g_mat
  pg_mat_ord0 = pg_mat[pg_maxes_ord,, drop = F]
  g_mat_ord0 = g_mat[pg_maxes_ord,, drop = F]
  for (i in 1:n_snp) {
    row_ordered_i = order(-pg_mat_ord0[i,])
    pg_mat_ord[i,] = pg_mat_ord0[i, row_ordered_i]
    g_mat_ord[i,] = g_mat_ord0[i, row_ordered_i]
  }
  S = min(trunc(log(combs_cutoff)/log(3))+1, n_snp)
  ord_list = vector(mode = "list", length = S)
  for (s in 1:S) {
    ord_list[[s]] = c(1:3)
  }
  ord_grid0 = as.matrix(expand.grid(ord_list))
  ord_grid0 = ord_grid0[, rev(1:S), drop = F]
  if(n_snp - S > 0) {
    ones = matrix(1, nrow = combs_cutoff, ncol = n_snp - S)
    ord_grid = cbind(ones, ord_grid0[1:combs_cutoff,, drop = F])
  } else {
    ord_grid = ord_grid0[1:combs_cutoff,, drop = F]
  }
  colnames(ord_grid) = paste0("s", 1:n_snp)
  g_combs0 = matrix(NA, combs_cutoff, n_snp)
  g_combs = matrix(NA, combs_cutoff, n_snp)
  pg_vec = rep(NA, combs_cutoff)
  for (g in 1:nrow(g_combs0)) {
    ind_g = cbind(seq_along(ord_grid[g,]), ord_grid[g,])
    g_combs0[g,] = g_mat_ord[ind_g]
    pg_vec[g] = prod(pg_mat_ord[ind_g])
  }
  for (i in 1:n_snp) {
    g_combs[, pg_maxes_ord[i]] = g_combs0[, i]
  }
  out = list(pg_vec = pg_vec, g_combs = g_combs)
}

solveB0_both = function(B0, nu, pg_vec, mu_k, type) {
  z_vec = B0 + nu
  if (type == "bin") {
    rhs =  sum(pg_vec * exp(z_vec)/(1+exp(z_vec)))
  }
  if (type == "cont") {
    rhs =  sum(pg_vec * z_vec)
  }
  delta = abs(mu_k - rhs)
  return(delta)
}

solveB0Est_both = function(mu_vec, gf_by_anc, B1, type, lower = -5, upper = 5, max_combs = 100000) {
  n_pop = length(mu_vec)
  B0Est = rep(NA, n_pop)
  pg_vec_out = vector(mode = "list", length = n_pop)
  g_combs_out = vector(mode = "list", length = n_pop)
  for(k in 1:n_pop) {
    mu_k = mu_vec[k]
    gen_k = gf_by_anc[,k]
    B0Est_init = 0
    p_g = gen_k
    if(max_combs == Inf) {
      L = length(p_g)/3
      g_vec = rep(c(0, 1, 2), times = L)
      g_list = vector(mode = "list", length = L)
      pg_list = vector(mode = "list", length = L)
      b = 0
      a = 0
      for (i in 1:L) {
        a = b + 1
        b = i * 3
        g_list[[i]] = c(0, 1, 2)
        pg_list[[i]] = p_g[a:b]
      }
      g_combs = as.matrix(expand.grid(g_list))
      pg_combs = expand.grid(pg_list)
      pg_vec = matrixStats::rowProds(as.matrix(pg_combs))
    } else {
      pg_mat = matrix(p_g, ncol = 3, byrow = T)
      mFG_out = mostFrequentGenos(pg_mat = pg_mat, max_combs = max_combs)
      pg_vec = mFG_out[["pg_vec"]]
      g_combs = mFG_out[["g_combs"]]
    }
    pg_vec_adj = pg_vec/sum(pg_vec)
    nu = g_combs%*%B1
    B0Est[k] = stats::optimize(solveB0_both, nu = nu, pg_vec = pg_vec_adj, mu_k = mu_k, type = type, lower = lower, upper = upper)$minimum
    pg_vec_out[[k]] = pg_vec
    g_combs_out[[k]] = g_combs
  }
  out = list(B0Est = B0Est, pg_vec_out = pg_vec_out, g_combs_out = g_combs_out)
  return(out)
}

solveB0Est_W = function(mu_vec, gf_by_anc, B1, type = "bin", lower = -8, upper = 8) {
  elig_pops = colnames(gf_by_anc)
  n_elig_pops = length(elig_pops)
  n_snp = length(B1)
  B0Est = matrix(NA, nrow = n_snp, ncol = n_elig_pops)
  a = 0
  b = 0
  for (s in 1:n_snp) {
    a = b + 1
    b = b + 3
    solveOut = solveB0Est_both(mu_vec = mu_vec, gf_by_anc = gf_by_anc[a:b,,drop = F], B1 = B1[s], type = type, lower = lower, upper = upper, max_combs = Inf)
    B0Est[s,] = solveOut[["B0Est"]]
    colnames(B0Est) = elig_pops
  }
  return(B0Est)
}

calc_pG_AD_s = function(B0Est_s, B1_s, spec_prev, gf_by_anc_s, adm_mat) {
  n_inds = nrow(adm_mat)
  elig_pops = colnames(gf_by_anc_s)
  ind_int = adm_mat[, elig_pops]%*%B0Est_s
  z_mat = matrix(NA, nrow = n_inds, ncol = 3)
  for (g in 1:3) {
    z_mat[,g] = ind_int + (g - 1)*B1_s
  }
  spec_prev_LO = log(spec_prev/(1-spec_prev))
  pG_A = adm_mat[, elig_pops]%*% t(gf_by_anc_s)
  pD_A_LO = matrix(as.numeric(adm_mat[, elig_pops]%*% spec_prev_LO), nrow = n_inds, ncol = 3)
  pD_A = exp(pD_A_LO)/(1+exp(pD_A_LO))
  pG_AD = exp(z_mat)/(1+exp(z_mat))*pG_A/pD_A
  return(pG_AD)
}

drawGenotypes = function(pG_AD_s) {
  n_inds = nrow(pG_AD_s)
  out = rep(NA, n_inds)
  for (i in 1:n_inds) {
    prob_i = pG_AD_s[i,]
    out[i] = sample.int(3, size = 1, replace = T, prob = prob_i)-1
  }
  return(out)
}

drawGenotypes_W = function(B0Est, B1, spec_prev, gf_by_anc, adm_mat) {
  snp = names(B1)
  n_snp = length(snp)
  elig_pops = colnames(gf_by_anc)
  n_inds = nrow(adm_mat)
  drawnGenos = matrix(NA, nrow = n_inds, ncol = n_snp)
  colnames(drawnGenos) = snp
  targAlleleCounts = matrix(NA, nrow = n_inds, ncol = n_snp)
  p_Alt_AD = matrix(NA, nrow = n_inds, ncol = n_snp)
  colnames(p_Alt_AD) = snp
  a = 0
  b = 0
  for (s in 1:n_snp) {
    a = b + 1
    b = b + 3
    B0Est_s = B0Est[s,]
    B1_s = B1[s]
    gf_by_anc_s = gf_by_anc[a:b,]
    pG_AD_s = calc_pG_AD_s(B0Est_s = B0Est_s, B1_s = B1_s, spec_prev = spec_prev, gf_by_anc_s = gf_by_anc_s, adm_mat = adm_mat)
    draw_s = drawGenotypes(pG_AD_s = pG_AD_s)
    drawnGenos[,s] = draw_s
    targAlleleCounts[,s] = rowSums(pG_AD_s)
    p_Alt_AD[,s] = pG_AD_s%*%c(0,1,2)/2
  }
  out = list(drawnGenos = drawnGenos, targAlleleCounts = targAlleleCounts, p_Alt_AD = p_Alt_AD)
  return(out)
}

createHap = function(dir, n_loci, adm_ind, start_pos, start_df, oth_df, recomb_prob, max_ma_count, read_type = "rds", recomb_sites = NULL,
                     hap_vec_mat = NULL) {
  if (is.null(recomb_sites) == TRUE) {
    recomb = stats::runif(n_loci - 1) < recomb_prob
    n_recomb = sum(recomb)
    if (n_recomb > 0) {
      recomb_sites = which(recomb == 1)
    } else {
      recomb_sites = n_loci
    }
  } else {
    n_recomb = length(recomb_sites)
  }
  recomb_sites = sort(unique(recomb_sites))
  n_piece = n_recomb + 1
  pos_mat = matrix(NA, nrow = n_piece, ncol = 2)
  colnames(pos_mat) = c("begin", "end")
  pos_df = as.data.frame(pos_mat)
  pos_df[n_piece, "end"] = n_loci
  pos_df[1, "begin"] = 1
  if (n_recomb > 0) {
    pos_df[1:(n_recomb), "end"] = recomb_sites
    pos_df[-1, "begin"] = recomb_sites + 1
  }
  pos_df[, "pop"] = sample(names(adm_ind), n_piece, replace = TRUE, prob = adm_ind)
  pos_df[, "start"] = 0
  pos_df[(pos_df[, "begin"] <= start_pos & pos_df[,"end"] >= start_pos), "start"] = 1
  pos_df[,"file"] = NA
  pos_df[,"ma_count"] = NA
  setwd(dir)
  hap_out = rep(NA, max_ma_count)
  q = 0
  r = 0
  for (m in 1:n_piece) {
    s = pos_df[m, "begin"]
    t = pos_df[m, "end"]
    pop_m = pos_df[m, "pop"]
    if (pos_df[m, "start"] == 1){
      files_m = start_df[start_df[, "pop"] == pop_m, "file"]
      weights_m = start_df[start_df[, "pop"] == pop_m, "weight"]
    } else {
      files_m = oth_df[oth_df[, "pop"] == pop_m, "file"]
      weights_m = oth_df[oth_df[, "pop"] == pop_m, "weight"]
    }
    sel_file_m = sample(files_m, size = 1, prob = weights_m)
    pos_df[m, "file"] = sel_file_m
    if (read_type == "ram") {
      hap_m = hap_vec_mat[, sel_file_m]
    } else {
      hap_m = readRDS(sel_file_m)
    }
    hap_m_piece = hap_m[hap_m >= s & hap_m <= t]
    hap_m_piece_len = length(hap_m_piece)
    pos_df[m,"ma_count"] = hap_m_piece_len
    if (hap_m_piece_len > 0) {
      q = r + 1
      r = q + hap_m_piece_len - 1
      hap_out[q:r] = hap_m_piece
    }
  }
  attr(hap_out, "pos_df") = pos_df
  return(hap_out)
}

createGeno = function(dis_chr, dir, n_loci, adm_ind, start_geno, start_pos, maj_sample_df, min_sample_df, all_sample_df, read_type = "rds",
                      hap_vec_mat = NULL, recomb_sites1 = NULL, recomb_sites2 = NULL, recomb_prob = NULL, max_ma_count) {
  start_df1 = all_sample_df
  start_df2 = all_sample_df
  oth_df = all_sample_df
  if (dis_chr == "yes") {
    if (start_geno == 0) {
      start_df1 = maj_sample_df
      start_df2 = maj_sample_df
    } else if (start_geno == 1) {
      start_df1 = min_sample_df
      start_df2 = maj_sample_df
    } else {
      start_df1 = min_sample_df
      start_df2 = min_sample_df
    }
  }
  hap1 = createHap(dir, n_loci, adm_ind, start_pos, start_df1, oth_df, recomb_prob, max_ma_count, read_type = read_type, recomb_sites = recomb_sites1,
                   hap_vec_mat = NULL)
  hap2 = createHap(dir, n_loci, adm_ind, start_pos, start_df2, oth_df, recomb_prob, max_ma_count, read_type = read_type, recomb_sites = recomb_sites2,
                   hap_vec_mat = NULL)
  hap1_long = rep(0, n_loci)
  hap2_long = rep(0, n_loci)
  hap1_long[hap1] = 1
  hap2_long[hap2] = 1
  geno_long = hap1_long + hap2_long
  return(geno_long)
}

calcPoisDiff = function(t, lambda){
  p0 = stats::dpois(0, lambda-t)
  abs_diff = abs((p0*(1) - t))
}

recomb_poisApprox = function(n_hap, recomb_prob) {
  lambda = sum(recomb_prob)
  optOut = stats::optimize(calcPoisDiff, interval = c(0, lambda), lambda = lambda)
  t = optOut$minimum
  rpoisOut = stats::rpois(n_hap, lambda - t)
  n_targ_recomb = pmax(rpoisOut, 1)
  recomb_sites = sample(1:length(recomb_prob), size = sum(n_targ_recomb), prob = recomb_prob, replace = T)
  attr(recomb_sites, "n_targ_recomb") = n_targ_recomb
  return(recomb_sites)
}

createGeno_W = function(dis_chr = "no", dir, chr_sample_df, adm_mat, start_genos, start_pos, recomb_prob, max_ma_count, chr, prefix, write_genos,
                        return_genos, loci_out = NULL, compress = F) {
  chrom_tag = paste0("chr", chr)
  n_ind = nrow(adm_mat)
  n_hap = n_ind*2
  n_loci = length(recomb_prob) + 1
  all_sample_df = chr_sample_df
  maj_sample_df = chr_sample_df[chr_sample_df[, "min_allele"] == 0,]
  min_sample_df = chr_sample_df[chr_sample_df[, "min_allele"] == 1,]
  all_recomb_sites = recomb_poisApprox(n_hap = n_hap, recomb_prob = recomb_prob)
  n_targ_recomb = attr(all_recomb_sites, "n_targ_recomb")
  recomb_index = rep(1:n_hap, times = n_targ_recomb)
  if (return_genos == T) {
    all_genos = matrix(NA, n_loci, n_ind)
  } else {
    all_genos = NULL
  }
  if (write_genos == T) {
    files = rep(NA, n_ind)
  } else {
    files = NULL
  }
  export_obj = c("createGeno", "createHap")
  foreach::foreach(i = 1:n_ind, .packages = c("matrixStats"), .export = c(export_obj)) %dopar% {
    try({
      adm_ind = adm_mat[i,]
      start_geno = start_genos[i]
      recomb_sites1 = sort(unique(all_recomb_sites[recomb_index == i*2-1]))
      recomb_sites2 = sort(unique(all_recomb_sites[recomb_index == i*2]))
      geno_long = createGeno(dis_chr, dir, n_loci, adm_ind, start_geno, start_pos, maj_sample_df, min_sample_df, all_sample_df, read_type = "rds",
                             hap_vec_mat = NULL, recomb_sites1 = recomb_sites1, recomb_sites2 = recomb_sites2, recomb_prob = NULL, max_ma_count = max_ma_count)
      if(return_genos == T) {
        all_genos[,i] = geno_long
      }
      if(write_genos == T) {
        max_char = nchar(as.character(n_ind))
        nchar_i = nchar(as.character(i))
        if (nchar_i < max_char) {
          prepend_zeros = paste0(rep("0", max_char-nchar_i), collapse = "")
          affix_i = paste0(prepend_zeros, as.character(i))
        } else {
          affix_i = as.character(i)
        }
        file_i = paste0(paste(prefix, chrom_tag, affix_i, sep = "_"), ".rds")
        if (is.null(loci_out)) {
          saveRDS(geno_long, paste0(dir,file_i), compress = compress)
        } else {
          saveRDS(geno_long[loci_out], paste0(dir,file_i))
        }
      }
    })
  }
  if (write_genos == T) {
    for (i in 1:n_ind) {
      max_char = nchar(as.character(n_ind))
      nchar_i = nchar(as.character(i))
      if(nchar_i < max_char) {
        prepend_zeros = paste0(rep("0", max_char-nchar_i), collapse = "")
        affix_i = paste0(prepend_zeros, as.character(i))
      } else {
        affix_i = as.character(i)
      }
      file_i = paste0(paste(prefix, chrom_tag, affix_i, sep = "_"), ".rds")
      files[i] = file_i
    }
  }
  out = list(genos = all_genos, files = files)
  return(out)
}

createGenoMult = function(dis_chr = "no", dir, loci_info, sample_ma_list, adm_mat, start_genos, snp_map, chr_recomb, hap_select_probs = NULL,
                          prefix = "geno", write_genos = T, return_genos = F, write_every = 5, n_cores = 4, compress = F) {
  n_snp = ncol(start_genos)
  chr_vec = as.numeric(names(loci_info))
  n_chr = length(chr_vec)
  genoList = vector(mode = "list", length = n_chr)
  lociList = vector(mode = "list", length = n_chr)
  fileList = vector(mode = "list", length = n_chr)
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  for (i in 1:n_chr) {
    chr_i = chr_vec[i]
    snp_i = snp_map[snp_map[, "chr"] == chr_i, "snp"][1]
    if (dis_chr[i] == "yes") {
      start_pos_i = snp_map[snp_map[, "chr"] == chr_i, "order"][1]
    } else {
      start_pos_i = sample(1:nrow(loci_info[[i]]), 1)
    }
    max_ma_count_i = max(sample_ma_list[[i]][, "ma_count"])
    chr_sample_df_i = sample_ma_list[[i]][, 1:5]
    if (is.null(hap_select_probs)) {
      chr_sample_df_i[, "weight"] = 1
    } else {
      hap_tag = gsub(".rds", "", chr_sample_df_i[, "file"], fixed = T)
      hap_tag = gsub(paste0("chr", chr_i, "_"), "", hap_tag, fixed = T)
      chr_sample_df_i[, "hap_tag"] = hap_tag
      chr_sample_df_i[, "weight"] = hap_select_probs[match(hap_tag, hap_select_probs[, "hap_tag"]), "weight"]
    }
    colnames(chr_sample_df_i)[5] = "min_allele"
    recomb_prob_i = chr_recomb[[i]]
    n_loci_i = length(recomb_prob_i)+1
    loci_range_i = 1:n_loci_i
    loci_sample_i = loci_range_i[loci_range_i %% write_every == 0]
    loci_out_i = sort(unique(c(start_pos_i, loci_sample_i)))
    names(loci_out_i) = loci_info[[i]][loci_out_i, "name"]
    cG_out = createGeno_W(dis_chr = dis_chr[i], dir = dir, chr_sample_df = chr_sample_df_i, adm_mat = adm_mat,
                          start_genos = start_genos[, snp_i], start_pos = start_pos_i, recomb_prob = recomb_prob_i, max_ma_count = max_ma_count_i, chr = chr_i,
                          prefix = prefix, write_genos = write_genos, return_genos = return_genos, loci_out = loci_out_i, compress = compress)
    geno_i = cG_out[["genos"]]
    files_i = cG_out[["files"]]
    genoList[[i]] = geno_i
    lociList[[i]] = loci_out_i
    fileList[[i]] = files_i
  }
  parallel::stopCluster(cl)
  out = list(genoList = genoList, lociList = lociList, fileList = fileList)
  return(out)
}

snpWindows = function(snp_pos, n_loci, windowSize, div = 1) {
  n_snp = length(snp_pos)
  windows_end = rep(NA, n_snp)
  windows_start = rep(NA, n_snp)
  ranges_list = vector(mode = "list", length = n_snp)
  window_num_list = vector(mode = "list", length = n_snp)
  for (i in 1:n_snp){
    windows_start[i] = max(snp_pos[i]-trunc(windowSize/2), 1)
    windows_end[i]= min(snp_pos[i]+trunc(windowSize/2), n_loci)
    window_range_i = windows_start[i]:windows_end[i]
    ranges_list[[i]] = sort(unique(c(snp_pos[i], window_range_i[window_range_i %% div == 0])))
    window_num_list[[i]] = rep(i, length(ranges_list[[i]]))
  }
  rangesOut = unlist(ranges_list)
  window_num = unlist(window_num_list)
  out_df = cbind(window = window_num, locus_pos = rangesOut)
  return(out_df)
}

snpWindows_W = function(chr_vec, snp, loci_info, windowSize, div = 1) {
  n_chr = length(chr_vec)
  n_snp = length(snp)
  out = as.data.frame(matrix(NA, nrow = (windowSize+2)*n_snp, ncol = 5))
  out_cols = c("chr", "window", "name", "locus_pos", "center_locus")
  colnames(out) = out_cols
  a = 0
  b = 0
  max_window = 0
  for (h in 1:n_chr){
    loci_info_h = loci_info[[h]]
    snp_h = snp[snp %in% loci_info_h[, "name"]]
    snp_pos_h = sort(match(snp_h, loci_info_h[, "name"]))
    snp_h = loci_info_h[snp_pos_h, "name"]
    out_h = snpWindows(snp_pos = snp_pos_h, n_loci = nrow(loci_info_h), windowSize = windowSize, div = div)
    a = b + 1
    b = b + nrow(out_h)
    out[a:b, "center_locus"] = snp_h[out_h[, "window"]]
    out[a:b, "chr"] = chr_vec[h]
    out[a:b, "name"] = loci_info_h[out_h[, "locus_pos"], "name"]
    out[a:b, "locus_pos"] = out_h[, "locus_pos"]
    out[a:b, "window"] = out_h[, "window"] + max_window
    max_window = max(out[a:b, "window"])
  }
  out = out[is.na(out[, "window"]) == F,]
  return(out)
}

make_inSampleList = function(in_df) {
  chr_vec = sort(unique(in_df[, "chr"]))
  n_chr = length(chr_vec)
  out = vector(mode = "list", length = n_chr)
  for (h in 1:n_chr){
    out[[h]] = in_df[in_df[, "chr"] == chr_vec[h], "locus_pos"]
    names(out[[h]]) = in_df[in_df[, "chr"] == chr_vec[h], "name"]
  }
  return(out)
}

subsettedGenoMat = function(dir, obj_cGM, inSampleList = NULL) {
  lociList = obj_cGM[["lociList"]]
  fileList = obj_cGM[["fileList"]]
  n_chr = length(lociList)
  n_ind = length(fileList[[1]])
  if(is.null(inSampleList)) {
    inSampleList = obj_cGM[["lociList"]]
  }
  n_sample = rep(NA, n_chr)
  for (h in 1:n_chr) {
    n_sample[h] = length(inSampleList[[h]])
  }
  sampleList = vector(mode = "list", length = n_chr)
  genoMat = matrix(NA, nrow = n_ind, ncol = sum(n_sample))
  colnames(genoMat) = paste0("v", 1:ncol(genoMat))
  a = 0
  b = 0
  for (h in 1:n_chr) {
    a = b + 1
    b = b + n_sample[h]
    length_h = length(lociList[[h]])
    sample_h = inSampleList[[h]]
    sampleList[[h]] = sample_h
    colnames(genoMat)[a:b] = names(sample_h)
    for (i in 1:n_ind){
      hap_h_i = readRDS(paste0(dir,fileList[[h]][i]))
      genoMat[i, a:b] = hap_h_i[sample_h]
    }
  }
  out = list(genoMat = genoMat, sampleList = sampleList)
}

gwas = function(gwas_df, y_name = "y", cov_names = NULL, s_names, type, n_cores = 1) {
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  n_gwas_loci = length(s_names)
  n_cov = length(cov_names)
  gwas_p = foreach::foreach (n = 1:n_gwas_loci, .combine = c) %dopar% {
    gwas_formula_n = paste0(y_name, " ~ ", paste(c(cov_names, s_names[n]), sep = "", collapse = " + "))
    if (type == "logit") {
      glm_n = stats::glm(formula = stats::as.formula(gwas_formula_n), data = gwas_df,
                         family = stats::binomial(link = "logit"))
    }
    if (type == "quant") {
      glm_n = stats::lm(formula = stats::as.formula(gwas_formula_n), data = gwas_df)
    }
    summary_n = summary(glm_n)$coefficients
    out_n = 1
    if(nrow(summary_n) > n_cov + 1) {
      out_n = summary_n[n_cov + 2, 4]
    }
    return(out_n)
  }
  parallel::stopCluster(cl)
  names(gwas_p) = s_names
  return(gwas_p)
}

calcLD = function(center_locus, quant_sub, windows) {
  n_snp = length(unique(center_locus))
  n_gwas_loci = nrow(quant_sub)
  corrs = rep(NA, n_gwas_loci)
  n_obs = ncol(quant_sub)
  names(corrs) = rownames(quant_sub)
  a = 0
  b = 0
  for (i in 1:n_snp){
    snp_i = center_locus[windows == i][1]
    a = b + 1
    n_corr_i = sum(windows == i)
    b = b + n_corr_i
    snp_i_mat = matrix(quant_sub[snp_i,], n_corr_i, n_obs, byrow = T)
    covs_i_mat = (quant_sub[a:b,] - rowMeans(quant_sub[a:b,])) * ((snp_i_mat) - rowMeans(snp_i_mat))
    corrs_i = rowSums(covs_i_mat)/(n_obs-1)/(sqrt(matrixStats::rowVars(quant_sub[a:b,]))*sqrt(matrixStats::rowVars(snp_i_mat)))
    corrs[a:b] = corrs_i
  }
  return(corrs)
}

makeSpacedGroups = function(values, group, spaceSize = 10) {
  u_groups = unique(group)
  n_groups = length(u_groups)
  n_values = length(values)
  values_spaced = rep(NA, n_values + (n_groups-1)*spaceSize)
  a = 0
  b = 0
  for (i in 1:n_groups){
    a = b + 1
    s = a + (i-1)*spaceSize
    n_values_i = sum(group == u_groups[i])
    b = b + n_values_i
    t = b + (i-1)*spaceSize
    values_spaced[s:t] = values[group == u_groups[i]]
    names(values_spaced)[s:t] = names(values[group == u_groups[i]])
  }
  return(values_spaced)
}

plot_gwas = function(p_values, denom, y_axe_max = NULL, cols, width = 1050, height = 350,
                     mar = c(3,3,.5,1), mgp = c(2,1,0), pch = 19, ps = 12, bty = "l", cex = 1, cex.axis = 1, cex.lab = 1,
                     line_col = "grey", lwd = 2, lty = 2, legend_pos = "topright", box.lty = 0) {
  if(is.null(y_axe_max)) {
    max_p = max(-log10(p_values), na.rm = T)
    y_axe_max = (trunc(max_p/10)+1)*10
  }
  graphics::par(mar = mar, ps = ps)
  plot(-log10(p_values), xlab = "", ylab = "-log10(p-values)", xaxt = "n", pch = pch,
       cex = cex, ylim = c(0, y_axe_max), col = cols, bty = bty,
       cex.axis = 1, cex.lab = cex.lab, mgp = mgp)
  graphics::abline(h = -log10(.05/ denom), col = line_col, lwd = lwd, lty = lty)
  graphics::legend(legend_pos, legend=c("1.0", "[0.8, 1.0)", "[0.6, 0.8)", "[0.4, 0.6)", "[0.2, 0.4)", "[0.0, 0.2)"),
         col=c("magenta", "red", "orange", "green", "cyan", "blue"), pch = pch, cex=cex,
         title="R-Sq.", box.lty = box.lty, bg = "transparent")
}

selectSNP = function(loci_info) {
  n_chr = length(loci_info)
  snp = rep(NA, n_chr)
  for (h in 1:n_chr) {
    n_loci_h = nrow(loci_info[[h]])
    pos_h = round(stats::median(1:n_loci_h), 0)
    snp[h] = loci_info[[h]][pos_h, "name"]
  }
  return(snp)
}

calcRV_weights = function(dir, chr_vec, maf_thresh = NULL, mean_effect, weight_type = "fixed") {
  n_chr = length(chr_vec)
  freq_list = vector(mode = "list", length = n_chr)
  for (h in 1:n_chr) {
    chr_h = chr_vec[h]
    chrom_tag = paste0("chr", chr_h)
    all_freq_file = paste0(chrom_tag, "_ALL_freqs.rds")
    infile = paste0(dir, all_freq_file)
    freqs_h = readRDS(infile)
    freq_list[[h]] = freqs_h
  }
  freqs = unlist(freq_list)
  maf = pmin(freqs, 1 - freqs)
  n_loci = length(freqs)
  if(weight_type == "MB") {
    weights = 1/sqrt(maf*(1-maf))
  } else {
    weights = as.numeric(maf < maf_thresh)
  }
  weights[is.finite(weights) == F] = 0
  weights_adj = n_loci*mean_effect*weights/sum(weights)
  weights_adj[is.nan(weights_adj)] = 0
  weights_adj_list = vector(mode = "list", length = n_chr)
  a = 0
  b = 0
  for (h in 1:n_chr) {
    n_loci_h = length(freq_list[[h]])
    a = b + 1
    b = b + n_loci_h
    weights_adj_list[[h]] = weights_adj[a:b]
  }
  out = list(freq_list = freq_list, weights = weights, weights_adj_list = weights_adj_list)
  return(out)
}

scoreRV_hap = function(dir, chr_vec, hap_tag, mean_effect, maf_thresh, freq_list, weights_adj_list) {
  n_chr = length(chr_vec)
  chrom_tag = paste0("chr", chr_vec)
  hap_files = paste0(dir, chrom_tag, "_", hap_tag, ".rds")
  hap_score = matrix(NA, nrow = 1, ncol = n_chr)
  hap_count = matrix(NA, nrow = 1, ncol = n_chr)
  for (h in 1:n_chr) {
    hap_file_h = hap_files[h]
    alts_h = readRDS(hap_file_h)
    freqs_h = freq_list[[h]]
    n_loci_h = length(freqs_h)
    weights_h = weights_adj_list[[h]]
    refs_h = (1:n_loci_h)[-alts_h]
    alt_rares_h = which(freqs_h <= .5)
    ref_rares_h = which(freqs_h > .5)
    hap_h = rep(0, n_loci_h)
    hap_h[intersect(ref_rares_h, refs_h)] = 1
    hap_h[intersect(alt_rares_h, alts_h)] = 1
    hap_score[h] = sum(hap_h*weights_h)
    hap_count[h] = sum(hap_h)
    hap_out = cbind(hap_score, hap_count)
  }
  return(hap_out)
}

scoreRV_haps = function(dir, chr_vec, sample_filename, maf_thresh, elig_pops, mean_effect, weight_type, n_cores) {
  n_chr = length(chr_vec)
  sample_info = utils::read.table(paste0(dir,sample_filename), sep = "\t", header = T, stringsAsFactors = FALSE)
  n_ind = nrow(sample_info)
  pops = sort(unique(sample_info[, "pop"]))
  pop_counts = table(sample_info[, "pop"])
  n_pops = length(pops)
  hap_tags = paste0(rep(sample_info[, "sample"], 2), "_", sort(rep(c(1, 2), each = n_ind)))
  hap_tags = sort(hap_tags)
  hap_scores = as.data.frame(matrix(NA, nrow = nrow(sample_info)*2, ncol = n_chr + 4))
  colnames(hap_scores) = c("pop", "hap_tag", "ind_name", "total", paste0("chr", chr_vec))
  hap_scores[, "hap_tag"] = hap_tags
  hap_scores[, "ind_name"] = substr(hap_tags, 1, nchar(hap_tags)-2)
  hap_scores[, "pop"] = sample_info[match(hap_scores[, "ind_name"], sample_info[, "sample"]), "pop"]
  hap_scores = hap_scores[hap_scores[,"pop"] %in% elig_pops,]
  hap_counts = hap_scores
  n_hap = nrow(hap_scores)
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  weightsOut = calcRV_weights(dir = dir, chr_vec = chr_vec,  maf_thresh = maf_thresh, mean_effect = mean_effect, weight_type = weight_type)
  weights_adj_list = weightsOut[["weights_adj_list"]]
  freq_list = weightsOut[["freq_list"]]
  exportObj = c("scoreRV_hap")
  hap_dopar = foreach::foreach (i = 1:n_hap, .combine = rbind, .export = exportObj) %dopar% {
    hap_tag_i = hap_scores[i, "hap_tag"]
    hap_score_i = scoreRV_hap(dir = dir, chr_vec = chr_vec, hap_tag = hap_tag_i, mean_effect = mean_effect, maf_thresh = maf_thresh, freq_list = freq_list,
                              weights_adj_list = weights_adj_list)
  }
  hap_scores[, -c(1:4)]= hap_dopar[, 1:n_chr]
  hap_scores[, "total"] = rowSums(hap_dopar[, 1:n_chr])
  hap_counts[, -c(1:4)] = hap_dopar[, -(1:n_chr)]
  hap_counts[, "total"] = rowSums(hap_dopar[, -c(1:n_chr)])
  ind_scores = as.data.frame(matrix(NA, nrow = n_hap/2, ncol = 3))
  colnames(ind_scores) = c("pop", "name", "total")
  ind_scores_agg = stats::aggregate(total ~ ind_name, hap_scores, FUN = sum)
  ind_scores[, c("name", "total")] = ind_scores_agg
  ind_scores[, "pop"] = sample_info[match(ind_scores[, "name"], sample_info[, "sample"]), "pop"]
  parallel::stopCluster(cl)
  return(out = list(hap_scores = hap_scores, hap_counts = hap_counts, ind_scores = ind_scores, weights_adj_list = weights_adj_list, freq_list = freq_list))
}

solveRV_B0 = function(B0Est, spec_val_k, ind_scores_k, score_col, type = type) {
  if(type == "bin") {
    z_k = B0Est + ind_scores_k[, score_col]
    ind_scores_k[, "val"] = exp(z_k)/(1+exp(z_k))
  } else {
    ind_scores_k[, "val"] = B0Est + ind_scores_k[, score_col]
  }
  diff = (mean(ind_scores_k[, "val"]) - spec_val_k)^2
  return(diff)
}

solveRV_B0_W = function(spec_val, ind_scores, score_col = "total", type = "bin", lower = -20, upper = 20) {
  pops = sort(unique(ind_scores[, "pop"]))
  n_pops = length(pops)
  out = ind_scores
  out[, c("B0Est", "pred_val")] = NA
  if(type == "bin") {
    out[, "case_share"] = NA
    out[, "cont_share"] = NA
  }
  for (k in 1:n_pops) {
    spec_val_k = spec_val[k]
    inds_k = ind_scores[, "pop"] == pops[k]
    ind_scores_k = ind_scores[inds_k,]
    optOut_k = stats::optimize(solveRV_B0, spec_val_k = spec_val_k, ind_scores_k = ind_scores_k, score_col = score_col, type = type, lower = lower, upper = upper)
    B0Est_k = optOut_k[["minimum"]]
    out[inds_k, "B0Est"] = B0Est_k
    if(type == "bin") {
      pred_z_k = out[inds_k, "B0Est"] + out[inds_k, "total"]
      out[inds_k, "pred_val"] = exp(pred_z_k)/(1+exp(pred_z_k))
      out[inds_k, "case_share"] = out[inds_k, "pred_val"]/sum(out[inds_k, "pred_val"])
      out[inds_k, "cont_share"] = (1-out[inds_k, "pred_val"])/sum(1-out[inds_k, "pred_val"])
    } else {
      out[inds_k, "pred_val"] = out[inds_k, "B0Est"] + out[inds_k, "total"]
    }
  }
  return(out)
}

calcRV_shares_cc = function(mean_effect, ind_probs, hap_scores, score_col = "total", prop_weighted = T) {
  hap_probs = hap_scores
  inds_match = match(hap_probs[, "ind_name"], ind_probs[, "name"])
  hap_probs[, "ind_total"] = ind_probs[inds_match, score_col]
  hap_probs[, c("ind_case_share", "ind_cont_share")] = ind_probs[inds_match, c("case_share", "cont_share")]
  if(prop_weighted == T) {
    if(sign(mean_effect) >= 1) {
      hap_probs[, c("case_share")] = (hap_probs[, "ind_case_share"] * hap_probs[, "total"] / hap_probs[, "ind_total"])
      hap_probs[, c("cont_share")] = (hap_probs[, "ind_cont_share"] * (1 - hap_probs[, "total"] / hap_probs[, "ind_total"]))
    } else {
      hap_probs[, c("case_share")] = (hap_probs[, "ind_case_share"] * (1-hap_probs[, "total"] / hap_probs[, "ind_total"]))
      hap_probs[, c("cont_share")] = (hap_probs[, "ind_cont_share"] * hap_probs[, "total"] / hap_probs[, "ind_total"])
    }
  } else {
    hap_probs[, c("case_share", "cont_share")] = hap_probs[, c("ind_case_share", "ind_cont_share")] / 2
  }
  return(hap_probs)
}

scoreRV_geno = function(dir, file_list_h, weights_adj_list_h, freq_list_h) {
  n_inds = length(file_list_h)
  ref_rares_h = freq_list_h > .5
  out = foreach::foreach (i = 1:n_inds, .combine = rbind) %dopar% {
    file_i = file_list_h[i]
    geno_i = readRDS(paste0(dir,file_i))
    rares_i = geno_i
    rares_i[ref_rares_h] = 2 - geno_i[ref_rares_h]
    score_i = sum(rares_i*weights_adj_list_h)
    count_i = sum(rares_i)
    out_i = c(score_i, count_i)
  }
  return(out)
}

scoreRV_genos = function(dir, file_list, weights_adj_list, freq_list, adm_mat, B0Est, type = "bin", n_cores = 1) {
  cl = parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  elig_pops = names(B0Est)
  n_chr = length(file_list)
  for (h in 1:n_chr) {
    file_list_h = file_list[[h]]
    weights_adj_list_h = weights_adj_list[[h]]
    freq_list_h = freq_list[[h]]
    if(h == 1) {
      n_inds = length(file_list_h)
      scoreMat = matrix(NA, n_inds, n_chr)
      countMat = matrix(NA, n_inds, n_chr)
    }
    out_h = scoreRV_geno(dir = dir, file_list_h = file_list_h, weights_adj_list_h = weights_adj_list_h, freq_list_h = freq_list_h)
    scoreMat[,h] = out_h[,1]
    countMat[,h] = out_h[,2]
  }
  scoreTotal = adm_mat[,elig_pops]%*%B0Est + rowSums(scoreMat)
  if(type == "bin") {
    scoreOut = exp(scoreTotal)/(1 + exp(scoreTotal))
  } else {
    scoreOut = scoreTotal
  }
  parallel::stopCluster(cl)
  countTotal = rowSums(countMat)
  out = cbind(scoreOut, countTotal)
  colnames(out) = c("score", "count")
  return(out)
}
