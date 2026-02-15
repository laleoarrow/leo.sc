.onAttach <- function(libname, pkgname) {
  logo_core <- c(
    "  \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2557          \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D \u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557         \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u2550\u255D  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557   \u2588\u2588\u2551   \u2588\u2588\u2551   \u2588\u2588\u2557   \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2551        ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255D   \u2588\u2588\u2551   \u2588\u2588\u2551   \u255A\u2550\u255D   \u255A\u2550\u2550\u2550\u2550\u2588\u2588\u2551 \u2588\u2588\u2551        ",
    "  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D         \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2551 \u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D\u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u255D          \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D   "
  )

  device_df <- tryCatch(.collect_device_info(), error = function(e) NULL)
  if (is.null(device_df) || !all(c("device_line", "ram_line") %in% names(device_df))) {
    device_df <- data.frame(device_line = "Device: Unknown", ram_line = "Memory: NA", stringsAsFactors = FALSE)
  }

  meta_core <- c(
    as.character(device_df$device_line[1]),
    as.character(device_df$ram_line[1]),
    paste0("Welcome to leo.sc package (v", utils::packageVersion("leo.sc"), ")"),
    paste0("System Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Have fun with single cell data"
  )

  meta_core <- meta_core[nzchar(meta_core)]
  meta_width <- max(nchar(meta_core, type = "width"))
  meta <- character(length(meta_core))
  for (i in seq_along(meta_core)) {
    pad <- max(0, meta_width - nchar(meta_core[i], type = "width"))
    left_pad <- floor(pad / 2)
    meta[i] <- paste0(">>> ", strrep(" ", left_pad), meta_core[i], strrep(" ", pad - left_pad), " <<<")
  }

  max_width <- max(nchar(c(logo_core, meta), type = "width"))
  border <- strrep("=", max_width)
  full_content <- c(border)
  for (i in seq_along(logo_core)) {
    pad <- max(0, max_width - nchar(logo_core[i], type = "width"))
    full_content <- c(full_content, paste0(logo_core[i], strrep(" ", pad)))
  }
  full_content <- c(full_content, border)
  for (i in seq_along(meta)) {
    pad <- max(0, max_width - nchar(meta[i], type = "width"))
    left_pad <- floor(pad / 2)
    full_content <- c(full_content, paste0(strrep(" ", left_pad), meta[i], strrep(" ", pad - left_pad)))
  }

  color_ready <- nzchar(Sys.getenv("RSTUDIO")) || Sys.getenv("COLORTERM") == "truecolor" || Sys.getenv("TERM") != "dumb"
  if (color_ready) {
    packageStartupMessage(.render_vignette_physical(full_content, 2, 10, max_width))
  } else {
    packageStartupMessage(paste(full_content, collapse = "\n"))
  }
}

.render_vignette_physical <- function(lines,
                                      vertical_padding,
                                      horizontal_padding,
                                      core_width,
                                      fg_start = "#e82020",
                                      fg_mid = "#FFFFFF",
                                      fg_end = "#2626e1",
                                      bg_center = "#000000",
                                      bg_halo = "#523111") {
  ext_lines <- c(rep(strrep(" ", core_width), vertical_padding), lines, rep(strrep(" ", core_width), vertical_padding))
  row_count <- length(ext_lines)
  max_dist <- sqrt(vertical_padding ^ 2 + horizontal_padding ^ 2)
  fg_rgb <- grDevices::col2rgb(grDevices::colorRampPalette(c(fg_start, fg_mid, fg_end))(core_width))
  bg_center_rgb <- grDevices::col2rgb(bg_center)
  bg_halo_rgb <- grDevices::col2rgb(bg_halo)
  out <- character(row_count)

  for (row_index in seq_along(ext_lines)) {
    chars <- strsplit(ext_lines[row_index], "", fixed = TRUE)[[1]]
    dist_y <- max(0, vertical_padding + 1 - row_index, row_index - (row_count - vertical_padding))
    dist_y2 <- dist_y ^ 2
    mix_y <- if (max_dist > 0) (dist_y / max_dist) ^ 1.4 else 0
    row_bg <- (1 - mix_y) * bg_center_rgb + mix_y * bg_halo_rgb
    left <- character(horizontal_padding)
    right <- character(horizontal_padding)
    core <- character(length(chars))

    for (col_index in seq_len(horizontal_padding)) {
      left_dist <- sqrt((horizontal_padding + 1 - col_index) ^ 2 + dist_y2)
      right_dist <- sqrt(col_index ^ 2 + dist_y2)
      if (left_dist <= max_dist) {
        left_w <- (left_dist / max_dist) ^ 1.4
        left_bg <- (1 - left_w) * bg_center_rgb + left_w * bg_halo_rgb
        left[col_index] <- sprintf("\033[38;2;0;0;0m\033[48;2;%d;%d;%dm ", as.integer(left_bg[1]), as.integer(left_bg[2]), as.integer(left_bg[3]))
      } else {
        left[col_index] <- " "
      }
      if (right_dist <= max_dist) {
        right_w <- (right_dist / max_dist) ^ 1.4
        right_bg <- (1 - right_w) * bg_center_rgb + right_w * bg_halo_rgb
        right[col_index] <- sprintf("\033[38;2;0;0;0m\033[48;2;%d;%d;%dm ", as.integer(right_bg[1]), as.integer(right_bg[2]), as.integer(right_bg[3]))
      } else {
        right[col_index] <- " "
      }
    }

    visual_pos <- 1
    for (i in seq_along(chars)) {
      fg_idx <- min(core_width, max(1, visual_pos))
      core[i] <- sprintf(
        "\033[48;2;%d;%d;%dm\033[38;2;%d;%d;%dm%s",
        as.integer(row_bg[1]),
        as.integer(row_bg[2]),
        as.integer(row_bg[3]),
        fg_rgb[1, fg_idx],
        fg_rgb[2, fg_idx],
        fg_rgb[3, fg_idx],
        chars[i]
      )
      visual_pos <- visual_pos + nchar(chars[i], type = "width")
    }

    out[row_index] <- paste0(paste0(left, collapse = ""), paste0(core, collapse = ""), paste0(right, collapse = ""), "\033[0m")
  }

  paste(out, collapse = "\n")
}

.collect_device_info <- function() {
  sys_info <- tryCatch(Sys.info(), error = function(e) NULL)
  system_name <- tolower(if (!is.null(sys_info) && !is.na(sys_info[["sysname"]])) sys_info[["sysname"]] else .Platform$OS.type)
  os_release <- if (!is.null(sys_info) && !is.na(sys_info[["release"]])) sys_info[["release"]] else ""
  os_name <- switch(system_name, darwin = "macOS", linux = "Linux", windows = "Windows", system_name)
  os_label <- trimws(paste(os_name, os_release))

  cmd_out <- function(cmd, args = character()) {
    out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE))
    if (length(out) == 0) return(character())
    out[!is.na(out)]
  }
  num1 <- function(x) {
    if (length(x) == 0) return(NA_real_)
    suppressWarnings(as.numeric(x[1]))
  }
  pick_num <- function(lines, pattern) {
    hit <- lines[grepl(pattern, lines)][1]
    if (length(hit) == 0 || is.na(hit)) return(NA_real_)
    suppressWarnings(as.numeric(gsub("[^0-9]", "", hit)))
  }

  cpu_name <- ""
  memory_total_gb <- NA_real_
  memory_free_gb <- NA_real_

  if (system_name == "darwin") {
    cpu_name <- paste(cmd_out("sysctl", c("-n", "machdep.cpu.brand_string")), collapse = " ")
    total_bytes <- num1(cmd_out("sysctl", c("-n", "hw.memsize")))
    if (is.finite(total_bytes)) memory_total_gb <- total_bytes / 1024^3
    vm_lines <- cmd_out("vm_stat")
    page_line <- vm_lines[grepl("page size of", vm_lines)][1]
    page_size <- suppressWarnings(as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", page_line)))
    free_pages <- pick_num(vm_lines, "^Pages free")
    spec_pages <- pick_num(vm_lines, "^Pages speculative")
    inactive_pages <- pick_num(vm_lines, "^Pages inactive")
    if (is.finite(page_size)) {
      avail_pages <- sum(c(free_pages, spec_pages, inactive_pages), na.rm = TRUE)
      if (is.finite(avail_pages)) memory_free_gb <- (avail_pages * page_size) / 1024^3
    }
  }

  if (system_name == "linux") {
    if (file.exists("/proc/cpuinfo")) {
      cpu_lines <- readLines("/proc/cpuinfo", warn = FALSE)
      cpu_hit <- cpu_lines[grepl("^model name\\s*:", cpu_lines)][1]
      if (!is.na(cpu_hit)) cpu_name <- sub("^[^:]*:\\s*", "", cpu_hit)
    }
    if (file.exists("/proc/meminfo")) {
      mem_lines <- readLines("/proc/meminfo", warn = FALSE)
      memory_total_kb <- pick_num(mem_lines, "^MemTotal\\s*:")
      memory_free_kb <- pick_num(mem_lines, "^MemAvailable\\s*:")
      if (!is.finite(memory_free_kb)) memory_free_kb <- pick_num(mem_lines, "^MemFree\\s*:")
      if (is.finite(memory_total_kb)) memory_total_gb <- memory_total_kb / 1024^2
      if (is.finite(memory_free_kb)) memory_free_gb <- memory_free_kb / 1024^2
    }
  }

  if (system_name == "windows") {
    cpu_lines <- cmd_out("wmic", c("cpu", "get", "Name", "/value"))
    cpu_hit <- cpu_lines[grepl("^Name=", cpu_lines)][1]
    if (!is.na(cpu_hit)) cpu_name <- sub("^Name=", "", cpu_hit)
    mem_lines <- cmd_out("wmic", c("OS", "get", "FreePhysicalMemory,TotalVisibleMemorySize", "/value"))
    total_kb <- num1(sub("^TotalVisibleMemorySize=", "", mem_lines[grepl("^TotalVisibleMemorySize=", mem_lines)][1]))
    free_kb <- num1(sub("^FreePhysicalMemory=", "", mem_lines[grepl("^FreePhysicalMemory=", mem_lines)][1]))
    if (is.finite(total_kb)) memory_total_gb <- total_kb / 1024^2
    if (is.finite(free_kb)) memory_free_gb <- free_kb / 1024^2
  }

  cpu_name <- gsub("\\s+", " ", trimws(cpu_name))
  if (!nzchar(cpu_name) && !is.null(sys_info) && !is.na(sys_info[["machine"]])) cpu_name <- sys_info[["machine"]]
  if (!nzchar(cpu_name)) cpu_name <- "CPU Unknown"
  if (nchar(cpu_name, type = "width") > 34) cpu_name <- paste0(substr(cpu_name, 1, 31), "...")

  core_physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  core_logical <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(core_physical) && system_name == "darwin") core_physical <- num1(cmd_out("sysctl", c("-n", "hw.physicalcpu")))
  if (!is.finite(core_logical) && system_name == "darwin") core_logical <- num1(cmd_out("sysctl", c("-n", "hw.logicalcpu")))
  if (!is.finite(core_logical) && system_name == "linux") core_logical <- num1(cmd_out("nproc"))
  if (!is.finite(core_logical)) core_logical <- suppressWarnings(as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS")))
  if (!is.finite(core_physical) && is.finite(core_logical)) core_physical <- core_logical
  core_label <- if (is.finite(core_physical) && is.finite(core_logical)) paste0(core_physical, "C/", core_logical, "T") else if (is.finite(core_logical)) paste0(core_logical, " cores") else "cores NA"

  memory_used_gb <- if (is.finite(memory_total_gb) && is.finite(memory_free_gb)) max(0, memory_total_gb - memory_free_gb) else NA_real_
  total_ram <- if (is.finite(memory_total_gb)) format(round(memory_total_gb, 1), nsmall = 1) else "NA"
  used_ram <- if (is.finite(memory_used_gb)) format(round(memory_used_gb, 1), nsmall = 1) else "NA"
  free_ram <- if (is.finite(memory_free_gb)) format(round(memory_free_gb, 1), nsmall = 1) else "NA"
  ram_line <- if (is.finite(memory_total_gb)) paste0("Memory: ", total_ram, "G total | ", used_ram, "G used | ", free_ram, "G free") else paste0("Memory: ", free_ram, "G free | total/used NA")

  data.frame(device_line = paste0("Device: ", os_label, " | ", cpu_name, " | ", core_label), ram_line = ram_line, stringsAsFactors = FALSE)
}
