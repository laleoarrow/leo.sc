
.onAttach <- function(libname, pkgname) {
  # Based on user-provided pixel-perfect draft
  # Translation: █=\u2588, ║=\u2551, ╗=\u2557, ═=\u2550, ╔=\u2554, ╝=\u255d, ╚=\u255a
  banner <- paste(
    "======================================================",
    "  \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557          \u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557         \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551   \u2588\u2588\u2551   \u2588\u2588\u2557   \u255A\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2551       ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255D  \u2588\u2588\u2551   \u2588\u2588\u2551   \u255A\u2550\u255D    \u255A\u2550\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2551       ",
    "  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D         \u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D  ",
    "  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D\u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D \u255A\u2550\u2550\u2550\u2550\u2550\u255D          \u255A\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u255D  ",
    "======================================================",
    sep = "\n"
  )
  
  pkg_version <- utils::packageVersion("leo.sc")
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  meta_text <- c(
    paste0(">>> Welcome to leo.sc package (v", pkg_version, ") <<<"),
    paste0(">>> System Time: ", current_time, " <<<"),
    ">>> Have fun with single cell data <<<"
  ) 
  
  # Padding for Vignette Halo
  h_pad <- 8; v_pad <- 2
  
  banner_lines <- strsplit(banner, "\n")[[1]]
  inner_width  <- max(nchar(banner_lines))
  full_width   <- inner_width + 2 * h_pad
  
  # Build center content with horizontal padding
  core_lines <- c(
    vapply(banner_lines, function(l) paste0(strrep(" ", h_pad), l, strrep(" ", inner_width - nchar(l)), strrep(" ", h_pad)), character(1)),
    vapply(meta_text, function(t) {
      n_t <- nchar(t)
      pl <- floor((inner_width - n_t) / 2)
      pr <- inner_width - n_t - pl
      paste0(strrep(" ", h_pad), strrep(" ", pl), t, strrep(" ", pr), strrep(" ", h_pad))
    }, character(1))
  )
  
  # Final full message with top/bottom buffers
  full_msg_lines <- c(rep(strrep(" ", full_width), v_pad), core_lines, rep(strrep(" ", full_width), v_pad))
  
  if (nzchar(Sys.getenv("RSTUDIO")) || .has_color()) {
     packageStartupMessage(.halo_render(full_msg_lines, v_pad, h_pad))
  } else {
    packageStartupMessage(paste(full_msg_lines, collapse = "\n"))
  }
}

.has_color <- function() {
  if (Sys.getenv("COLORTERM") == "truecolor") return(TRUE)
  if (Sys.getenv("TERM") == "dumb") return(FALSE)
  TRUE
}

.halo_render <- function(lines, v_p, h_p, 
                        fg_start = "#e82020", fg_mid = "#FFFFFF", fg_end = "#2626e1",
                        bg_center = "#000000", bg_halo = "#4d3422") {
  
  rows <- length(lines); cols <- nchar(lines[1])
  r_start <- v_p + 1; r_end <- rows - v_p
  c_start <- h_p + 1; c_end <- cols - h_p
  max_dist <- sqrt(v_p^2 + h_p^2)
  
  vapply(seq_along(lines), function(r) {
    chars <- strsplit(lines[r], "")[[1]]
    fgs <- grDevices::colorRampPalette(c(fg_start, fg_mid, fg_end))(length(chars))
    fg_rgb <- grDevices::col2rgb(fgs)
    
    res <- vapply(seq_along(chars), function(c) {
      dx <- max(0, c_start - c, c - c_end)
      dy <- max(0, r_start - r, r - r_end)
      dist <- sqrt(dx^2 + dy^2)
      
      if (dist == 0) {
        bg_rgb <- grDevices::col2rgb(bg_center)
      } else if (dist <= max_dist) {
        w <- dist / max_dist
        bg_rgb <- (1-w) * grDevices::col2rgb(bg_center) + w * grDevices::col2rgb(bg_halo)
      } else {
        return(chars[c]) # Outside glow
      }
      sprintf("\033[48;2;%d;%d;%dm\033[38;2;%d;%d;%dm%s", 
              as.integer(bg_rgb[1]), as.integer(bg_rgb[2]), as.integer(bg_rgb[3]),
              fg_rgb[1,c], fg_rgb[2,c], fg_rgb[3,c], chars[c])
    }, character(1))
    paste0(paste(res, collapse = ""), "\033[0m")
  }, character(1)) |> paste(collapse = "\n")
}
