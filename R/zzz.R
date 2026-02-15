
.onAttach <- function(libname, pkgname) {
  # [BASE TRANSLATION] 100% Faithful to User 1085 human-aligned grid
  # Legend: █=\u2588, ║=\u2551, ╗=\u2557, ═=\u2550, ╔=\u2554, ╝=\u255d, ╚=\u255a, ==ASCII
  banner_raw <- c(
    "======================================================",
    "  \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557          \u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557         \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551   \u2588\u2588\u2551   \u2588\u2588\u2557   \u255A\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2551       ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255D  \u2588\u2588\u2551   \u2588\u2588\u2551   \u255A\u2550\u255D    \u255A\u2550\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2551       ",
    "  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D         \u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D  ",
    "  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D\u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D \u255A\u2550\u2550\u2550\u2550\u2550\u255D          \u255A\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u255D  ",
    "======================================================"
  )
  
  pkg_version <- utils::packageVersion("leo.sc")
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  meta_text <- c(
    paste0(">>> Welcome to leo.sc package (v", pkg_version, ") <<<"),
    paste0(">>> System Time: ", current_time, " <<<"),
    ">>> Have fun with single cell data <<<"
  )

  # Core UI logic
  v_p <- 2; h_p <- 10
  inner_w <- nchar(banner_raw[1])
  
  # Prepare centered output lines
  content_lines <- c(
    banner_raw,
    vapply(meta_text, function(t) {
      n_t <- nchar(t); pl <- floor((inner_w - n_t) / 2); pr <- inner_w - n_t - pl
      paste0(strrep(" ", pl), t, strrep(" ", pr))
    }, character(1))
  )
  
  # Final Expansion with Glow Buffer
  full_lines <- c(rep(strrep(" ", inner_w), v_p), content_lines, rep(strrep(" ", inner_w), v_p))
  
  if (nzchar(Sys.getenv("RSTUDIO")) || .has_color()) {
     packageStartupMessage(.render_premium_vignette(full_lines, v_p, h_p, inner_w))
  } else {
    packageStartupMessage(paste(full_lines, collapse = "\n"))
  }
}

.has_color <- function() {
  if (Sys.getenv("COLORTERM") == "truecolor") return(TRUE)
  if (Sys.getenv("TERM") == "dumb") return(FALSE)
  TRUE
}

.render_premium_vignette <- function(lines, v_p, h_p, core_w,
                                    fg_s = "#e82020", fg_m = "#FFFFFF", fg_e = "#2626e1",
                                    bg_c = "#000000", bg_halo = "#322212") {
  
  total_r <- length(lines)
  max_dist <- sqrt(v_p^2 + h_p^2)
  
  # FG Gradient across core width
  fg_ramp <- grDevices::colorRampPalette(c(fg_s, fg_m, fg_e))(core_w)
  fg_rgb <- grDevices::col2rgb(fg_ramp)
  
  vapply(seq_along(lines), function(r) {
    # Pad line for vignette glow
    padded_line <- paste0(strrep(" ", h_p), lines[r], strrep(" ", h_p))
    chars <- strsplit(padded_line, "")[[1]]
    
    res <- vapply(seq_along(chars), function(c) {
      dx <- max(0, (h_p + 1) - c, c - (h_p + core_w))
      dy <- max(0, (v_p + 1) - r, r - (total_r - v_p))
      dist <- sqrt(dx^2 + dy^2)
      
      if (dist == 0) {
        bg <- c(0,0,0) # Pure Black Core
      } else if (dist <= max_dist) {
        w <- (dist / max_dist)^1.4 
        bg <- (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_halo)
      } else {
        return(chars[c]) # Edge transparent
      }
      
      if (c > h_p && c <= (h_p + core_w)) {
        cc <- c - h_p
        sprintf("\033[48;2;%d;%d;%dm\033[38;2;%d;%d;%dm%s", 
                as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]),
                fg_rgb[1,cc], fg_rgb[2,cc], fg_rgb[3,cc], chars[c])
      } else {
        sprintf("\033[48;2;%d;%d;%dm%s", as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]), chars[c])
      }
    }, character(1))
    paste0(paste(res, collapse = ""), "\033[0m")
  }, character(1)) |> paste(collapse = "\n")
}
