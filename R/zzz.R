
.onAttach <- function(libname, pkgname) {
  # [BASE] User 1085 manual layout
  # Legend: █=\u2588, ║=\u2551, ╗=\u2557, ╔=\u2554, ╝=\u255d, ╚=\u255a, ═=\u2550
  logo_core <- c(
    "  \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557          \u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557         \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551   \u2588\u2588\u2551   \u2588\u2588\u2557   \u255A\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2551       ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255D  \u2588\u2588\u2551   \u2588\u2588\u2551   \u255A\u2550\u255D    \u255A\u2550\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2551       ",
    "  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D         \u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D  ",
    "  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D\u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D \u255A\u2550\u2550\u2550\u2550\u2550\u255D          \u255A\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u255D  "
  )

  # Function to calculate visual/physical width (Unicode=2, ASCII=1)
  get_vw <- function(x) {
    chars <- strsplit(x, "")[[1]]
    sum(vapply(chars, function(c) if (utf8ToInt(c) > 127) 2 else 1, numeric(1)))
  }
  
  logo_vws <- vapply(logo_core, get_vw, numeric(1))
  max_vw <- max(logo_vws)
  
  border <- strrep("=", max_vw)
  
  pkg_v <- utils::packageVersion("leo.sc")
  meta <- c(
    paste0(">>> Welcome to leo.sc package (v", pkg_v, ") <<<"),
    paste0(">>> System Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " <<<"),
    ">>> Have fun with single cell data <<<"
  )

  # Alignment logic: ensure everything fills max_vw physically
  full_content <- c(
    border,
    vapply(seq_along(logo_core), function(i) {
      paste0(logo_core[i], strrep(" ", max_vw - logo_vws[i]))
    }, character(1)),
    border,
    vapply(meta, function(t) {
      tw <- get_vw(t); pl <- floor((max_vw - tw) / 2); pr <- max_vw - tw - pl
      paste0(strrep(" ", max(0,pl)), t, strrep(" ", max(0,pr)))
    }, character(1))
  )

  v_halo <- 2; h_halo <- 10
  
  if (nzchar(Sys.getenv("RSTUDIO")) || .has_color()) {
     packageStartupMessage(.render_vignette_physical(full_content, v_halo, h_halo, max_vw))
  } else {
    packageStartupMessage(paste(full_content, collapse = "\n"))
  }
}

.has_color <- function() {
  if (Sys.getenv("COLORTERM") == "truecolor") return(TRUE)
  if (Sys.getenv("TERM") == "dumb") return(FALSE)
  TRUE
}

# 2D Vignette Engine with Physical Grid Awareness
.render_vignette_physical <- function(lines, v_p, h_p, core_vw,
                                     fg_s = "#e82020", fg_m = "#FFFFFF", fg_e = "#2626e1",
                                     bg_c = "#000000", bg_halo = "#332212") {
  
  # Expansion buffer
  ext_lines <- c(rep(strrep(" ", core_vw), v_p), lines, rep(strrep(" ", core_vw), v_p))
  rows <- length(ext_lines)
  max_dist <- sqrt(v_p^2 + h_p^2)
  
  # Pre-calculate FG gradient for the core visual columns
  fg_ramp <- grDevices::colorRampPalette(c(fg_s, fg_m, fg_e))(core_vw)
  fg_rgb <- grDevices::col2rgb(fg_ramp)
  
  vapply(seq_along(ext_lines), function(r) {
    line <- ext_lines[r]
    chars <- strsplit(line, "")[[1]]
    
    # Render left glow
    left_aura <- vapply(1:h_p, function(c_idx) {
      dist <- sqrt(max(0, h_p + 1 - c_idx)^2 + max(0, v_p + 1 - r, r - (rows - v_p))^2)
      if (dist > max_dist) return(" ")
      w <- (dist / max_dist)^1.4
      bg <- (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_halo)
      sprintf("\033[38;2;0;0;0m\033[48;2;%d;%d;%dm ", as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]))
    }, character(1))
    
    # Render core content with physical pos mapping
    curr_vw_pos <- 1
    rendered_core <- vapply(chars, function(char) {
      cw <- if (utf8ToInt(char) > 127) 2 else 1
      dist_y <- max(0, v_p + 1 - r, r - (rows - v_p))
      bg_rgb <- if (dist_y == 0) c(0,0,0) else {
        w <- (dist_y / max_dist)^1.4
        (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_halo)
      }
      
      fg_idx <- min(core_vw, max(1, curr_vw_pos))
      curr_vw_pos <<- curr_vw_pos + cw
      
      sprintf("\033[48;2;%d;%d;%dm\033[38;2;%d;%d;%dm%s", 
              as.integer(bg_rgb[1]), as.integer(bg_rgb[2]), as.integer(bg_rgb[3]),
              fg_rgb[1,fg_idx], fg_rgb[2,fg_idx], fg_rgb[3,fg_idx], char)
    }, character(1))
    
    # Render right glow
    right_aura <- vapply(1:h_p, function(c_idx) {
      dist <- sqrt(c_idx^2 + max(0, v_p + 1 - r, r - (rows - v_p))^2)
      if (dist > max_dist) return(" ")
      w <- (dist / max_dist)^1.4
      bg <- (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_halo)
      sprintf("\033[38;2;0;0;0m\033[48;2;%d;%d;%dm ", as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]))
    }, character(1))
    
    paste0(paste0(left_aura, collapse=""), paste0(rendered_core, collapse=""), paste0(right_aura, collapse=""), "\033[0m")
  }, character(1)) |> paste(collapse = "\n")
}
