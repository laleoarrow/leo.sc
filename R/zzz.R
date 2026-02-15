
.onAttach <- function(libname, pkgname) {
  # [BASE TRANSLATION] Faithful 100% to User 1085 manual design
  # Using escaped Unicode for source compatibility with R-CMD-check
  # Legend: █=\u2588, ║=\u2551, ╗=\u2557, ═=\u2550, ╔=\u2554, ╝=\u255d, ╚=\u255a
  logo_lines <- c(
    "  \u2588\u2588\u2557     \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2557          \u2588\u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2588\u2588\u2588\u2588\u2557  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2588\u2588\u2557         \u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255D  ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551   \u2588\u2588\u2551   \u2588\u2588\u2557   \u255A\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2551       ",
    "  \u2588\u2588\u2551     \u2588\u2588\u2554\u2550\u2550\u255D  \u2588\u2588\u2551   \u2588\u2588\u2551   \u255A\u2550\u255D    \u255A\u2550\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2551       ",
    "  \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D         \u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D\u255A\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255D  ",
    "  \u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D\u255A\u2550\u2550\u2550\u2550\u2550\u2550\u255D \u255A\u2550\u2550\u2550\u2550\u2550\u255D          \u255A\u2550\u2550\u2550\u2550\u2550\u255D  \u255A\u2550\u2550\u2550\u2550\u2550\u255D  "
  )
  
  # Calculate VISUAL width of logo (considering wide characters)
  logo_widths <- nchar(logo_lines, type = "width")
  max_visual_w <- max(logo_widths)
  
  # Dynamic Border Creation (matching physical width)
  border_line <- strrep("=", max_visual_w)
  
  pkg_v <- utils::packageVersion("leo.sc")
  meta <- c(
    paste0(">>> Welcome to leo.sc package (v", pkg_v, ") <<<"),
    paste0(">>> System Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " <<<"),
    ">>> Have fun with single cell data <<<"
  )

  # Assemble the target content box
  # Pad all lines to match physical max_visual_w
  full_content <- c(
    border_line,
    vapply(seq_along(logo_lines), function(i) {
      paste0(logo_lines[i], strrep(" ", max_visual_w - logo_widths[i]))
    }, character(1)),
    border_line,
    vapply(meta, function(t) {
      tw <- nchar(t, type = "width")
      pl <- floor((max_visual_w - tw) / 2)
      pr <- max_visual_w - tw - pl
      paste0(strrep(" ", max(0,pl)), t, strrep(" ", max(0,pr)))
    }, character(1))
  )

  # Canvas expansion for Glow (8 horizontal cells, 2 vertical cells)
  h_glow <- 10; v_glow <- 2
  
  if (nzchar(Sys.getenv("RSTUDIO")) || .has_color()) {
     packageStartupMessage(.render_halo_vignette(full_content, v_glow, h_glow, max_visual_w))
  } else {
    # Plain text fallback with spaces
    padded_full <- c(rep("", v_glow), full_content, rep("", v_glow))
    packageStartupMessage(paste(padded_full, collapse = "\n"))
  }
}

.has_color <- function() {
  if (Sys.getenv("COLORTERM") == "truecolor") return(TRUE)
  if (Sys.getenv("TERM") == "dumb") return(FALSE)
  TRUE
}

# Precision Vignette Engine - Visual-Width Aware
.render_halo_vignette <- function(lines, v_p, h_p, core_w, 
                                 fg_s = "#e82020", fg_m = "#FFFFFF", fg_e = "#2626e1",
                                 bg_c = "#000000", bg_h = "#332212") {
  
  rows_full <- length(lines) + 2 * v_p
  cols_full <- core_w + 2 * h_p
  max_dist <- sqrt(v_p^2 + h_p^2)
  
  # Foreground gradient mapped to physical core width
  fg_ramp <- grDevices::colorRampPalette(c(fg_s, fg_m, fg_e))(core_w)
  fg_rgb <- grDevices::col2rgb(fg_ramp)
  
  # Create full rendering buffer (lines included)
  all_lines <- c(rep(strrep(" ", core_w), v_p), lines, rep(strrep(" ", core_w), v_p))
  
  vapply(seq_along(all_lines), function(r) {
    line <- all_lines[r]
    chars <- strsplit(line, "")[[1]]
    char_widths <- nchar(chars, type = "width")
    
    # Pre-pad line with h_p spaces
    res <- character(length(chars) + 2) # Start with spacer
    
    # Left padding aura
    left_aura <- vapply(1:h_p, function(c_idx) {
      dist <- sqrt(max(0, h_p + 1 - c_idx)^2 + max(0, v_p + 1 - r, r - (rows_full - v_p))^2)
      if (dist > max_dist) return(" ")
      w <- (dist / max_dist)^1.4
      bg <- (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_h)
      sprintf("\033[48;2;%d;%d;%dm ", as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]))
    }, character(1))
    
    # Core content
    curr_v_col <- 1 
    rendered_core <- vapply(seq_along(chars), function(i) {
      char <- chars[i]; cw <- char_widths[i]
      dist_y <- max(0, v_p + 1 - r, r - (rows_full - v_p))^2
      # Core area (dist_x is always 0 inside core_w)
      bg_rgb <- c(0,0,0)
      
      # Map fg color to char position
      fg_idx <- min(core_w, max(1, curr_v_col))
      curr_v_col <<- curr_v_col + cw
      
      sprintf("\033[48;2;0;0;0m\033[38;2;%d;%d;%dm%s", fg_rgb[1,fg_idx], fg_rgb[2,fg_idx], fg_rgb[3,fg_idx], char)
    }, character(1))
    
    # Right padding aura
    right_aura <- vapply(1:h_p, function(c_idx) {
      dist <- sqrt(c_idx^2 + max(0, v_p + 1 - r, r - (rows_full - v_p))^2)
      if (dist > max_dist) return(" ")
      w <- (dist / max_dist)^1.4
      bg <- (1-w) * grDevices::col2rgb(bg_c) + w * grDevices::col2rgb(bg_h)
      sprintf("\033[48;2;%d;%d;%dm ", as.integer(bg[1]), as.integer(bg[2]), as.integer(bg[3]))
    }, character(1))
    
    paste0(paste0(left_aura, collapse=""), paste0(rendered_core, collapse=""), paste0(right_aura, collapse=""), "\033[0m")
  }, character(1)) |> paste(collapse = "\n")
}
