
.onAttach <- function(libname, pkgname) {
  banner <- paste(
    "  ██╗     ███████╗ ██████╗          ██████╗  ██████╗ ",
    "  ██║     ██╔════╝██╔═══██╗         ██╔════╝██╔════╝ ",
    "  ██║     █████╗  ██║   ██║   ██╗   ╚█████╗ ██║      ",
    "  ██║     ██╔══╝  ██║   ██║   ╚═╝    ╚═══██╗██║      ",
    "  ███████╗███████╗╚██████╔╝         ██████╔╝╚██████╗ ",
    "  ╚══════╝╚══════╝ ╚═════╝          ╚═════╝  ╚═════╝ ",
    sep = "\n"
  )
  
  # Prepare metadata lines (User format)
  pkg_version <- packageVersion("leo.sc")
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  meta_text <- c(
    paste0(">>> Welcome to leo.sc package (v", pkg_version, ") <<<"),
    paste0(">>> System Time: ", current_time, " <<<"),
    ">>> Have fun with single cell data <<<"
  ) 
  
  # Calculate width for centering based on the banner
  banner_lines <- strsplit(banner, "\n")[[1]]
  banner_width <- max(nchar(banner_lines))
  
  # Center metadata lines
  meta_lines <- vapply(meta_text, .center_line, width = banner_width, FUN.VALUE = character(1))
  
  # Combine
  full_msg <- paste(c(banner, meta_lines), collapse = "\n")
  
  # Display
  if (nzchar(Sys.getenv("RSTUDIO")) || .has_color()) {
     packageStartupMessage(.gradient(full_msg))
  } else {
    packageStartupMessage(full_msg)
  }
}

# Helper functions
.center_line <- function(text, width) {
  if (width <= 0) return(text)
  pad <- floor((width - nchar(text)) / 2)
  if (pad < 0) pad <- 0
  paste0(strrep(" ", pad), text)
}

.has_color <- function() {
  if (Sys.getenv("COLORTERM") == "truecolor") return(TRUE)
  if (Sys.getenv("TERM") == "dumb") return(FALSE)
  TRUE
}

.gradient <- function(text, start = c(160, 0, 0), mid = c(255, 255, 255), end = c(0, 0, 160)) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  colored <- vapply(lines, function(line) {
    if (!nzchar(line)) return("")
    chars <- strsplit(line, "")[[1]]
    if (length(chars) == 0) return("")
    
    n <- length(chars)
    # Map t from 0..1
    t <- if (n == 1) 0 else seq(0, 1, length.out = n)
    
    # Interpolate
    rgb_vals <- vapply(t, function(val) {
      if (val < 0.5) {
        # Interpolate start -> mid
        # Normalize val to 0..1 range for this segment: val * 2
        p <- val * 2
        start + (mid - start) * p
      } else {
        # Interpolate mid -> end
        # Normalize val to 0..1 range for this segment: (val - 0.5) * 2
        p <- (val - 0.5) * 2
        mid + (end - mid) * p
      }
    }, numeric(3))
    
    # rgb_vals is 3 x n
    r <- round(rgb_vals[1, ])
    g <- round(rgb_vals[2, ])
    b <- round(rgb_vals[3, ])
    
    paste0(sprintf("\033[38;2;%d;%d;%dm%s", r, g, b, chars), collapse = "")
  }, character(1))
  paste(paste0(colored, "\033[0m"), collapse = "\n")
}

