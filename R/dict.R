#' Canonical gene markers for T-cell subtypes (Zhang 2018)
#' @keywords internal
.tcell_genes <- list(
  # CD8 T cells
  CD8_Tn          = c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1"),
  CD8_Tcm         = c("CCR7","SELL","IL7R","CD27","CD28","PRF1","GZMA","CCL5",
                      "GPR183","S1PR1"),
  CD8_Temra_Teff  = c("KLRG1","CX3CR1","FCGR3A","FGFBP2","PRF1","GZMH","TBX21",
                      "EOMES","S1PR1","S1PR5"),
  CD8_Tem         = c("GZMK","CXCR4","CXCR3","CD44"),
  CD8_Trm         = c("CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2",
                      "NR4A3","CD69","ITGAE"),
  CD8_IEL         = c("CD160","KIR2DL4","TMIGD2","KLRC1","KLRC2","KLRC3",
                      "NR4A1","NR4A2","NR4A3","IKZF2","ENTPD1","CD69","ITGAE"),
  CD8_TEX         = c("HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","GZMB",
                      "MIR155HG","TNFRSF9","ITGAE"),
  MAIT            = c("SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA"),
  # CD4 T cells
  CD4_Tn          = c("CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1"),
  CD4_Blood_Tcm   = c("CCR7","SELL","PTGER2","ICAM2","ANXA1","ANXA2","S1PR1"),
  CD4_TemRA_TEFF  = c("KLRG1","CX3CR1","NKG7","PRF1","GNLY","GZMH","TBX21",
                      "CTSW","S1PR1","S1PR5"),
  CD4_Normal_Tcm  = c("CCR7","TCF7","RGS1","CD69"),
  CD4_Trm         = c("CD69","KLRB1","PTGER4","IL7R","CXCR6","NR4A1","NR4A2",
                      "NR4A3","MYADM"),
  TFH             = c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST","MAGEH1",
                      "BTLA","ICOS","PDCD1","CD200"),
  CD4_Tem_Th1like = c("GZMK","GZMA","CCL5","IFNG","RUNX3","EOMES","CXCR3",
                      "CXCR4","CD44"),
  Th17            = c("IL23R","RORC","IL17A","FURIN","CTSH","CCR6","KLRB1",
                      "CAPG","ITGAE"),
  Th1_like        = c("CXCL13","IFNG","CXCR3","BHLHE40","GZMB","PDCD1","HAVCR2",
                      "ICOS","IGFLR1","ITGAE"),
  Blood_Treg      = c("FOXP3","IL2RA","IL10RA","IKZF2","RTKN2","CDC25B","S1PR4"),
  TFR             = c("FOXP3","IL2RA","CXCR5","PDCD1","IL10","CCR4","CD69")
)

.t_nk_marker <- list(
  `T` = c("CD3D","CD3E"), # CD3+ Complex
  CD4 = c("CD4","CD40LG"), # CD4+; Helper T cells
  CD8 = c("CD8A","CD8B"), # CD8+; CD8A~α; CD8B~β
  Tn = c("CCR7","SELL","LRRN3"), # CCR7-Naive
  Tem = c("GZMK","CCL5"),
  Tcm = c("CD27","CD28"),
  Temra = c("FGFBP2", "KLRG1"), # CD8+CD45RA-CCR7-
  Treg = c("FOXP3"), # CD4+CD25+FOXP3+; Regulatory T cells
  CTL = c("GZMB","PRF1"), # Cytotoxic T cells
  # DPT = c("TYROBP"), # Not so sure
  DNT = c("TRDV2"), # `CD4-8-` TRDV2 can also be gdT's marker
  gdT =c("TRGV9"), # Integrating single-cell RNA and T cell/B cell receptor sequencing with mass cytometry reveals dynamic trajectories of human peripheral immune cells from birth to old age
  MAIT = c("KLRB1"), # Also known as `T-mito`;KLRB1~CD161
  # NKT = c("S1PR1", "IL32"), # Not so sure，perhaps just check NK markers
  # NK markers below
  CD16_NK = c("FCGR3A","FCER1G","KLRF1"),
  CD56_NK = c("NCAM1","XCL1","CMC1")
)

# Marker hub --------
#' Marker hub
#'
#' A minimal nested list for quick access to all kinds of markers collected in `leo.sc`.
#'
#' @format List with two elements: \code{$data}, \code{$source}.
#' @export
#' @examples
#' devtools::load_all("/Users/leoarrow/project/mypackage/leo.sc")
#' t_nk_marker <- leo.marker$t_nk_marker$data
leo.marker <- list(
  t_nk_marker = list(
    data = .t_nk_marker,
    name = "Leo curated T/NK cell markers",
    note = "Some can refer to Azimuth (https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC)"
  ),
  tcell_genes = list(
    data = .tcell_genes,
    name = "Zhang2018",
    title = "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer",
    source = "https://www.nature.com/articles/s41586-018-0694-x#Fig9",
    doi = "https://doi.org/10.1038/s41586-018-0694-x",
    note = "Canonical maybe a little bit old?"
  ),
  leo.note = list(
    name = "sc note for Leo.sc marker hub. Come check this when you have no clues.",
    msg1 = "IL7R can be a good marker expressed in Tn, Tcm and Tem. But it is mostly in Tcm in VKH2024 case as it is essential to T cell survival.",
    note = "You can use `leo.marker$t_nk_marker$data` to access the T/NK cell markers."
  )
)
