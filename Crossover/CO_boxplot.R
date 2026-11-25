library(ggplot2)

script_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    if (identical(sys.function(i), base::source)) {
      source_frame <- sys.frame(i)
      path <- source_frame$ofile
      if (is.character(path) && length(path) == 1L) {
        if (!is.null(source_frame$owd)) path <- basename(path)
        return(dirname(normalizePath(path, mustWork = TRUE)))
      }
    }
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    stop("Open the script and use Source, or run it with Rscript.", call. = FALSE)
  }
  path <- sub("^--file=", "", file_arg[[1L]])
  path <- gsub("~+~", " ", path, fixed = TRUE)  # Decode spaces in --file.
  dirname(normalizePath(path, mustWork = TRUE))
})
input_file <- file.path(script_dir, "CO-results.xlsx")
# Experimental gamete totals include samples with zero crossovers.
expected_n <- c(WT_female = 275L, WT_male = 225L, ebs_female = 277L, ebs_male = 244L)
group_order <- names(expected_n)
group_colors <- c("#0089BB", "#00B51B", "#FF0034", "#FF8000")

co_raw <- as.data.frame(readxl::read_excel(input_file, sheet = 1))
required_columns <- c("group", "sample_id", "co_num")
if (!all(required_columns %in% names(co_raw))) {
  stop("The workbook must contain group, sample_id and co_num.")
}
if (anyNA(co_raw[required_columns]) ||
    any(!co_raw$group %in% group_order) ||
    any(!nzchar(as.character(co_raw$sample_id)))) {
  stop("The workbook contains missing values, an unknown group or an empty sample ID.")
}
if (!is.numeric(co_raw$co_num) || any(!is.finite(co_raw$co_num)) ||
    any(co_raw$co_num < 0 | co_raw$co_num != floor(co_raw$co_num))) {
  stop("co_num must contain non-negative integers.")
}

co_per_gamete <- aggregate(co_num ~ group + sample_id, co_raw, sum)
names(co_per_gamete)[names(co_per_gamete) == "co_num"] <- "CO_count"

observed_n <- table(factor(co_per_gamete$group, levels = group_order))
missing_n <- expected_n - as.integer(observed_n)
if (any(missing_n < 0)) {
  stop("Sample count exceeds the experimental total for: ",
       paste(group_order[missing_n < 0], collapse = ", "))
}
# Use placeholder IDs for zero-CO gametes missing from the table.
zero_groups <- rep(group_order, missing_n)
if (length(zero_groups)) {
  zero_rows <- data.frame(
    group = zero_groups,
    sample_id = paste0("__zero_placeholder__", zero_groups, "_", sequence(missing_n)),
    CO_count = 0L,
    stringsAsFactors = FALSE
  )
  if (any(zero_rows$sample_id %in% co_per_gamete$sample_id)) {
    stop("A zero-CO placeholder ID matches an existing sample ID.")
  }
  co_per_gamete <- rbind(co_per_gamete, zero_rows)
}
cat("Zero-CO samples added by group:\n")
print(missing_n)
co_per_gamete$group <- factor(co_per_gamete$group, levels = group_order)
co_per_gamete <- co_per_gamete[
  order(co_per_gamete$group, co_per_gamete$sample_id),
]
row.names(co_per_gamete) <- NULL
stopifnot(sum(co_per_gamete$CO_count) == sum(co_raw$co_num))

co_summary <- do.call(rbind, lapply(seq_along(group_order), function(i) {
  values <- co_per_gamete$CO_count[co_per_gamete$group == group_order[i]]
  if (!length(values)) stop("No gametes found for group: ", group_order[i])

  q <- quantile(values, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
  iqr <- q[3] - q[1]
  data.frame(
    group = group_order[i], x = i, n = length(values),
    total_CO = sum(values), mean_CO = mean(values),
    q1 = q[1], median = q[2], q3 = q[3],
    low = min(values[values >= q[1] - 1.5 * iqr]),
    high = max(values[values <= q[3] + 1.5 * iqr]),
    group_color = group_colors[i]
  )
}))
print(co_summary[c("group", "n", "total_CO", "mean_CO")], row.names = FALSE)

set.seed(20261013)
co_per_gamete$x <- as.integer(co_per_gamete$group) +
  runif(nrow(co_per_gamete), min = -0.29, max = 0.29)
co_per_gamete$point_color <- group_colors[as.integer(co_per_gamete$group)]

axis_labels <- as.expression(list(
  bquote(atop("WT female", .(sprintf("(n=%d)", co_summary$n[1])))),
  bquote(atop("WT male", .(sprintf("(n=%d)", co_summary$n[2])))),
  bquote(atop(italic(ebs) ~ "female", .(sprintf("(n=%d)", co_summary$n[3])))),
  bquote(atop(italic(ebs) ~ "male", .(sprintf("(n=%d)", co_summary$n[4]))))
))

co_plot <- ggplot() +
  geom_rect(
    data = co_summary,
    aes(xmin = x - 0.22, xmax = x + 0.22, ymin = q1, ymax = q3,
        fill = group_color, color = group_color),
    alpha = 0.16, linewidth = 0.6
  ) +
  geom_segment(data = co_summary, aes(x = x, xend = x, y = low, yend = q1),
               color = "#444444", linewidth = 0.6) +
  geom_segment(data = co_summary, aes(x = x, xend = x, y = q3, yend = high),
               color = "#444444", linewidth = 0.6) +
  geom_segment(data = co_summary, aes(x = x - 0.11, xend = x + 0.11, y = low, yend = low),
               color = "#444444", linewidth = 0.6) +
  geom_segment(data = co_summary, aes(x = x - 0.11, xend = x + 0.11, y = high, yend = high),
               color = "#444444", linewidth = 0.6) +
  geom_point(data = co_per_gamete, aes(x = x, y = CO_count, color = point_color),
             size = 2, shape = 16, alpha = 0.6) +
  geom_segment(
    data = co_summary,
    aes(x = x - 0.22, xend = x + 0.22, y = median, yend = median, linetype = "Median"),
    color = "#333333", linewidth = 0.9
  ) +
  geom_point(data = co_summary, aes(x = x, y = mean_CO, shape = "Mean"),
             size = 4.5, stroke = 0.8, fill = "black", color = "white") +
  geom_text(
    data = co_summary,
    aes(x = x, y = 10.3, label = sprintf("%.2f", mean_CO), color = group_color),
    inherit.aes = FALSE, size = 9
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_manual(name = NULL, values = c(Mean = 23)) +
  scale_linetype_manual(name = NULL, values = c(Median = "solid")) +
  scale_x_continuous(breaks = 1:4, labels = axis_labels, expand = c(0, 0)) +
  scale_y_continuous(breaks = 0:11, minor_breaks = NULL,
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-0.4, 11)) +
  labs(x = NULL, y = "Number of crossovers per gamete",
       title = "A  Individual gametes and summary") +
  guides(shape = guide_legend(order = 1, override.aes = list(color = "black", size = 4)),
         linetype = guide_legend(order = 2)) +
  theme_classic(base_size = 20, base_family = "sans") +
  theme(
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 20, margin = margin(t = 10)),
    axis.text.y = element_text(size = 17, margin = margin(r = 7)),
    axis.title.y = element_text(size = 24, margin = margin(r = 12)),
    plot.title = element_text(size = 24, hjust = 0, margin = margin(b = 35)),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 15),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.text = element_text(size = 17),
    legend.key.width = grid::unit(1.2, "cm"),
    legend.margin = margin(t = 25),
    legend.box.spacing = grid::unit(10, "pt")
  )

print(co_plot)
