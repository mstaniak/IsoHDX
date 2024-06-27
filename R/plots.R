#' @export
#' @importFrom ggplot2 ggplot geom_bar facet_grid theme_bw theme scale_fill_discrete
plotSpectraComparison = function(fitted_iso_model) {
  comp_long = data.table::melt(fitted_iso_model[["FinalComparison"]],
                               measure.vars = c("Intensity", "ExpectedPeak"),
                               variable.factor = FALSE,
                               variable.name = "Type", value.name = "Intensity")
  comp_long[, Type := ifelse(Type == "Intensity", "observed", "predicted")]
  ggplot(comp_long,
         aes(x = IntDiff, y = Intensity, fill = Type)) +  
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(Time~paste(Peptide, Charge)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_discrete(name = "intensity")
}

#' @export
#' @importFrom ggplot2 geom_linerange geom_point facet_wrap geom_line geom_linerange
plotPredictedProbabilities = function(fitted_iso_model, selected_times = NULL,
                                      type = "by_exchanged") {
  if (is.null(selected_times)) {
    selected_times = unique(fitted_iso_model[["FittedProbabilities"]][["Time"]])
  }
  if (type == "by_exchanged") {
    ggplot(fitted_iso_model[["FittedProbabilities"]][Time %in% selected_times],
           aes(x = reorder(as.character(NumExchanged), NumExchanged), y = Probability)) +
      geom_point() +
      geom_linerange(aes(ymin = Lower, ymax = Upper)) +
      facet_grid(Segment~Time, scales = "free_x") +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    ggplot(fitted_iso_model[["FittedProbabilities"]][Time %in% selected_times],
           aes(x = Time, y = Probability, 
               group = reorder(as.character(NumExchanged), NumExchanged), 
               color = reorder(as.character(NumExchanged), NumExchanged))) +
      geom_point() +
      geom_line() +
      geom_linerange(aes(ymin = Lower, ymax = Upper)) +
      facet_grid(~Segment, scales = "free_x") +
      theme_bw() +
      theme(legend.position = "bottom")
  }
}
