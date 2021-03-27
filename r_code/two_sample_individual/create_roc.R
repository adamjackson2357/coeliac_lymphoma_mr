## Functions to create roc curves

get_fpr_fnr <- function(case_covars, model, Y, label){
  
  pred <- predict(model, type="response")
  rocr <- prediction(pred, case_covars[,Y], label.ordering = NULL)
  print(performance(rocr, measure = "auc")@y.values)
  roc <- performance(rocr, measure="tpr", x.measure="fpr")
  values = data.frame(fpr=roc@x.values[[1]],
                      fnr=roc@y.values[[1]])
  values$Variable <- label
  return(values)
}

create_roc_curve <- function(fpr_fnr){
  prs_roc <- fpr_fnr %>%
    ggplot(aes(x = fpr, y = fnr, colour=Variable)) +
    geom_line() +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") +
    theme_minimal() +
    xlim(0, 1) + ylim(0, 1) +
    xlab("False Positive Rate") + ylab("False Negative Rate")
  return(prs_roc)
}