## Functions to visualise the PRS

create_hist <- function(case_covars, variable, legend_title){
  return(case_covars %>%
           ggplot(aes_string(x="PRS", colour=variable)) +
           geom_density() +
           theme_classic() +
           theme(legend.position = c(0.8, 0.8)) +
           labs(color = legend_title) +
           ylab("Density"))
}

get_rocr <- function(case_covars, model, Y){
  pred <- predict(model)
  rocr <- prediction(pred, case_covars[,Y], label.ordering = NULL)
  return(rocr)
}

get_fpr_fnr <- function(rocr, label){
  roc <- performance(rocr, measure="tpr", x.measure="fpr")
  values = data.frame(fpr=roc@x.values[[1]],
                      fnr=roc@y.values[[1]])
  values$Variable <- label
  return(values)
}

plot_roc_curve <- function(fpr_fnr, legend_title){
  prs_roc <- fpr_fnr %>%
    ggplot(aes(x = fpr, y = fnr, colour=Variable)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, color="black", 
                linetype="dashed")+
    theme_classic() +
    labs(color = legend_title) +
    theme(legend.position = c(0.8, 0.2)) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("False Positive Rate") + ylab("False Negative Rate")
  return(prs_roc)
}

create_roc_curve <- function(case_covars, model, Y, label){
  rocr <- get_rocr(case_covars, model, Y)
  auc <- paste0("AUC = ", format(round(performance(rocr, measure = "auc")@y.values[[1]][1], 3), nsmall = 3))
  values <- get_fpr_fnr(rocr, label)
  roc <- plot_roc_curve(values, auc)
}


get_prec_rec <- function(rocr, label){
  prec_rec <- performance(rocr, measure="prec", x.measure="rec")
  values = data.frame(rec=prec_rec@x.values[[1]],
                      prec=prec_rec@y.values[[1]])
  values$Variable <- label
  return(values)
}

plot_prec_rec_curve <- function(prec_rec){
  prs_roc <- prec_rec %>%
    ggplot(aes(x = rec, y = prec, colour=Variable)) +
    geom_line() +
    theme_classic() +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.2)) +
    xlab("Recall") + ylab("Precision")
  return(prs_roc)
}

create_prec_rec_curve <- function(case_covars, model, Y, label){
  rocr <- get_rocr(case_covars, model, Y)
  values <- get_prec_rec(rocr, label)
  prec_rec <- plot_prec_rec_curve(values)
}
