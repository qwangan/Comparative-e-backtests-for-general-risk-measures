write_numerical_tables <- function() {
  dir.create("tables",showWarnings=FALSE)
  rates <- read.csv("results/traditional_iid_rejections.csv",check.names=FALSE)
  scenarios <- c("Baseline","-5% VaR","-10% VaR","-5% ES","-10% ES","-5% Both")
  rows <- vapply(c(2,5,10,20),function(threshold) {
    values <- rates$rejection[match(paste(scenarios,threshold),paste(rates$scenario,rates$threshold))]
    paste0(threshold," & ",paste(sprintf("%.4f",values),collapse=" & ")," \\\\")
  },character(1))
  writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}","\\toprule",
    "Threshold & Baseline & -5\\% VaR & -10\\% VaR & -5\\% ES & -10\\% ES & -5\\% Both \\\\",
    "\\midrule",rows,"\\bottomrule","\\end{tabular}"),"tables/traditional_iid.tex")
  rates <- read.csv("results/comparative_iid_rejections.csv",check.names=FALSE)
  models <- c("-50%","-20%","Noisy","+20%","+50%","Exact")
  labels <- c("-50\\%","-20\\%","truth","+20\\%","+50\\%","opt")
  make_rows <- function(thresholds) vapply(seq_along(models),function(j) {
    entries <- vapply(seq_along(models),function(i) {
      if(i==j) return("--")
      cells <- rates[rates$internal==models[i] & rates$standard==models[j],]
      values <- cells$red[match(thresholds,cells$threshold)]*100
      paste(formatC(values,format="f",digits=1),collapse="/")
    },character(1))
    paste0(labels[j]," & ",paste(entries,collapse=" & ")," \\\\")
  },character(1))
  for(version in c("legacy","05")) {
    thresholds <- if(version=="legacy") c(2,5,10) else 40
    writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}","\\toprule",
      paste0("Standard / Internal & ",paste(labels,collapse=" & ")," \\\\") ,
      "\\midrule",make_rows(thresholds),"\\bottomrule","\\end{tabular}"),
      paste0("tables/comparative_iid_",version,".tex"))
  }
  invisible(TRUE)
}
if(sys.nframe()==0L) write_numerical_tables()
