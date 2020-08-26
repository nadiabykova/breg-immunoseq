wilcoxon_pair_test = function(cols,df,test_name,g1,g2){
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  res = data.frame()
  for (i in c(1:nrow(cols_data))){
    i_line = cols_data %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    wt.pair = tryCatch(
      wilcox.test(df.i[,g1],df.i[,g2],paired = TRUE,conf.int=T),
      error = function(e) return(NA))
    if (is.na(wt.pair[1])){
      pv = NA
      est = NA
    } else {
      pv = wt.pair$p.value
      est = wt.pair$estimate
    }
    res = rbind(res,cbind(i_line,data.frame(pvalue=pv,estimate=est)))
  }
  res$test_type="Wilcoxon paired"
  res$test_name=test_name
  res$test_group1=g1
  res$test_group2=g2
  return(res)
}

wilcoxon_general_test = function(cols,df,test.column,g1,g2){
  variable = c(test.column)
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  res = data.frame()
  for (i in c(1:nrow(cols_data))){
    i_line = cols_data %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    df.i1= df.i %>% filter(get(variable) == g1)
    df.i2 = df.i %>% filter(get(variable)== g2)
    wt = tryCatch(
      wilcox.test(df.i1$value,df.i2$value,conf.int=T),
      error = function(e) return(NA))
    if (is.na(wt[1])){
      pv = NA
      est = NA
    } else {
      pv = wt$p.value
      est = wt$estimate
    }
    res = rbind(res,cbind(i_line,data.frame(pvalue=pv,estimate=est)))
  }
  res$test_type="Wilcoxon"
  res$test_name=test.column
  res$test_group1=g1
  res$test_group2=g2
  return(res)
}

wilcox_general_comb_test = function(cols,df,test.column){
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  groups = as.data.frame(t(combn(unique(df[,test.column]), 2)))
  colnames(groups)=c("test_group1","test_group2")
  cols_data_group = merge(cols_data,groups)
  res = data.frame()
  for (i in c(1:nrow(cols_data_group))){
    i_line = cols_data_group %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    df.i.g1 = df.i %>% filter(get(test.column)==unique(df.i$test_group1))
    df.i.g2 = df.i %>% filter(get(test.column)==unique(df.i$test_group2))
    wt = tryCatch(
      wilcox.test(df.i.g1$value,df.i.g2$value,conf.int=T),
      error = function(e) return(NA))
    if (is.na(wt[1])){
      pv = NA
      est = NA
    } else {
      pv = wt$p.value
      est = wt$estimate
    }
    res = rbind(res,cbind(i_line,data.frame(pvalue=pv,estimate=est)))
  }
  res$test_type="Wilcoxon"
  res$test_name=test.column
  return(res)
}

t_test_general_test = function(cols,df,test.column,g1,g2){
  variable = c(test.column)
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  res = data.frame()
  for (i in c(1:nrow(cols_data))){
    i_line = cols_data %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    df.i1 = df.i %>% filter(get(variable)==g1)
    df.i2 = df.i %>% filter(get(variable)==g2)
    t = tryCatch({
      t.test(df.i1$value,df.i2$value)
    }, error = function(error_condition) {NA})
    if (is.na(t[3])) {
      est=NA
      pv=NA
    } else{
      est=t$estimate[1]-t$estimate[2]
      pv=t$p.value
    }
    res = rbind(res,cbind(i_line,data.frame(pvalue=pv,estimate=est)))
  }
  res$test_type="t-test"
  res$test_name=test.column
  res$test_group1=g1
  res$test_group2=g2
  return(res)
}

anova_general_test = function(cols,df,test.column){
  variable = c(test.column)
  df = as.data.frame(df)
  cols_data = df %>% select(cols) %>% unique
  res = data.frame()
  for (i in c(1:nrow(cols_data))){
    i_line = cols_data %>% slice(i)
    df.i = merge(df,i_line,by=cols)
    aov.pval = tryCatch({
        summary(aov(value~get(variable),df.i))[[1]][["Pr(>F)"]][1]
    }, error = function(error_condition) {NA})
    res = rbind(res,cbind(i_line,data.frame(pvalue = aov.pval, estimate = NA)))
  }
  res$test_type="anova"
  res$test_name=test.column
  res$test_group1=NA
  res$test_group2=NA
  return(res)
}

get.adjusted = function(adj.cols,df,m="holm"){
  df = as.data.frame(df)
  padj_data = df %>% select(adj.cols) %>% unique
  res.adj = data.frame()
  for (i in c(1:nrow(padj_data))){
    i_line = padj_data %>% slice(i)
    df.i = merge(df,i_line,by=adj.cols)
    df.i$adjusted = p.adjust(df.i$pvalue,method=m)
    df.i$adj.method = m
    df.i$adj.cols = paste(adj.cols,collapse ="|")
    df.i$adj.group.size = nrow(df.i)
    res.adj = rbind(res.adj,df.i)
  }
  return(res.adj)
}
