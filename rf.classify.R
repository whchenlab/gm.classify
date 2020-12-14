## -- 
## a pseudo-objective pipeline for building random forest models based on case-control data
## see rf.system.chenlab.rmd for more details ... 
## -- created on: May 30, 2020
## -- last modified: May 31, 2020
## -- last modified: June 2, 24, 25, 2020;  --
## -- last modified: Sept 25, 2020;  --

## -- data cleansing -- 
require(tidyverse);

## -- modeling, ROC and helper libraries 
require(ROCR);
require(randomForestSRC); ## use randomforestSRC instead ... 
require(caret);
require(pROC);

## -- system helper libs --
require(progress);
require(tictoc); ## show elapsed time

## -- libraries for parallel computing --;
require( doSNOW );


## -- 


### ################################################
### --- set data -- 
### --- 

## -- input: --
## -- feat.data: should be a data.frame
## --        each row is a sample, each column is a feature, 
## --        rownames are sample names
## --        the values should be relative abundances, i.e. range from 0 to 1 
## --        this table should ONLY contain relative abundances, not any meta data
## --        此表只包含相对丰度信息，不要包括 任何 meta- data -- 

## -- meta.data: should be a data.frame 
## --        each row is a sample, each column is a meta data, e.g. gender, age, bmi, group 
## --        rownames are sample names
## --        rownames should contain all sample IDs in the feat.data, it may contain extra sample names 

## -- grouping: column name in metadata in which grouping information is contained, i.e. which samples are control, and which are cases 
## -- control: to indicate which value in the grouping column indicates control 

## -- output: ---
## -- a list contains the following values --
##    "feat.raw"      "meta.raw"      "grouping.column"
##    "feat.data"     "control.group" "case.group"   
rf.setdata <- function( feat.data, meta.data, grouping = "Group", control = "CTR" ){

  ## -- elapsed time --
  tic( "Set data " );
  
  ## -- create a list to contain all information ...
  rf <- list();
  
  ## ========================================================================
  ## -- sanity check --
  ## -- 1. if meta.data contains info for all samples;
  sample_names.all <- rownames( feat.data );
  sample_names.not_in_meta_data <- sample_names.all[ !sample_names.all %in% rownames( meta.data ) ];
  if( length( sample_names.not_in_meta_data ) > 0 ){
    stop( "San check failed: the following samples in the feature data are not in the metadata: ",  paste( sample_names.not_in_meta_data, collapse = ", " ), "." );
  }
  
  ## -- 2. if relative abundances --
  feat.range <- range( feat.data, na.rm = TRUE );
  if( feat.range[1] < 0 | feat.range[2] > 1 ){
    warning( "San check warning: feature data may not be relative abundances (ranging from 0 to 1): ", "min = ", feat.range[1], ", max = ", feat.range[2], "\n" );
  }
  
  ## -- 3. if grouping column exists --
  if( ! grouping %in% colnames( meta.data ) ){
    stop( "San check failed: meta data does not contain column with name of '", grouping, "'." );
  }
  
  ## -- 4. if the grouping column contain not two values --
  groups <- table( meta.data[, grouping] );
  control.count <- groups[ control ];
  if( length( table( meta.data[, grouping] ) ) != 2 ){
    stop( "San check failed: the grouping column in meta data should contain EXACTLY TWO categories." );
  }
  
  if( is.na( control.count ) ){
    stop( "San check failed: the control group '", control , "' does not seem to exist in your meta data." );
  }
  
  ## ======================================================================== 
  ## -- set data --
  ## -- 1. back up all raw data --
  rf[["feat.raw"]] <- feat.data;
  rf[["meta.raw"]] <- meta.data;
  
  ## -- 2. add grouping data to feat.data --
  dat <- as.data.frame(feat.data);
  dat[, "Group"] <- as.character( meta.data[ rownames( feat.data ), grouping ] );
  
  ## -- 3. get case and controls 
  group.names <- names(groups);
  case <- group.names[ group.names != control ];
  
  rf[["feat.data"]] <- dat;
  rf[["control.group"]] <- control;
  rf[["case.group"]] <- case;
  rf[["grouping.column"]] <- grouping;
  
  ## ======================================================================== 
  ## -- wrap up and return rf --
  cat( "All look great.", "\n");
  
  cat( "-----------------------------------------------------------", "\n");
  cat( "  In total, the feature data contains ", nrow(feat.data), " samples, with ", ncol(feat.data), " features.\n", sep = "" );
  cat( "  Case group: '", case, "', w/ ", groups[case], " samples, \n", sep = "");
  cat( "  Control group: '", control, "', w/ ", groups[control], " samples, \n", sep = "");
  cat( "-----------------------------------------------------------", "\n\n");
  
  cat( " Data size: ", format( object.size(rf), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc(); ## show elapsed time 
  ## -- return --
  return( rf );
}

### ################################################
### --- filter features and samples based on preset roles -- 
### --- input: 
### --- rf contains: feat.raw, meta.raw, feat.data, grouping.column, control.group, case.group 
### --- feat.except : a vector contains features that do not participate the filtering 
### --- In each dataset, species with max abundance <0.1% in all samples and species whose average abundance across all samples below 0.01% were removed. 
### --- Pathways with maximum relative abundances below 1×10-6 in all samples and with zero value in over 15% samples were removed.

### --- input:
##      - rf: a rf object, output of rf.setdata --
##      - datatype: taxon or pathway, default is taxon --
##      - max.percentage: relative abundances are presented as percentage (0-100,max.percentage = 100) or not (0-1,max.percentage = 1), default is 100
##      - feat.except: a vector contains features that should not be excluded --
rf.filter <- function( rf, feat.except = NULL, datatype = "taxon", max.percentage = 100 ){
    
    ## -- elapsed time --
    tic( "Filter data " );
    
    cat( "Parameters:", "\n");
    cat( " feat.except  :", feat.except, "\n");
    cat( " datatype     :", datatype, "\n");
    cat( " percentage   : 0 -", max.percentage, "\n\n");
    
    ## ========================================================================
    ## -- sanity check --
    ## -- 1. if feat.data is available --
    if( ! "feat.data" %in% names(rf) ){
        stop( "San check failed: input list 'rf' does not contain 'feat.data', please check your input.");
    }
    
    ## ========================================================================
    ## -- get data --
    feat.data <- rf[["feat.data"]];
    
    ## -- temporarily remove certain columns from the data --
    feat.not.to.be.filtered <- intersect( c( feat.except, "Group" ), colnames( feat.data ) ); ## 不参与 filter的 feature 列表；
    feat.tmp <- feat.data[, !colnames( feat.data ) %in% feat.not.to.be.filtered ];
    
    if( datatype == "taxon" ){
        ## species with max abundance <0.1% in all samples and species whose average abundance across all samples below 0.01% were removed. 
        abun.cutoff.max <- ifelse( max.percentage==100, 0.1, 0.001 );    ## max cutoff是 0.1% 
        abun.cutoff.mean <- ifelse( max.percentage==100, 0.01, 0.0001 ); ## mean cutoff是 0.01% 
        
        ## -- 
        col.max <- apply( feat.tmp, 2, max );
        col.mean <- colMeans( feat.tmp );
        
        ## -- do filtering --
        feat.tmp.filtered <- feat.tmp[,  col.max > abun.cutoff.max & col.mean > abun.cutoff.mean ];
    } else {
        ## -- Pathways with maximum relative abundances below 1×10-6 in all samples and with zero value in over 15% samples were removed.
        abun.cutoff.max <- ifelse( max.percentage==100, 1e-4, 1e-6 ); ## abundance cutoff是 1e-6 
        
        samples.count <- nrow( feat.tmp );
        
        ## -- do filtering --
        feat.tmp.filtered <- feat.tmp[ , colSums( feat.tmp > abun.cutoff.max ) > samples.count * 0.15 ];
    }
    
    ## -- discarded features --
    discard.feat.count <- ncol( feat.tmp ) - ncol( feat.tmp.filtered );
    
    ## -- put features back --
    feat.tmp.filtered2 <- cbind( feat.tmp.filtered, feat.data[ feat.not.to.be.filtered] ); ## -- great --
    
    ## -- print --
    cat( "-----------------------------------------------------------", "\n");
    cat( " # of total features     :", ncol( feat.data )-1, "\n");
    cat( " # of features removed   :", discard.feat.count, "\n");
    cat( " # of features retained  :", ncol( feat.tmp.filtered2 )-1, "\n");
    cat( "-----------------------------------------------------------", "\n");
    
    ## 
  
    cat( " Data size: ", format( object.size(rf), units = "MB" ), "\n", sep = "");
    
    ## ========================================================================
    ## get data to return --
    rf[["feat.data.before.filtering"]] <- feat.data;
    rf[["feat.data"]] <- feat.tmp.filtered2;
    rf[["feat.data.filtered"]] <- TRUE;
    
    ## -- return ... 
    toc(); ## show elapsed time 
    ## -- return --
    return( rf );
}

### ################################################
### --- filter features and samples based on preset roles -- 
### --- input: 
### --- rf: a rf object, output of rf.setdata --
### --- a: Significance threshold for wilcox test--
##3 --- LDA.score: Significance threshold for lda--
rf.featureselect.diff <- function(rf, a = 0.05, LDA.score = 2){
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  require(MASS);#lda
  require(dplyr);
  require(reshape2);
  require(ggpubr);
  require(ggplot2);
  rFromWilcox<-function(wilcoxModel,N)
  {
    z<-qnorm(wilcoxModel$p.value/2)
    r<-z/sqrt(N)
    return(r)
  }
  ## -- elapsed time --
  tic( "Backwards Feature Selection" );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  
  if( !"feat.data" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'feat.data', please check if you have trained the models." );
  }
  if( !"control.group" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'control.group', please check if you have trained the models." );
  }
  if( !"case.group" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'case.group', please check if you have trained the models." );
  }
  ## -- 
  cat("Get raw feature data ... \n", sep = "");
  
  ## -- data 
  feat.data <- rf[["feat.data"]];
  feat.data$Group <- as.factor(feat.data$Group);
  features <- feat.data %>% subset(select=(-Group)) %>% names;
  control <- rf$control.group;
  case <- rf$case.group
  wilcox.result <- data.frame(features = features,
                              p.value = NA,
                              greater= NA,
                              less = NA,
                              lda = 0,
                              effect_size = 0)
  rownames(wilcox.result) <- features;
  for (i in features) {
    #H0:case=control H1:case!=control
    w <- wilcox.test(feat.data[which(feat.data$Group == case),i],feat.data[which(feat.data$Group == control),i]);
    wilcox.result[i, "p.value"] <- w$p.value;
    wilcox.result[i, "effect_size"] <- abs(rFromWilcox(w,nrow(feat.data)))
    #H0:case=control H1:case>control
    w <- wilcox.test(feat.data[which(feat.data$Group == case),i],feat.data[which(feat.data$Group == control),i], alternative = "greater");
    wilcox.result[i, "greater"] <- w$p.value;
    #H0:case=control H1:case<control
    w <- wilcox.test(feat.data[which(feat.data$Group == case),i],feat.data[which(feat.data$Group == control),i], alternative = "less");
    wilcox.result[i, "less"] <- w$p.value;
    #LDA
    try({lda <- lda(Group ~ feat.data[,i], data = feat.data)
    ldmean<-(mean(lda$scaling[,"LD1"]*feat.data[which(feat.data$Group==case),i])-mean(lda$scaling[,"LD1"]*feat.data[which(feat.data$Group==control),i]))
    realmean<-mean(feat.data[which(feat.data$Group==case),i])-mean(feat.data[which(feat.data$Group==control),i])
    lda2real<-(ldmean+realmean)/2
    wilcox.result[i, "lda"] <- lda2real},silent = TRUE)
  }
  wilcox.result$fdr.pvalue <- p.adjust(wilcox.result$p.value, method = "fdr");
  wilcox.result$fdr.greater <- p.adjust(wilcox.result$greater, method = "fdr");
  wilcox.result$fdr.less <- p.adjust(wilcox.result$less, method = "fdr");
  wilcox.feature <- wilcox.result[which(wilcox.result$fdr.pvalue < a),"features"]%>%as.vector();
  wilcox.featdata <- feat.data[, c(wilcox.feature,"Group")];
  wilcox.boxplot <- melt(wilcox.featdata, value.name = "rel.abun") %>%
    mutate(log.rela = log10(rel.abun + 0.00001)) %>%
    mutate(Group2 = ifelse(Group==case,"case","control")) %>% 
    mutate(Group2 = factor(Group2, levels = c("control", "case"))) %>%
    ggplot(aes(x=Group2,y=log.rela)) +
    geom_boxplot(aes(fill=Group2)) +
    facet_wrap(~ variable, scales="free") +
    labs(y = "relative abundance(log10)",x="") +
    scale_fill_manual(values=c("#1F6D41", "#DA001B"))+
    stat_compare_means(method = "wilcox.test");
  
  wilcox.result$significant <- NULL;
  wilcox.result[which((wilcox.result$fdr.pvalue < a) & (abs(wilcox.result$lda) > LDA.score)), "significant"] <- "significant";
  wilcox.result$trend <- NULL;
  wilcox.result[which(wilcox.result$fdr.greater < a & abs(wilcox.result$lda) > LDA.score), "trend"] <- "case";
  wilcox.result[which(wilcox.result$fdr.less < a & abs(wilcox.result$lda) > LDA.score), "trend"] <- "control";
  sign.feature <- wilcox.result[which(wilcox.result$significant == "significant"),"features"]%>%as.vector();
  sign.featdata <- feat.data[, c(sign.feature,"Group")];
  sign.boxplot <- melt(sign.featdata, value.name = "rel.abun") %>%
    mutate(log.rela = log10(rel.abun + 0.00001)) %>%
    mutate(Group2 = ifelse(Group==case,"case","control")) %>% 
    mutate(Group2 = factor(Group2, levels = c("control", "case"))) %>%
    ggplot(aes(x=Group2,y=log.rela)) +
    geom_boxplot(aes(fill=Group2)) +
    facet_wrap(~ variable, scales="free") +
    labs(y = "relative abundance(log10)",x="") +
    scale_fill_manual(values=c("#1F6D41", "#DA001B"))+
    stat_compare_means(method = "wilcox.test");
  #+geom_jitter(aes(x=Group,y=rel.abun))##散点
  #p+stat_compare_means( label = "p.signif")
  #sign.boxplot;
  
  sign.result <- wilcox.result[which(wilcox.result$significant=="significant"),c("features","lda","trend")];
  sign.lda.plot <- 
    ggplot(sign.result,aes(x=factor(reorder(features,lda)),y=lda,fill=trend))+
    geom_bar(stat = "identity",width = 0.78, position = position_dodge(0.7))+
    labs(x="",y = "LDA SCORE")+
    coord_flip()+#翻转
    theme_bw()+#删除灰背景
    theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+#删除纵向网格
    theme(axis.text.y = element_blank())+#删除刻度标签
    theme(axis.ticks = element_blank())+#删除刻度标签
    theme(panel.border = element_blank())+#去掉边框
    geom_text(aes(y=ifelse(lda>0,-1,1),label = features))+
    scale_fill_manual(values=c("case" = "#DA001B", "control" = "#1F6D41"));
  
  
  ## -- print some statistics on screen --
  cat("  # ",length(sign.feature)," features selected by wilcox.test (FDR-P.value < ",a," and LDA score > ",LDA.score," )\n", sep = "");
  cat("  NOTE: a new list will be created with necessary data to run 'rf.train()' function. \n");
  rf.new <- list( "feat.raw" = rf[["feat.raw"]],
                  "meta.raw" = rf[["meta.raw"]],
                  "control.group" = rf[["control.group"]],
                  "case.group" = rf[["case.group"]],
                  "grouping.column" = rf[["grouping.column"]],
                  "feat.data" = sign.featdata,
                  "sign.feature" = sign.feature,
                  "wilcox.lda.result" = wilcox.result,
                  "sign.lda.plot" = sign.lda.plot,
                  "sign.boxplot" = sign.boxplot,
                  "wilcox.boxplot" = wilcox.boxplot);
  
  cat( "  Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc();
  return( rf.new );
}

### ################################################
### --- filter features and samples based on preset roles -- 
### --- input: 
### --- rf: a rf object, output of rf.setdata --
### --- number: either the number of folds or number of resampling iterations
### --- repeats: for repeated k-fold cross-validation: the number of complete sets of folds to compute#指定抽样组数
### --- sizes: a numeric vector of integers corresponding to the number of features that should be retained
### --- eg, number = 10, repeats = 10: do 10 folds 10 repeats
rf.featureselect.rfe <- function(rf, number = 10, repeats = 10, sizes = c(1:20) ){
  
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  require(caret);
  
  ## -- elapsed time --
  tic( "Backwards Feature Selection" );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  
  if( !"feat.data" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'feat.data', please check if you have trained the models." );
  }
  if( !"meta.raw" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'meta.raw', please check if you have trained the models." );
  }
  if( !"grouping.column" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'grouping.column', please check if you have trained the models." );
  }
  
  ## -- 
  cat("Get raw feature data and metadata... \n", sep = "");
  
  ## -- data 
  feat.data <- rf[["feat.data"]];
  feat.tmp <- feat.data %>% dplyr::select(-"Group");
  meta.tmp <- rf[["meta.raw"]] %>% dplyr::select(rf[["grouping.column"]]);
  meta.data <- meta.tmp[,1];
  ## ========================================================================
  ## -- generates a control object and do a simple backwards selection ...
  ctrl <- rfeControl(functions = rfFuncs, method = "repeatedcv", number = number, repeats = repeats, verbose = FALSE, returnResamp = "final");
  num <- ncol(feat.tmp);
  if(num>20 && num%%5!=0){
    size <- c(1:20,5*(5:(num%/%5)),num);
  }else if(num>20 && num%%5==0){
    size <- c(1:20,5*(5:(num%/%5)));
  }else{
    size<-c(1:num);}
  if(!is.null( sizes )){size=sizes};
  cat("  Parameter used to select features by backwards(caret::rfe): \n  --------functions = rfFuncs, method = repeatedcv, number = ",number, ", repeats = ",repeats,", sizes :","\n", sep = "");
  print(size);
  Profile <- rfe(feat.tmp, meta.data,  rfeControl = ctrl, sizes = size);
  ref.featdata = feat.data[, c( Profile$optVariables, "Group" )]
  ## -- print some statistics on screen --
  cat("  # ",Profile$optsize," features selected \n", sep = "");
  cat("  NOTE: a new list will be created with necessary data to run 'rf.train()' function. \n");
  rf.new <- list( "feat.raw" = rf[["feat.raw"]],
                  "meta.raw" = rf[["meta.raw"]],
                  "control.group" = rf[["control.group"]],
                  "case.group" = rf[["case.group"]],
                  "grouping.column" = rf[["grouping.column"]],
                  "feat.data" = ref.featdata,
                  "rfe.Profile" = Profile);
  
  cat( "  Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc();
  return( rf.new );
}
### ################################################
### --- modeling ...  -- 
## -- June 24, 2020; added parameter: classwt
rf.train <- function( rf, fold = 10, resample = 10, parallel = T, num.cores = 10, class.weights = NULL ){
  
  ## -- requires ... 
  require(purrr); 
  
  ## -- elapsed time --
  tic( "Model training " );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  if( sum( c("feat.data", "control.group", "case.group") %in% names(rf)  ) < 3  ){
    stop( "San check failed: input data does not contain all required slots: feat.data, control.group, case.group, please try from the begining." );
  }
  
  
  ## ========================================================================
  ## -- get data from the input list --
  dat <- rf[["feat.data"]]; 
  dat$Group <- as.factor( dat$Group );
  
  meta.data <- rf[["meta.raw"]];
  grouping <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
  control.group <- rf[["control.group"]];
  case.group <- rf[["case.group"]];
  
  ## -- 
  cat("Training models on ", nrow( dat ), " samples, with ", ncol( dat ) - 1, " features.\n", sep = "" );
  cat("Case group: ", case.group, ", control group: ", control.group, "\n", sep = "" );
  cat("  - running ", resample, " resample and ", fold , " fold modeling ... \n", sep = "");
  
  ## -- 
  ## -- class weights --
  if( isTRUE( class.weights ) ){ 
      cat("  - using class weight option of randomForest classifier ... \n" );
  } else {
      cat("  - NOT using class weight option of randomForest classifier ... \n" );
  }
  
  ## ========================================================================
  ## -- modeling --
  
  ##require group name is "Group"
  ## -- create confusion matrix --
  sample.split <- createMultiFolds( as.vector( dat$Group ), times = resample, k = fold  );
  
  ## -- create a few lists to save final data 
  rf.models <- list();
  pred.matrix <- list();
  confusion.result <- list();
  
  ## --- ----
  if( parallel ){
    ## ========================================================================
    ##    parallel version 
    ## ========================================================================
    ## 1. decide number of cpu cores to use 
    available.cores <- parallel::detectCores() - 1;
    if( num.cores > available.cores  ){
      num.cores = available.cores;
    }
    
    if( num.cores < 1 ){
      num.cores = 1;
    }
    
    cat( "  - running parallel modeling using ", num.cores, " cpu cores ... \n", sep = "" );
    
    ## -- 2. modeling using designated CPU cores  --
    cl <- makeSOCKcluster(num.cores);
    registerDoSNOW(cl);
    
    pb <- txtProgressBar(min=1, max=length( sample.split ), style=3);
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    model.results <- 
      foreach( split = sample.split, .packages="randomForestSRC", .options.snow=opts ) %dopar% {
        ## -- get data ready --
        train.data <- dat[ split , ];
        test.data <- dat[ -split , ];
        
        ## -- train, and get model prediction results ... 
        if(  isTRUE( class.weights ) ){
            this.model <- rfsrc( Group ~ .,train.data, proximity = T, importance = T, ntree = 500,
                                 case.wt = randomForestSRC:::make.wt( train.data$Group ),
                                 sampsize = randomForestSRC:::make.size( train.data$Group ));
        } else {
            this.model <- rfsrc(Group ~ .,train.data, proximity = T, importance = T, ntree = 500);
        }
        ## 
        pred <- predict( this.model, test.data, type = "prob");
        this.pred.matrix <- pred$predicted ;
        rownames( this.pred.matrix ) <- rownames( test.data );
        
        ## 
        this.confusion.result <- predict( this.model, test.data, type= "response");
        
        return( list( "model" = this.model,
                      "pred.matrix" = this.pred.matrix,
                      "confusion.result" = this.confusion.result) );
      }  ## end of foreach ... 
    close(pb);
    stopCluster(cl);
    
    ## -- 3. post process to get required data ... 
    iters <- names( sample.split );
    for( idx in 1:length(model.results) ){
      iter <- iters[ idx ];
      result.list <- model.results[[idx]];
      
      ## -- get data from list and save them to other lists ... 
      rf.models[[iter]] <- result.list[["model"]];
      pred.matrix[[iter]] <- result.list[["pred.matrix"]];
      confusion.result[[iter]] <- result.list[["confusion.result"]];
    }
    
  } else {
    ## ========================================================================
    ##    regular version 
    ## ========================================================================
    ## -- make a progress bar for modeling -- 
    pb1 <- progress_bar$new( format = "Modeling [:bar] :percent in :elapsed", clear = FALSE, total =  fold * resample );
    
    ## -- do modeling ... 
    for (iter in names(sample.split)) {
      ## -- get data ready --
      train.data <- dat[sample.split[[iter]], ]
      test.data <- dat[-sample.split[[iter]], ]
      
      ## -- train, and predict ... 
      if(  isTRUE( class.weights ) ){
          this.model <- rfsrc( Group ~ .,train.data, proximity = T, importance = T, ntree = 500,
                               case.wt = randomForestSRC:::make.wt(train.data$Group),
                               sampsize = randomForestSRC:::make.size(train.data$Group));
      } else {
          this.model <- rfsrc(Group ~ .,train.data, proximity = T, importance = T, ntree = 500);
      }
      
      ## -- predict --
      pred <- predict( this.model, test.data, type = "prob");
      this.pred.matrix <- pred$predicted;
      rownames( this.pred.matrix ) <- rownames( test.data );
      
      ## -- confusion matrix --
      this.confusion.result <- predict(this.model, test.data, type= "response")
      
      ## -- save results ...
      rf.models[[iter]] = this.model;
      pred.matrix[[iter]] <- this.pred.matrix;
      confusion.result[[iter]] <- this.confusion.result;
      
      ## -- tick --
      pb1$tick();
    }
  }
  
  ## --- importance scores --- ##
  importance <- lapply(rf.models, function(data){ imp <-  data$importance[, c("all")] * 100 }); ## -- ... 
  importance.df <- purrr::reduce(importance, cbind)
  importance.df <- data.frame(importance.df, stringsAsFactors = F);
  names(importance.df) <- names(rf.models);
  
  
  ## -- assemble rf --
  ## -- modeling results ...
  rf[["sample.split"]] <- sample.split;
  rf[["rf.models"]] <- rf.models;
  rf[["pred.matrix"]] <- pred.matrix;
  rf[["confusion.result"]] <- confusion.result;
  rf[["importance.df"]] <- importance.df;
  
  cat( " Data size: ", format( object.size(rf), units = "MB" ), "\n", sep = "");
  
  ## -- return ... 
  toc(); ## show elapsed time 
  ## -- return --
  return( rf );
}

### ################################################
### --- evaluate ; created on May 30, 2020 --
### --- last modified : June 1, 2020 ; June 24, 2020 --
rf.evaluate <- function( rf ){
  
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  require(tidyverse);
  
  ## -- modeling, ROC and helper libraries 
  require(ROCR);
  require(pROC);
  
  ## -- elapsed time --
  tic( "Model performance evaluation " );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  if( sum( "pred.matrix" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: pred.matrix, please check if you have trained the models." );
  }
  
  if( sum( "meta.raw" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: meta.raw, please check if you have the meta data correctly assigned." );
  }
  
  if( sum( "grouping.column" %in% names(rf) ) < 1  ){
    stop( "San check failed: input data does not contain required slot: grouping.column, please check if you have the meta data correctly assigned." );
  }
  
  ## -- 
  cat("Evaluating model performance ... \n", sep = "");
  
  ## ========================================================================
  ## -- get data from previous results --
  pred.matrix       <- rf[["pred.matrix"]];
  meta.data         <- rf[["meta.raw"]];
  grouping          <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
  control.group     <- rf[["control.group"]];
  case.group        <- rf[["case.group"]];
  
  ## ========================================================================
  ## -- get model performance ... 
  ## -- 1. iterate the prediction matrix, combine all predictions, and calculate mean scores for each sample --
  res.pred.df <- pred.matrix %>% 
    map_df( function(x) { tibble( Sample = rownames(x), Control = x[, control.group], Case = x[, case.group ] ) } ) %>% 
    group_by( Sample ) %>% 
    summarize( Control = mean(Control), Case = mean(Case) ) %>% as.data.frame( . );
  
  colnames( res.pred.df ) <- c("Sample", control.group, case.group );
  
  ## -- 2. summarize prediction results and add the original grouping data ... 
  res.pred.df$Predicted <- apply(res.pred.df, 1, function(data){ names(which.max(data[ c( control.group, case.group ) ]))} );
  res.pred.df$Group <-  meta.data[ res.pred.df$Sample, grouping ];
  
  ## -- 3. calculate AUC & CI --
  rf.new <- evaluate.internal.use( res.pred.df, control.group, case.group );
  
  toc(  ); ## show elapsed time ... 
  return( rf.new );
}

### ################################################
### --- feature selection -- 
rf.featureselect <- function(rf, top = NA, importance.cutoff = 1 ){
  
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  
  ## -- elapsed time --
  tic( "Feature selection " );
  
  ## ========================================================================
  ## -- sanity check --
  ## 1. check if input contains necessary data 
  if( !"importance.df" %in% names(rf)  ){
    stop( "San check failed: input data does not contain required slot: 'importance.df', please check if you have trained the models." );
  }
  
  if( !"feat.data" %in% names(rf) ){
    stop( "San check failed: input data does not contain required slot: 'feat.data', please check if you have trained the models." );
  }
  
  ## -- 
  cat("Get top features ... \n", sep = "");
  
  ## -- data 
  importance.df <- rf[["importance.df"]];
  feat.data <- rf[["feat.data"]];
  
  ## ========================================================================
  ## -- arrange median relative importance [0-100]; relative importance values can be negative ... 
  feat.relaimpo <- importance.df %>% 
    rowMeans(.) %>% 
    tibble( relaimpo = ., species = names(.) ) %>% 
    arrange( desc (relaimpo ) );
  
  if( is.na(top) & max( feat.relaimpo$relaimpo ) < importance.cutoff ){
      cat("------------------------------------  ERROR -----------------------------------------\n");
      cat("There is no feature with importance scores higher than your input cutoff:", importance.cutoff, ", please revise accordingly.\n" );
      cat(" # of features:", ncol( feat.data ) - 1, "\n" );
      feat.impo.range <- range( feat.relaimpo$relaimpo );
      cat(" feature importance range", feat.impo.range[1], "to", feat.impo.range[2], "\n" );
      cat(" feature importance distribution:\n");
      print( summary( feat.relaimpo$relaimpo  ) );
      
      stop();
  }
  
  ## ========================================================================
  ## -- 1. select top features
  ## -- if top is set, choose the number of top features as specified by 'top'
  ## -- otherwise, choose top features according to the 'importance.cutoff' 
  if( !is.na( top ) ){
    cat("The top parameter is set, get top ", top, " features ... \n");
    feat.select <- feat.relaimpo %>% dplyr::slice( 1:top );
    param_used <- "top";
  } else {
    feat.select <- feat.relaimpo %>% filter( relaimpo > importance.cutoff );
    param_used <- "importance.cutoff";
  }
  
  ## -- 2. filter data by selected features --
  feat.data.select <- feat.data[, c( feat.select$species, "Group" )];
  
  ## -- 3. relative importance plot --
  feat.importance.select <- as.data.frame( t( importance.df[ feat.select$species, ] ) ) %>% gather( feature, relaimpo );
  
  feat.importance.select.plot <- 
    feat.importance.select %>% ggplot( aes( reorder( feature, relaimpo, median) ,  relaimpo) ) + 
    geom_boxplot( aes( fill = feature ) ) + coord_flip() + theme( legend.position = "none" );
  
  ## ========================================================================
  ## -- print some statistics on screen --
  ## -- 
  cat("  Parameter used to select top features: ", param_used, ", \n", sep = "");
  cat("  # of top features selected: ", nrow( feat.select ) , ", \n");
  cat("  Max importance: ", max( feat.select$relaimpo ) , ", \n");
  cat("  Min importance: ", min( feat.select$relaimpo ) , ", \n\n");
  
  cat("NOTE: a new list will be created with necessary data to run 'rf.train()' function. \n");
  
  rf.new <- list( "feat.raw" = rf[["feat.raw"]],
                  "meta.raw" = rf[["meta.raw"]],
                  "control.group" = rf[["control.group"]],
                  "case.group" = rf[["case.group"]],
                  "grouping.column" = rf[["grouping.column"]],
                  "feat.data" = feat.data.select,
                  "feat.select" = feat.select,
                  "feat.select.importance.plot" = feat.importance.select.plot,
                  "feat.rank.median" = feat.relaimpo );
  
  cat( " Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc();
  return( rf.new );
}

### ################################################
### --- predict ---
##  --- input, :
##      -- a rf.train (a rf.train result) object 
##      -- a data.frame contains relative abundances
##         each row is a sample, each column is a feature (taxon)
##         it should NOT contain grouping information 
rf.predict <- function(rf,  test.data.frame ){
    
    ## -- elapsed time --
    tic( "Prediction using external data " );
    
    ## ========================================================================
    ## -- sanity check --
    ## 1. check if input contains necessary data 
    if( sum( "rf.models" %in% names(rf) ) < 1  ){
        stop( "San check failed: input data does not contain required slot: 'rf.models', please check if you have trained the models." );
    }
    
    ## -- 2. input data should be a data frame ; 输入数据必需是 data.frame 
    if( !is.data.frame( test.data.frame )  ){
        stop( "San check failed: input data should be a data.frame. Please double check!!." );
    }
    
    ## 3. check if contain non-numeric columns -- 输入数据不能包含 非数值 列；
    if( sum( !unlist(lapply( test.data.frame , is.numeric))) > 0  ){
        stop( "San check failed: input data.frame contains non-numeric column(s), please remove them or re-prepare your data." );
    }
    
    ## ---select features needed by models ----------## 
    ### 从模型中找到建模所需的 feature  (taxon names ) ... 
    feat.names <- rf$rf.models[[1]]$xvar.names; 
    
    ## 注意，rf.data$feat.data里除了feature外，还有Grouping列
    useless.features<-setdiff( colnames( test.data.frame ),  feat.names );
    lacking.features<-setdiff( feat.names, colnames( test.data.frame ));
    
    ## 移除无用的feature
    test.feat.data <- test.data.frame[ !colnames(test.data.frame) %in% useless.features ];
    
    ## 插入缺失的feature 填入0
    for(f in lacking.features){
        test.feat.data[f] = vector(mode = "numeric", nrow(test.feat.data) );
    }
    
    ## warning messages ...
    if( length(lacking.features) > 0 ){
        cat("\n--------------------------------------- WARNING ---------------------------------------\n");
        cat("# of required features: ", length( feat.names ), "\n", sep = "");
        cat("# of features NOT found in your input: ", length( lacking.features ), "\n", sep = "" );
        cat("  the predicted results may NOT be reliable due to these missing features!!\n");
        cat("---------------------------------------------------------------------------------------\n\n");
    }
    
    ## -- reorder columns according to those in the models ... 
    test.feat.data <- test.feat.data[, feat.names ]; 
    
    ## 预测时，feature的顺序重要不 ???? 
    
    ## --
    pb1 <- progress_bar$new( format = "Predicting [:bar] :percent in :elapsed", clear = FALSE, total =  length( rf$rf.models ) );
    
    ## ---for external validation -----------##
    pred.result <- list();
    
    for (iter in names(rf$rf.models)) {
        pred <- predict(rf$rf.models[[iter]], test.feat.data, type = "prob");
        pred.df <- data.frame(pred$predicted, stringsAsFactors = F);
        pred.df$Sample <- rownames(test.feat.data);
        pred.result[[iter]] <- pred.df;
        
        ## -- tick --
        pb1$tick();
    }
    
    pred.result.df <- pred.result %>%
        purrr::reduce(rbind) %>% group_by(Sample);
    
    ## -- 注意 gather 与 spread 的用法 ... 
    pred.result.df <- pred.result.df %>% 
        gather( Verdict, Prob, -Sample ) %>% 
        group_by( Sample, Verdict ) %>% 
        summarise( Prob = mean( Prob ) ) %>% 
        spread( Verdict, Prob );
    
    pred.result.df <- as.data.frame( pred.result.df );
    pred.result.df$Predicted <- apply(pred.result.df, 1, function(data){ names(which.max(data[ c( rf$control.group, rf$case.group ) ]))} );
    
    toc(  ); ## show elapsed time ... 
    
    return( list(
        "userdata" = test.data.frame,
        "userdata.reformatted" = test.feat.data, 
        "predicted.results.df" = pred.result.df
    ) );
}


### ################################################
### created on June 24, 2020 --
### -- input:
###   - rf: a rf.train object, a list that contains :
###         'rf.models', control.group, case.group, grouping.column 
###   - ext.feat.data
###   - ext.meta.data, feat and meta data for external validation; the same as the two input files for tf.setdata --
### -- output: a list containing the following objects: 
###   - 
rf.external.validation <- function( rf, ext.feat.data, ext.meta.data ){
    
    ## -- elapsed time --
    tic( "External validation " );
    
    ## ========================================================================
    ## -- sanity check & run predictions--
    ## 1. check if input contains necessary data 
    if( sum( c("rf.models", "control.group", "case.group", "grouping.column" ) %in% names( rf )  ) < 4  ){
        stop( "San check failed: input data does not contain all required slots: 'rf.models', 'control.group', 'case.group', 'grouping.column', please check your input." );
    }
    
    ## -- get data from tf ... 
    grouping <- rf[["grouping.column"]]; ## the column in meta.data for the grouping information, i.e. CASE, CONTROL
    control.group <- rf[["control.group"]];
    case.group <- rf[["case.group"]];
    
    ## -- 2. try to validate feat and meta data by build a new rf object using rf.setdata --
    rf.ext.data <- rf.setdata( ext.feat.data, ext.meta.data, grouping = grouping, control = control.group );
    
    ## -- 3. get feat data with and without groups, and use rf.predict to do predictions --
    feat.w.grp <- rf.ext.data$feat.data;
    feat.w.o.grp <- feat.w.grp[ , !colnames(feat.w.grp) %in% "Group" ];
    
    ext.data.pred <- rf.predict( rf, feat.w.o.grp );
    
    ## -- 4. evaluate performance --
    res.pred.df <- ext.data.pred$predicted.results.df;
    res.pred.df$Group <-  ext.meta.data[ res.pred.df$Sample, grouping ];
    
    #### ======= THE FOLLOWING CODES ARE THE SAME FROM rf.evaluate =======
    rf.new <- evaluate.internal.use( res.pred.df, control.group, case.group );
    
    toc(  ); ## show elapsed time ... 
    return( rf.new );
}

### ################################################
### --- plot multiple ROC and compare -- 
## -- input, a list of rf.evaluation results ... 
rf.plot.multiple.roc <- function( ..., labels = NULL, colours = NULL ){ 
  
  ## -- requires ... 
  require(tictoc); ## show elapsed time
  require(tidyverse);
  
  ## -- system helper libs --
  require(progress);
  
  ## -- elapsed time --
  tic( "Merge multiple AUC curves and plot " );
  
  ## ========================================================================
  ## sanity check ... 
  rocs <- list( ... );
  roc_count <- length(rocs);
  
  if( !is.null( labels ) & length(labels) < roc_count  ){
    stop( "San check failed: 'labels' is not empty but its length is fewer than the number of rf.eval objects, please check your data." );
  }
  
  roc.df.all <- tibble( fpr = numeric(), tpr = numeric(), label = character() );
  for( idx in 1:length(rocs) ){
    roc <- rocs[[idx]];
    label <- if_else( is.null( labels ), paste0( "Dataset ", idx )  , labels[idx] );
    data.roc.df <- rocs[[idx]][["roc.df"]] %>% dplyr::select( "fpr", "tpr" ) %>% mutate( label = paste0( label, ": ", roc[["auc.string"]] ) );
    roc.df.all <- bind_rows( roc.df.all, data.roc.df );
  }
  
  roc.plot <- 
    ggplot( roc.df.all, aes(x = fpr, y = tpr, colour = label) ) + geom_line() + theme_bw() + 
    labs(x = "False Positve Rate", y = "True Positive Rate", title = "ROC plot")  + 
    theme( plot.title = element_text(hjust = .5, face = "bold", size = 11), 
           axis.text = element_text(size = 9), 
           legend.justification = c(0.99,0.01), legend.position = c(0.99,0.01), legend.title = element_blank() );
  
  if( !is.null( colours ) ){
    roc.plot <- roc.plot + scale_color_manual( values = colours );
  }
  
  print( roc.plot );
  
  return( roc.plot );
}


### ################################################
### a function used only internally --
evaluate.internal.use <- function( res.pred.df, control.group, case.group ){
    
    forestpred <- ROCR::prediction( res.pred.df[, case.group ], res.pred.df[, "Group" ],
                                    label.ordering = c( control.group, case.group ));
    
    forestperf <- ROCR::performance(forestpred, "tpr", "fpr");
    data.roc.df <- data.frame( fpr = forestperf@x.values[[1]],  tpr = forestperf@y.values[[1]] );
    
    ## -- 4. use pROC to calculate AUC with CI -
    auc <- pROC::roc( res.pred.df[, "Group"], res.pred.df[, case.group ], levels = c( control.group, case.group ), ci = TRUE, plot = FALSE );
   
    auc.ci.values <- round( as.numeric( auc$ci ), digits = 2 ); ## a vector of three
    auc.string <- paste0( auc.ci.values[2], " (95%CI: ",  auc.ci.values[1], "-",  auc.ci.values[3], ")")
    
    ## -- 5. plot --
    data.roc.plot <- 
        ggplot(data.roc.df, aes(x = fpr, y = tpr) ) + geom_line() + theme_bw() + 
        labs(x = "False Positve Rate", y = "True Positive Rate", title = "ROC plot")  + 
        annotate( "text", x = .75, y = .75, label = paste0( "AUC:", auc.string ), size = 4) +
        theme( plot.title = element_text(hjust = .5, face = "bold", size = 11), 
               axis.text = element_text(size = 9) );
    
    print( data.roc.plot );
    
    ## ========================================================================
    ### calculate confusion matrix --
    confusion <- table( res.pred.df[,  c( "Group", "Predicted" ) ] );
    errors <- round( 1 - diag( confusion ) / rowSums( confusion ), digits = 3 );
    errors.overall <- round(  1 - ( sum(diag(confusion)) / sum( confusion ) ) , digits = 3 );
    
    ## ========================================================================
    ## -- print on-screen statistics ...  
    cat("\n--------------------------------------------\n");
    cat( "Overall performance (AUC): ", auc.string, "\n\n" );
    cat( "Predicted vs. Truth (Group): \n" );
    cat( "Control group: ", control.group, "\n" );
    cat( "Case group: ", case.group, "\n\n" );
    print( confusion );
    cat( "---------\n" );
    cat( "Error rates: \n" );
    print( errors );
    cat("\nOverall error rate: ", errors.overall, "\n", sep = "");
    cat("--------------------------------------------\n\n");
    ## ========================================================================
    ## -- 
    
    rf.new <- list( "roc.df" = data.roc.df, 
                    "roc.plot" = data.roc.plot, 
                    "auc" = auc.ci.values[2], 
                    "auc.with.ci" = auc.ci.values, 
                    "auc.obj" = auc, 
                    "auc.string" = auc.string,
                    "pred.df" = res.pred.df, 
                    "error.rates" = errors,
                    "overall.error.rate" = errors.overall);
    
    cat( "Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
    cat("---\n");
    
    return(rf.new);
}

rf.master <- function(feat.data, meta.data, grouping = "Group", control = "CTR",datatype = "taxon", max.percentage = 100, feat.except = NULL, fold = 10, resample = 10, parallel = T, num.cores = 10, class.weights = TRUE, sizes =NULL ){
  ## -- requires ... 
  require(tictoc); 
  
  ## -- elapsed time --
  tic( "Modeling " );
  ## -- Prepare data --
  cat( "...... 01 ......\n...... do rf.setdata() ......\n");
  rowdata.01 <- rf.setdata(feat.data, meta.data, grouping = grouping, control = control);
  cat( "\n...... 02 ......\n...... do rf.filter() ......\n");
  filter.data.02 <- rf.filter( rowdata.01, feat.except = feat.except, datatype = datatype, max.percentage = max.percentage );
  #cat( "\n...... 03 ......\n...... do rf.featureselect.wilcox() ......\n");
  #wilcox.data.03 <- rf.featureselect.wilcox(filter.data.02, a=0.05);
  #num <- ncol(wilcox.data.03[["feat.data"]])-1;
  num <- ncol(filter.data.02[["feat.data"]])-1;
  if(num>20 && num%%5!=0){
    size <- c(1:20,5*(5:(num%/%5)),num);
  }else if(num>20 && num%%5==0){
    size <- c(1:20,5*(5:(num%/%5)));
    }else{
      size<-c(1:num);}
  if(!is.null( sizes )){size=sizes};
  cat( "\n...... 04 ......\n...... do rf.featureselect.rfe() ......\n");
  ref.data.04 <- rf.featureselect.rfe(filter.data.02, number = fold, repeats = resample, sizes = size );
  num <- ref.data.04[["rfe.Profile"]]$optsize;
  #ref.data.04 <- rf.featureselect.rfe(wilcox.data.03, number = fold, repeats = resample, sizes = size );
  
  ## -- modeling --
  cat( "\n...... 05 ......\n...... do rf.train() ......\n");
  rf <- rf.train( ref.data.04, fold = fold, resample = resample, parallel = parallel, num.cores = num.cores, class.weights = class.weights );
  cat( "\n...... 06 ......\n...... do rf.evaluate() ......\n");
  rf.eval <- rf.evaluate(rf);
  cat( "\n...... 07 ......\n...... do rf.featureselect() ......\n");
  feat.importance <- rf.featureselect(rf, top = num);
  feat.importance.top20 <- rf.featureselect(rf, top = 20);
  rf.new <- list( "rf" = rf,
                  #"feat.wilcox" = wilcox.data.03[["wilcox.result"]],
                  "feat.rfe" = ref.data.04[["rfe.Profile"]],
                  "feat.select" = ref.data.04[["rfe.Profile$optVariables"]],
                  "feat.select.importance.plot" = feat.importance[["feat.select.importance.plot"]],
                  "feat.select.importance.plot.top20" = feat.importance.top20[["feat.select.importance.plot"]],
                  "roc.plot" = rf.eval[["roc.plot"]],
                  "rf.evaluate" = rf.eval);
  
  cat( " Data size: ", format( object.size(rf.new), units = "MB" ), "\n", sep = "");
  
  ## -- 
  toc();
  return( rf.new );
}

