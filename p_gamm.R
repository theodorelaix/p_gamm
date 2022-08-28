p_gamm <- function(model, series, comparison, rm.ranef = TRUE, series_length = 100, title = NULL, col = NULL, lty = NULL, width = NULL, r_range = NULL, plot_diff = FALSE, diff_pair = NULL, diff_legend = TRUE, legend_pos = "vertical", remove_level = NULL, return_data = FALSE){
  # 2022/4/28 Theodore Lai v1.2
    # add lty (linetype)
    # add width
    # add diff_pair
  # 2022/1/6 Theodore Lai v1.1
    # add control of legend position 
    # add control of difference legend
  # 2021/11/23 Theodore Lai v1
    # function created
  
  # This function is based on plot_smooth() and plot_diff() in 'itsadug' package. 
  
  # Arguments
    # model (list): A bam/gam model
    # series (character): X-axis
    # comparison (character): The comparison of the model
    # rm.ranef (logical): passed to itsadug::get_predictions()
    # series_length (numeric): The number of predicted points for the model
    # title (character): The string for the title in the plotly
    # col (vector): Customized colors for the comparison. If length of vector is not equal to the number of comparison, col will be reassigned with rainbow() with warning message.
    # lty (vector): Customized line types for the comparison.
    # r_range (2-element vector): Customized range for the radius. If not provided, the range will be calculated based on maximum radius.
    # plot_diff (logical): If FALSE (default), the plotly plot will not show the area with significant differences.
    # legend_pos (character): horizontal or vertical.
    # remove_level (vector): If provided, the comparisons matched in with elements in remove_level will be removed.
    # return_data (logical): If TRUE (default), the predicted data (in polar coordinate) will be returned for further plotting. If FALSE, the plotly plot will be returned.

  require(tidyverse)
  require(mgcv)
  require(itsadug)
  require(plotly)
  
  if(!(legend_pos %in% c("vertical", "v", "horizontal", "h"))){
    stop("undefined legend_pos")
  }
  
  # Get condition
  if(model$family$family != "gaussian"){
    condition_list <- model$xlevels[[1]] # family = "scat"
  }else{
    condition_list <- model$G$xlevels[[1]] # family = "gaussian"
  }
  

  
  if(!is.null(remove_level)){
    condition_list <- condition_list[!condition_list %in% remove_level]
  }
  
  if(is.null(col)){
    col <- rainbow(length(condition_list), start = 0, end = 0.8)
  }else{
    if(length(col) != length(condition_list)){
      print(paste0("Length of conditions (", length(condition_list), ") is not equal to length of 'col' (", length(col), "). rainbow() is used instead."))
      col <- rainbow(length(condition_list), start = 0, end = 0.8)
    }
  }
  
  if(is.null(lty)){
    lty <- rep('line', length(condition_list))
  }else{
    if(length(lty) != length(condition_list)){
      print(paste0("Length of conditions (", length(condition_list), ") is not equal to length of 'lty' (", length(lty), "). 'line' is used instead."))
      lty <- rep('line', length(condition_list))
    }
  }
  
  if(is.null(width)){
    width <- rep(1, length(condition_list))
  }else{
    if(length(width) != length(condition_list)){
      print(paste0("Length of conditions (", length(condition_list), ") is not equal to length of 'width' (", length(width), "). No adjustment is made."))
      width <- rep(1, length(condition_list))
    }
  }
  
  r_data <- tibble()
  r_diff <- tibble()
  p <- plot_ly(type="scatterpolar", mode="lines")
  for(i in 1:length(condition_list)){
    temp_list <- list(condition_list[i], seq(min(model$var.summary[[2]]), max(model$var.summary[[2]]), length = series_length))
    names(temp_list) <- c(comparison, series)
    temp <- itsadug::get_predictions(model, cond = temp_list, rm.ranef = rm.ranef, print.summary = FALSE)
    temp <- temp %>% 
      dplyr::mutate(ul = fit + CI, ll = fit-CI)
    
    r_data <- bind_rows(r_data, temp)
    
    # plotting
    p <- p %>% 
      add_trace(theta = temp$X*180/pi + 180, 
                r = temp$fit, 
                line = list(color = col[i], width = 1.5*width[i], dash = lty[i]), legendgroup = condition_list[i], name=condition_list[i]) %>% 
      add_trace(theta = temp$X*180/pi + 180, 
                r = temp$ul, 
                line = list(color = col[i], dash = "dot", width = 0.5), legendgroup = condition_list[i], showlegend=FALSE) %>% 
      add_trace(theta = temp$X*180/pi + 180, 
                r = temp$ll, 
                line = list(color = col[i], dash = "dot", width = 0.5), legendgroup = condition_list[i], showlegend=FALSE)
  }
  
  if(is.null(r_range)){
    r_range <- c(0, max(r_data$ul) + mean(r_data$ul)/10)
  }
  dtick <- round(r_range[2]/5, digits = -1)
  
  # Plot difference
  if(plot_diff){
    if(is.null(diff_pair)){
      # Make pairs
      condition_list <- as.vector(condition_list)
      col1 <- c()
      col2 <- c()
      for(i in 1:(length(condition_list)-1)){
        temp <- rep(condition_list[i], length(condition_list)-i)
        col1 <- c(col1, temp)
        temp <- rep(condition_list[(i+1):length(condition_list)], 1)
        col2 <- c(col2, temp)
      }
      comb <- tibble(Col1 = col1, Col2 = col2)
      comb <- comb %>% 
        dplyr::mutate(Col1.2 = interaction(Col1, Col2, sep = " vs. "))
    }else{
      comb <- diff_pair %>% 
        dplyr::mutate(Col1.2 = interaction(Col1, Col2, sep = " vs. "))
    }
    
    for(i in 1:nrow(comb)){
      temp_list <- list(c(comb[[i,1]], comb[[i,2]]))
      names(temp_list) <- comparison
      diff <- itsadug::plot_diff(model = model, view = series, comp = temp_list, plot = FALSE, print.summary = FALSE)
      diff <- itsadug::find_difference(mean = diff$est, se = diff$CI, diff$X)
      
      temp <- tibble(start = diff$start, end = diff$end, comparison = comb[[i,3]])
      r_diff <- bind_rows(r_diff, temp)
      
      if(is.null(diff)){
        print(paste0("There is no difference in ", comb[[i,3]]))
      }else{
        for(j in 1:length(diff$start)){
          if(j == 1){
            p <- p %>% 
              add_trace(theta = c(diff$start[j]*180/pi + 180, seq(diff$start[j]*180/pi + 180, diff$end[j]*180/pi + 180, length.out = 20), diff$end[j]*180/pi + 180), 
                        r = c(0, rep(r_range[2], 20), 0), 
                        line = list(color = "black", width=0.5), fill = "toself", fillcolor = rgb(0,0,0,max=255,alpha=25), legendgroup = comb[[i,3]], showlegend = diff_legend, name = comb[[i,3]])
          }else{
            p <- p %>% 
              add_trace(theta = c(diff$start[j]*180/pi + 180, seq(diff$start[j]*180/pi + 180, diff$end[j]*180/pi + 180, length.out = 20), diff$end[j]*180/pi + 180), 
                        r = c(0, rep(r_range[2], 20), 0), 
                        line = list(color = "black", width = 0.5), fill = "toself", fillcolor = rgb(0,0,0,max=255,alpha=25), legendgroup = comb[[i,3]], showlegend = FALSE)
          }
        }
      }
    }
  }
  
  if(legend_pos %in% c("vertical", "v")){
    legend = list(orientation = 'v', y = 0.5, yanchor = "center")
  }else{
    legend = list(orientation = 'h', x = 0.5, xanchor = "center")
  }
  m <- list(l = 40, r = 40, b = 0, t = 40, pad = 0)
  p <- layout(p, legend = legend, polar = list(sector = c(15,165), radialaxis = list(angle = 15, range = r_range, dtick = dtick), angularaxis = list(thetaunit = "degrees", direction = "clockwise", rotation = 180)), title = title, margin = m)
  
  if(return_data){
    if(plot_diff){
      print(r_diff)
    }
    return(r_data)
  }else{
    return(p)
  }
}