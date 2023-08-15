#' FDsingle(data, distance, tau, q) Functional Diversity of a single site for specified values of tau and q of Hill
#'
#' This function calculates functional diversity for a single sample
#' at a specific value of tau and Hill's q.
#'
#' @param data a vector of species sample frequencies.
#' @param distance a matrix of species-pairwise distances.
#' @param tau a numeric for a specified level of threshold distinctiveness.
#' @param q a numeric for a specified diversity order q.
#' @return a numeric value of Functional Diversity.
#' @export
FDsingle <- function(data, distance, tau, q){
  distance <- as.matrix(distance)
  distance[which(distance>tau,arr.ind = T)] <- tau
  w <- as.vector((1 - distance/tau) %*% data )
  data <- data[w!=0]
  w <- w[w!=0]
  b <- data/w
  np <- sum(data)
  if(q==1){
    exp(sum(-b*w/np*log(w/np)))
  }else{
    (sum(b*(w/np)^q))^(1 / (1-q))
  }
}

#' FDchao(data, distance, tau, q, boot) Functional Diversity of N sites for various values of tau and q
#'
#' This function calculate Functional Diversity of N sites
#' for various values of tau and q
#'
#' @param data a list with N sites; each element of list is species abundances.
#' @param distance a matrix of species-pairwise distances.
#' @param tau a numeric or a vector of levels of threshold distinctiveness.
#' @param q a numeric or a vector of diversity orders; the suggested range for q is [0, 2].
#' @param boot a numeric of number of bootstrap replications.
#' @return two matrices of FD; the first matrix is the q-profile information, the second matrix is the tau profile information.
#' @export
FDchao <- function(data, distance, tau, q, boot){
  EstiBootComm.Func = function(data, distance){n = sum(data)
  distance = as.matrix(distance)
  distance = distance[data!=0, data!=0]
  X = data[data>0]
  f1 <- sum(X == 1) ; f2 <- sum(X == 2)
  f0.hat <- ceiling(ifelse(f2>0, ((n-1)/n)*f1^2/2/f2, ((n-1)/n)*f1*(f1-1)/2))
  C1 = ifelse(f2>0, 1-f1*(n-1)*f1/n/((n-1)*f1+2*f2), 1-f1*(n-1)*(f1-1)/n/((n-1)*(f1-1)+2))
  W <- (1 - C1)/sum(X/n*(1-X/n)^n)
  Prob.hat.Unse <- rep((1-C1)/f0.hat, f0.hat)
  Prob.hat <- X/n*(1-W*(1-X/n)^n)
  Prob <- c(Prob.hat, Prob.hat.Unse)
  F.1 <- sum(distance[, X==1]) ; F.2 <- sum(distance[, X==2])
  F11 <- sum(distance[X==1, X==1]) ; F22 <- sum(distance[X==2, X==2])
  #
  F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
  F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  if (f0.hat==0) {
    d=distance
  } else if (f0.hat==1) {
    random_distance = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_distance*F.0hat, length(X), f0.hat, byrow = T)
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(distance, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    random_distance = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_distance*F.0hat, length(X), f0.hat, byrow = T)

    fo.num = (f0.hat * (f0.hat-1) )/2
    random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)*random_d00
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(distance, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  return(list("pi" = Prob,"distance" = d))
  }
  distance <-  as.matrix(distance)
  out <- as.vector(distance)
  out <- out[out!=0]
  dmin <- min(out)
  if(length(data)!=1){
    tmp <- apply(do.call(cbind,lapply(data, FUN = function(x) x/sum(x))), 1, mean)
    dmean <-  sum ( (tmp %*% t(tmp) ) * distance)
  }else{
    tmp <- data[[1]]/sum(data[[1]])
    dmean <-  sum ( (tmp %*% t(tmp) ) * distance)
  }
  FD.CI = function(data, distance, tau, q){
    qFun = FDsingle(data, distance, tau, q)
    if(boot!=0){
      BT = EstiBootComm.Func(data, distance)
      p_hat = BT[[1]]
      distance_boot = BT[[2]]
      distance_boot <-  as.matrix(distance_boot)
      distance_boot <- replace(distance_boot, distance_boot==0, 10^(-10))
      for (i in seq_len(nrow(distance_boot))) {
        distance_boot[i, i] <- 0
      }
      n=sum(data)
      Boot.X = rmultinom(boot, n, p_hat)
    } else {
      n=data[1]
      Boot.X = t(sapply(p_hat,function(i) rbinom(boot, n, i)))
    }
    qFun_sd = sd(sapply(seq_len(ncol(Boot.X)), function(i) {
      FDsingle(Boot.X[, i], distance_boot, tau, q)
    }))
    LCL = max(0, qFun - qnorm(0.975) * qFun_sd)
    UCL = qFun + qnorm(0.975) * qFun_sd
    a = round(c(qFun, qFun_sd, LCL, UCL), 4)
    a
  }
  Funq <- function(data){
    dminFDforq <- t(sapply(q, FUN = function(q) FD.CI(data, distance, dmin, q) ))
    dmeanFDforq <-t(sapply(q, FUN = function(q) FD.CI(data, distance, dmean, q) ))
    out <- data.frame(rep(q,2), rbind(dminFDforq,dmeanFDforq),rep(c("dmin","dmean"),each=length(q)))
  }
  if(length(data)!=1){
    name = names(data)
    Outputforq <- data.frame(do.call(rbind,lapply(data, Funq)), rep(name, each=2*length(q)), row.names = NULL)
  }else{
    name = names(data)
    Outputforq <- data.frame(Funq(data[[1]]), name, row.names = NULL)
  }
  colnames(Outputforq) <- c("q","estimate", "s.e.", "LCL", "UCL", "tau","site")

  Output <- list(forq = Outputforq)
  return(Output)
}

#'Div(data)  Data.frame for value of diversity required to plot the results
#'
#' This function prepares a "data.frame" object where it takes 5 variables, namely
#' (qEix, TauMin, TauMed, RedFunChao, Color).
#' These variables are required to plot the graphs (using the PlotDiv function)
#' depicting taxonomic, functional, and functional redundancy diversity results.
#'
#' @param data This parameter must necessarily be an object containing the outcome of the FDchao function; the div function will filter these outcomes to prepare the graph.
#' @return data.frame for value of taxonomic diversity, functional diversity and functional redundance.
#' @export
Div <- function(data) {
  data_plot <- data.frame(
    qEix = subset(data$forq, tau == "dmin")[, "q"],
    TauMin = subset(data$forq, tau == "dmin")[, "estimate"],
    TauMed = subset(data$forq, tau == "dmean")[, "estimate"],
    RedFunChao = 1 - (subset(data$forq, tau == "dmean")[, "estimate"] / subset(data$forq, tau == "dmin")[, "estimate"]),
    Color = subset(data$forq, tau == "dmin")[, "site"]
  )
  return(data_plot)
}

#' plotDiv(data, tog, cap) function to plot the value of diversity
#'
#' This is a function to plot the results of the FDchao function
#' filtered by the Div function, generating three graphs:
#' Taxonomic Diversity, Functional Diversity, and Functional Redundancy.
#' @param data This parameter must necessarily be an object containing the outcome of the Div function to prepare the graph.
#' @param tog This parameter is an abbreviation of the word "together" and is a logical object. If it is set to TRUE, the function will plot the three graphs together (side by side).
#' @param cap This parameter is an abbreviation of the word "captions" and is a logical object. If it is set to TRUE, the function will ask you to set the labels for the graph.
#' @return Three graphs, namely (functional diversity, taxonomic diversity, and functional redundancy).
#' @export
plotDiv <- function(data, tog = TRUE,  cap = FALSE){
  plotTax <- ggplot(data)+
    geom_line(aes(x=qEix, y = TauMin, color = Color, linetype = Color), linewidth = 1.2)+
    scale_color_brewer(palette = "Set1") +
    theme(panel.border = element_rect(colour = "black",fill = NA,size = .5), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 12, face = "plain"), axis.text.x = element_text(colour = "black", size = 10), axis.text.y = element_text(colour = "black", size = 10), legend.key = element_rect(fill = "white", colour = "white"), legend.title = element_blank(), panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text())+
    labs(x="Ordem q", y="Diversidade Taxonômica") +
    scale_y_continuous(limits = c(0, ceiling(max(data$TauMin)))) +
    ggtitle("Perfil q (tau = Mínimo)")

  plotFun <- ggplot(data)+
    geom_line(aes(x=qEix, y = TauMed, color = Color, linetype = Color), linewidth = 1.2)+
    scale_color_brewer(palette = "Set1") +
    theme(panel.border = element_rect(colour = "black",fill = NA,size = .5), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 12, face = "plain"), axis.text.x = element_text(colour = "black", size = 10), axis.text.y = element_text(colour = "black", size = 10), legend.key = element_rect(fill = "white", colour = "white"), legend.title = element_blank(), panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text())+
    labs(x="Ordem q", y="Diversidade Funcional") +
    scale_y_continuous(limits = c(0, ceiling(max(data$TauMed)))) +
    ggtitle("Perfil q (tau = Médio)")

  plotRedFun <- ggplot(data)+
    geom_line(aes(x = qEix, y = RedFunChao, color = Color, linetype = Color), linewidth=1.2)+
    scale_color_brewer(palette = "Set1") +
    theme(panel.border = element_rect(colour = "black",fill = NA,size = .5), legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 12, face = "plain"), axis.text.x = element_text(colour = "black", size = 10), axis.text.y = element_text(colour = "black", size = 10), legend.key = element_rect(fill = "white", colour = "white"), legend.title = element_blank(), panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(), panel.background = element_blank(), text = element_text())+
    labs(x="Ordem q", y="Redundância Funcional") +
    scale_y_continuous(limits = c(0, max(data$RedFunChao)*1.1)) +
    ggtitle("Perfil q (Redundância Funcional)")

  if (cap) {
    tax_title <- readline("Digite o título para o gráfico de Diversidade Taxonômica: ")
    tax_x_label <- readline("Digite a etiqueta para o eixo x do gráfico de Diversidade Taxonômica: ")
    tax_y_label <- readline("Digite a etiqueta para o eixo y do gráfico de Diversidade Taxonômica: ")

    fun_title <- readline("Digite o título para o gráfico de Diversidade Funcional: ")
    fun_x_label <- readline("Digite a etiqueta para o eixo x do gráfico de Diversidade Funcional: ")
    fun_y_label <- readline("Digite a etiqueta para o eixo y do gráfico de Diversidade Funcional: ")

    red_title <- readline("Digite o título para o gráfico de Redundância Funcional: ")
    red_x_label <- readline("Digite a etiqueta para o eixo x do gráfico de Redundância Funcional: ")
    red_y_label <- readline("Digite a etiqueta para o eixo y do gráfico de Redundância Funcional: ")

    # Alterar títulos, etiquetas de eixos x e y dos gráficos
    plotTax <- plotTax + ggtitle(tax_title) + labs(x = tax_x_label, y = tax_y_label)
    plotFun <- plotFun + ggtitle(fun_title) + labs(x = fun_x_label, y = fun_y_label)
    plotRedFun <- plotRedFun + ggtitle(red_title) + labs(x = red_x_label, y = red_y_label)
  }
  if (tog) {
    return(ggpubr::ggarrange(plotTax, plotFun, plotRedFun, ncol = 3, nrow = 1))
  } else {
    choice <- menu(c("Plot Diversidade Taxonômica", "Plot Diversidade Funcional", "Plot Redundância Funcional"), title = "Escolha um gráfico:")

    if (choice == 1) {
      return(plotTax)
    } else if (choice == 2) {
      return(plotFun)
    } else if (choice == 3) {
      return(plotRedFun)
    } else {
      stop("Escolha inválida.")
    }
  }
}
