## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})
##
### Load packages here:
##
### Figure counter:
(
  function() {
    log <- list(
      labels = character(),
      captions = character()
    )
    list(
      register = function(label, caption) {
        log$labels <<- c(log$labels, label)
        log$captions <<- c(log$captions, caption)
        invisible(NULL)
      },
      getNumber = function(label) {
        which(log$labels == label)
      },
      getCaption = function(label) {
        a <- which(log$labels == label)
        cap <- log$captions[a]
        cat(sprintf("Fig. %d. %s\n\n---\n",a,cap))
        invisible(NULL)
      }
    )
  }
)() -> figCounter

## ----load_package-------------------------------------------------------------
library(MPSEM)

## ----load_data----------------------------------------------------------------
data(perissodactyla,package="caper")

## ----plot_phylogeny, echo=FALSE, fig.height=4, fig.width = 7------------------
par(mar=c(2,2,2,2))
plot(perissodactyla.tree)
par(mar=c(5,4,4,2))
figCounter$register(
  "theTree",
  "The phylogenetic tree used for this example."
)

## ----plot_phylogeny_cap, echo=FALSE, results='asis'---------------------------
figCounter$getCaption("theTree")

## ----data_table, echo=FALSE, results="latex"----------------------------------
knitr::kable(perissodactyla.data[,c(1L,2L,4L)])

## ----droping_species----------------------------------------------------------
match(
  perissodactyla.tree$tip.label,
  perissodactyla.data[,1L]
) -> spmatch

drop.tip(
  perissodactyla.tree,
  perissodactyla.tree$tip.label[is.na(spmatch)]
) -> perissodactyla.tree

## ----check_order--------------------------------------------------------------
cbind(perissodactyla.tree$tip.label, perissodactyla.data[,1L])

## ----re-order_species---------------------------------------------------------
match(
  perissodactyla.tree$tip.label,
  perissodactyla.data[,1L]
) -> spmatch

perissodactyla.data[spmatch,] -> perissodactyla.data

all(perissodactyla.tree$tip.label == perissodactyla.data[,1L])

## ----change_rownames----------------------------------------------------------
perissodactyla.data[,1L] -> rownames(perissodactyla.data)
perissodactyla.data[,-1L] -> perissodactyla.data

## ----re-arranged_data, echo=FALSE, results="latex"----------------------------
knitr::kable(perissodactyla.data[,c(1L,3L)])

## ----training_testing_datasets------------------------------------------------
perissodactyla.data[-1L,,drop=FALSE] -> perissodactyla.train
perissodactyla.data[1L,,drop=FALSE] -> perissodactyla.test
drop.tip(
  perissodactyla.tree,
  tip = "Dicerorhinus sumatrensis"
) -> perissodactyla.tree.train

## ----display_weighting, echo=FALSE, fig.height=5, fig.width = 7---------------
par(mar=c(4.5,4.5,1,7) + 0.1)
d <- seq(0, 2, length.out=1000)
a <- c(0,0.33,0.67,1,0.25,0.75,0)
psi <- c(1,1,1,1,0.65,0.65,0.4)
cc <- c(1,1,1,1,1,1,1)
ll <- c(1,2,2,2,3,3,3)
trial <- cbind(a, psi)
colnames(trial) <- c("a","psi")
ntrials <- nrow(trial)
nd <- length(d)
matrix(
  NA,
  ntrials,
  nd,
  dimnames=list(paste("a=", trial[,"a"], ", psi=", trial[,"psi"], sep=""),
                paste("d=", round(d,4), sep=""))
) -> w
for(i in 1:ntrials)
  w[i,] <- MPSEM::PEMweights(d, trial[i,"a"], trial[i,"psi"])
plot(NA, xlim=c(0,2), ylim=c(0,1.6), ylab="Weight", xlab="Distance", axes=FALSE)
axis(1, at=seq(0,2,0.5), label=seq(0,2,0.5))
axis(2, las=1)
text(expression(paste(~~~a~~~~~~~psi)),x=2.2,y=1.57,xpd=TRUE,adj=0)
for(i in 1:ntrials) {
  lines(x=d, y=w[i,], col=cc[i], lty=ll[i])
  text(paste(sprintf("%.2f", trial[i,1]), sprintf("%.2f",trial[i,2]), sep="  "),
       x=rep(2.2,1), y=w[i,1000], xpd=TRUE, adj=0)
}
rm(d,a,psi,cc,ll,trial,ntrials,nd,w,i)
figCounter$register(
  "edgeWeighting",
  paste(
    "Output of the edge weighting function for different sets of parameters",
    "$a$ and $\\psi$."
  )
)

## ----display_weighting_cap, echo=FALSE, results='asis'------------------------
figCounter$getCaption("edgeWeighting")

## ----convert_to_graph---------------------------------------------------------
Phylo2DirectedGraph(
  perissodactyla.tree.train
) -> perissodactyla.pgraph

## ----graph_storage,echo=FALSE-------------------------------------------------
str(perissodactyla.pgraph)

## ----tree_labelled, fig.height=5, fig.width = 7-------------------------------
perissodactyla.tree.train -> tree
paste("N",1L:tree$Nnode) -> tree$node.label

par(mar=c(2,2,2,2))
plot(tree,show.node.label=TRUE)

edgelabels(
  1L:nrow(tree$edge),
  edge=1L:nrow(tree$edge),
  bg="white",
  cex=0.75
)

## ----tree_labelled_cap, echo=FALSE, results='asis'----------------------------
figCounter$register(
  "trainingTree",
  "The labelled training species tree for this example."
)
figCounter$getCaption("trainingTree")

## ----set_param----------------------------------------------------------------
rep(0,attr(perissodactyla.pgraph,"ev")[1L]) -> steepness
rep(1,attr(perissodactyla.pgraph,"ev")[1L]) -> evol_rate

steepness[15L:21] <- 0.25
evol_rate[15L:21] <- 2
steepness[9L:13] <- 0.8
evol_rate[9L:13] <- 0.5

## ----calculate_PEM------------------------------------------------------------
PEM.build(
  perissodactyla.pgraph,
  d="distance",
  sp="species",
  a=steepness,
  psi=evol_rate
) -> perissodactyla.PEM

## ----Eigenvector_example, fig.height=4, fig.width=7.5-------------------------
layout(matrix(c(1,1,1,2,2,3,3),1L,7L))
par(mar=c(5.1,2.1,4.1,2.1))

## Singular vectors are extracted using the as.data.frame method:
as.data.frame(perissodactyla.PEM) -> perissodactyla.U

plot(perissodactyla.tree.train, x.lim=60, cex=1.5)
plot(y = 1L:nrow(perissodactyla.train), ylab="", xlab = "Loading",
     x = perissodactyla.U[,1L], xlim=0.5*c(-1,1),
     axes=FALSE, main = expression(bold(v)[1]), cex=1.5)
axis(1)
abline(v=0)

plot(y = 1L:nrow(perissodactyla.train), ylab="", xlab = "Loading",
     x = perissodactyla.U[,5L], xlim=0.5*c(-1,1),
     axes=FALSE, main = expression(bold(v)[5]), cex=1.5)
axis(1)
abline(v=0)

## ----Eigenvector_example_cap, echo=FALSE, results='asis'----------------------
figCounter$register(
  "eigenvectorExample",
  "Example of two eigenvectors obtained from the training species phylogeny."
)
figCounter$getCaption("eigenvectorExample")

## ----PEM_opt1-----------------------------------------------------------------
PEM.fitSimple(
  y = perissodactyla.train[,"log.neonatal.wt"],
  x = NULL,
  w = perissodactyla.pgraph,
  d = "distance",
  sp="species",
  lower = 0,
  upper = 1
) -> perissodactyla.PEM_opt1

## ----PEM_opt2-----------------------------------------------------------------
PEM.fitSimple(
  y = perissodactyla.train[,"log.neonatal.wt"],
  x = perissodactyla.train[,"log.female.wt"],
  w = perissodactyla.pgraph,
  d = "distance",
  sp="species",
  lower = 0,
  upper = 1
) -> perissodactyla.PEM_opt2

## ----build_PEM_models---------------------------------------------------------
lmforwardsequentialAICc(
  y = perissodactyla.train[,"log.neonatal.wt"],
  object = perissodactyla.PEM_opt1
) -> lm1
summary(lm1)

lmforwardsequentialAICc(
  y = perissodactyla.train[,"log.neonatal.wt"],
  x = perissodactyla.train[,"log.female.wt",drop=FALSE],
  object = perissodactyla.PEM_opt2
) -> lm2
summary(lm2) 

## ----make_prediction----------------------------------------------------------
getGraphLocations(
  perissodactyla.tree,
  targets = "Dicerorhinus sumatrensis"
) -> perissodactyla.loc

predict(
  object = perissodactyla.PEM_opt2,
  targets = perissodactyla.loc,
  lmobject = lm2,
  newdata = perissodactyla.test,
  interval = "prediction",
  level = 0.95) -> pred

## ----cross-validation---------------------------------------------------------
data.frame(
  perissodactyla.data,
  predictions = NA,
  lower = NA,
  upper = NA
) -> perissodactyla.data

jackinfo <- list()
for(i in 1L:nrow(perissodactyla.data)) {
  jackinfo[[i]] <- list()
  getGraphLocations(
    perissodactyla.tree,
    targets = rownames(perissodactyla.data)[i]
  ) -> jackinfo[[i]][["loc"]]
  PEM.fitSimple(
    y = perissodactyla.data[-i,"log.neonatal.wt"],
    x = perissodactyla.data[-i,"log.female.wt"],
    w = jackinfo[[i]][["loc"]]$x
  ) -> jackinfo[[i]][["PEM"]]
  lmforwardsequentialAICc(
    y = perissodactyla.data[-i,"log.neonatal.wt"],
    x = perissodactyla.data[-i,"log.female.wt",drop=FALSE],
    object = jackinfo[[i]][["PEM"]]
  ) -> jackinfo[[i]][["lm"]]
  predict(
    object = jackinfo[[i]][["PEM"]],
    targets = jackinfo[[i]][["loc"]],
    lmobject = jackinfo[[i]][["lm"]],
    newdata = perissodactyla.data[i,"log.female.wt",drop=FALSE],
    interval = "prediction",
    level = 0.95
  ) -> predictions
  unlist(predictions) -> perissodactyla.data[i, c("predictions", "lower", "upper")]
}
rm(i, predictions)

## ----plot_pred_obs, echo=FALSE, fig.height=7, fig.width=7---------------------
par(mar=c(5,5,2,2)+0.1)
range(
  perissodactyla.data[,"log.neonatal.wt"],
  perissodactyla.data[,c("predictions","lower","upper")]
) -> rng

plot(NA, xlim = rng, ylim = rng, xlab = "observed", ylab = "Predicted",
     asp = 1, las = 1)
points(
  x = perissodactyla.data[,"log.neonatal.wt"],
  y = perissodactyla.data[,"predictions"]
)

abline(0,1)
arrows(
  x0 = perissodactyla.data[,"log.neonatal.wt"],
  x1 = perissodactyla.data[,"log.neonatal.wt"],
  y0 = perissodactyla.data[,"lower"],
  y1 = perissodactyla.data[,"upper"],
  length = 0.05,
  angle = 90,
  code = 3
)

figCounter$register(
  "crossPreds",
  paste(
    "Leave-one-out crossvalidated prediction of the neonatal weight for",
    nrow((perissodactyla.data)),
    "odd-toed ungulate species."
  )
)

## ----plot_pred_obs_cap, echo=FALSE, results='asis'----------------------------
figCounter$getCaption("crossPreds")

## ----influence_matrix---------------------------------------------------------
PEMInfluence(perissodactyla.pgraph) -> res
res

## ----PEM_updater--------------------------------------------------------------
PEM.updater(object = perissodactyla.PEM, a = 0, psi = 1) -> res
res

## ----forcedSimple-------------------------------------------------------------
PEM.forcedSimple(
  y = perissodactyla.train[,"log.neonatal.wt"],
  x = perissodactyla.train[,"log.female.wt"],
  w = perissodactyla.pgraph,
  a = steepness,
  psi = evol_rate
) -> res
res

## ----get_scores---------------------------------------------------------------
Locations2PEMscores(
  object = perissodactyla.PEM_opt2,
  gsc = perissodactyla.loc
) -> scores
scores

