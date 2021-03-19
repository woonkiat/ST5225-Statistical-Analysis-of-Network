
library(igraph)
g <- read_graph("polblogs.gml", format = "gml")
summary(g)
g <- simplify(g)


####### 1 ########

set.seed(888)
plot(g,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(g),vertex.cex=1)

V(g)$color <- vertex_attr(g,"value",index=V(g))+1

#plot the graph with nodes without edge removed
sg <- induced.subgraph(g, which(degree(g,mode="total")>=1))
set.seed(6)
plot(sg,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(sg))



##### 1a #####

diameter(g, directed=TRUE, weights=NA)
get_diameter(g, directed=TRUE, weights=NA)
mean_distance(g, directed=TRUE)
edge_density(g)



##### 1b ######

#degree distribution plot
deg = degree(g, mode="total")
deg.dist <- degree.distribution(g, mode = "total",cumulative = FALSE)
plot( x=0:max(deg), y=deg.dist, cex=0.7, col=1, xlab="Degree", ylab="f(d)")
#hist(degree(g, mode="total"))
cdeg.dist <- degree.distribution(g, mode = "total",cumulative = TRUE)
plot( x=0:max(deg), y=1-cdeg.dist, cex=0.7, col=1, xlab="Degree", ylab="F(d)")


fit_power_law1 = function(graph) {
  d = degree(graph, mode = "total") 
  dd = degree.distribution(graph, mode = "total", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  degree = (0:max(d))+1
  probability = dd
  
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("C =", round(cozf[[1]], 3)))
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "f(d) (log)", 
       col = 1, main = "log(f(d)) versus log(d)")
  
  legend(5,0.09,paste("f(d) = exp(", round(cozf[[2]], 3),"*log(d)", round(cozf[[1]], 3),")"))
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law1(g)


fit_power_law2 = function(graph) {
  d = degree(graph, mode = "total") 
  dd <- 1-degree.distribution(graph, mode = "total", cumulative = TRUE)
  dl <- (0:max(d))+1
  
  reg = lm(log(1-dd) ~ log(dl))
  cozf = coef(reg)
  alpha = -cozf[[2]]+1
  R.square = summary(reg)$r.squared
  print(paste("C =", round(cozf[[1]], 3)))
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  
  plot(dl,1-dd,log = "xy",xlab="Degree (log)", ylab="1-F(d) (log)",col = 1, main = "log(1-F(d)) versus log(d)")
  legend(1,0.002,paste("1-F(d) = exp(", round(cozf[[2]], 3),"*log(d) + ", round(cozf[[1]], 3),")"))
  
  fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  curve(fit, col = "red", add = T, n = length(d))
}
fit_power_law2(g)


#built in function for finding alpha
power.law.fit(degree(g,mode="all")+1)




##### 1c #####

degin = degree(g, mode="in")
degintop <- order(degin, decreasing=TRUE)[1:6]
t1 <- cbind.data.frame(degintop, degin[degintop], V(g)[degintop]$value, V(g)[degintop]$label)
names(t1) <- c("node","in deg","value","label")
t1
#plot(sg,edge.arrow.size=.3,vertex.size=degin/10,vertex.label=NA,layout=layout_with_fr(sg))


degout = degree(g, mode="out")
degouttop <- order(degout, decreasing=TRUE)[1:6]
t2 <- cbind.data.frame(degouttop, degout[degouttop], V(g)[degouttop]$value, V(g)[degouttop]$label)
names(t2) <- c("node","out deg","value","label")
t2


clo = closeness(g,mode="out")*(length(V(g)) - 1)
clotop <- order(clo, decreasing=TRUE)[1:6]
t3 <- cbind.data.frame(clotop, clo[clotop], V(g)[clotop]$value, V(g)[clotop]$label)
names(t3) <- c("node","closeness","value","label")
t3


btw = betweenness(g, directed = TRUE)
btwtop <- order(btw, decreasing=TRUE)[1:6]
t4 <- cbind.data.frame(btwtop, btw[btwtop], V(g)[btwtop]$value, V(g)[btwtop]$label)
names(t4) <- c("node","betweenness","value","label")
t4





########### 2 ###########

hs <- hub.score(g)$vector
hstop <- order(hs, decreasing=TRUE)[1:10]
htable <- cbind.data.frame(hstop, hs[hstop], V(g)[hstop]$value, V(g)[hstop]$label)
htable
#plot(degout,hs,col=V(g)$value+1)

as <- authority.score(g)$vector
astop <- order(as, decreasing=TRUE)[1:10]
atable <- cbind.data.frame(astop, as[astop], V(g)[astop]$value, V(g)[astop]$label)
atable
#plot(degin,as,col=V(g)$value+1)


pr <- page.rank(g, damping = 0.99999)$vector
prtop <- order(pr, decreasing=TRUE)[1:10]
prtable <- cbind.data.frame(prtop, pr[prtop], V(g)[prtop]$value, V(g)[prtop]$label)
prtable

prd <- page.rank(g, damping = 0.85)$vector 
prdtop <- order(prd, decreasing=TRUE)[1:10]
prdtable <- cbind.data.frame(prdtop, prd[prdtop], V(g)[prdtop]$value, V(g)[prdtop]$label)
prdtable



# Page rank -- manual way
A = get.adjacency(g); 
A = as.matrix(A);

deg = degree(g, mode = "out")  

P = diag(1./deg)%*%A; 
#P[is.na(P)] <- 0
#for nodes with zero out degree, modified to point to all nodes! this is the same approach as the built in page.rank function
P[is.na(P)] <- 1/length(V(g))

r = rep(1/length(V(g)), length(V(g)));

for (i in 1:10000) {
  r = r%*%P
  r = r/(sum(r))
}
r[1:10]

names(r) <- V(g)$id
sort(r, decreasing=TRUE)[1:10]

vcolor <- (r>0.1) +1
plot(g,edge.arrow.size=.3,vertex.color =vcolor, vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(g))

sg_dangling <- induced.subgraph(g, c(1293, neighbors(g,1293,"total")))
plot(sg_dangling,edge.arrow.size=0.4, vertex.size=degree(sg_dangling,mode = "in"),vertex.label=V(sg_dangling)$id)

sg_dangling <- induced.subgraph(g, c(1159, neighbors(g,1159,"total")))
plot(sg_dangling,edge.arrow.size=0.4, vertex.size=degree(sg_dangling,mode = "in"),vertex.label=V(sg_dangling)$id)



# Page rank with damping -- manual way
r = rep(1/length(V(g)), length(V(g)));
l = length(r)

for (i in 1:5000) {
  r = r%*%P
  r = r/(sum(r))
  r = 0.85*r + 0.15*(1/l)
}
r[1:10]

names(r) <- V(g)$id
sort(r, decreasing=TRUE)[1:10]






######################3

ug = as.undirected(g)
ug = simplify(ug)

vcount(ug)
ecount(ug)
diameter(ug, directed=FALSE, weights=NA)
get_diameter(ug, directed=FALSE, weights=NA)
mean_distance(ug, directed=FALSE)
edge_density(ug)
is_connected(ug)
transitivity(ug, type = "global")


deg = degree(ug, mode="total")
deg.dist <- degree.distribution(ug, mode = "total",cumulative = FALSE)
plot( x=0:max(deg), y=deg.dist, cex=0.7, col=1, xlab="Degree", ylab="f(d)")

cdeg.dist <- degree.distribution(ug, mode = "total",cumulative = TRUE)
plot( x=0:max(deg), y=1-cdeg.dist, cex=0.7, col=1, xlab="Degree", ylab="F(d)",xlim=c(0,400))


#fitting regression with f(d) and d
fit_power_law1(ug)

#fitting regression with 1-F(d) and d
fit_power_law2(ug)

#built in function for finding alpha
power.law.fit(degree(ug,mode="all")+1)



is.connected(ug)
clu <- components(ug) #return the nodes in each component
clu
groups(clu)


comp <- decompose.graph(ug) #return the components, each as a graph
# Select the giant component
giantIndex <- which.max(sapply(comp, vcount))
GiantComp <- comp[[giantIndex]]; 
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))


coreness(ug)
maxk = max(coreness(ug))
kvertex = which(coreness(ug) == maxk)

vsize = (coreness(ug)+1)/10
plot(ug,edge.arrow.size=.3,vertex.size=vsize,vertex.label=NA,layout=layout_with_graphopt(ug))

plot(ug,edge.arrow.size=.3,vertex.size=vsize,vertex.label=NA,layout=layout_nicely)

#induced a subgraph of kcore 
#kcore <- induced.subgraph(ug, which(degree(ug)>= maxk))
kcore <- induced.subgraph(ug, kvertex)
plot(kcore,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(kcore),asp=9/9,vertex.label.dist=1)
vcount(kcore)

################################4

### 1. Induced-Subgraph Sampling
set.seed(525)
sv <- sample(1:vcount(ug), 500)
sg1 <- induced.subgraph(ug, sv)
plot(sg1,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(sg1),asp=9/9,vertex.label.dist=1)


#par(mai=c(0,0,1,0), mfrow = c(1, 2))
ug1 <- ug
V(ug1)$prop <- "pop" #Now all of the nodes are not chosen
V(ug1)[sv]$prop <- "sample" # Some of them are chosen as sample
V(ug1)$color <- ifelse(V(ug1)$prop == "sample", "red", "grey")

ev <- get.edges(ug1, E(ug1)); 
se <- V(ug1)[ev[,1]]$prop == "sample" & V(ug1)[ev[,2]]$prop == "sample"; 
E(ug1)$color <- ifelse(se, "red", "grey")
set.seed(888)
plot(ug1,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(ug1),asp=9/9,vertex.label.dist=1)



### 2. Incident-Subgraph Sampling
set.seed(525)
ne = 50
nv = 0
e <- get.edges(ug, E(ug));
while (nv < 500) {
  ne <- ne+1    #add one edge if nodes less than 500
  se <- sample(1:length(E(ug)), ne)    #sample edges
  sg2 <- graph.edgelist(e[se,], directed=FALSE)     #create a graph based on sampled edges,however if edge start with, say 3, node 1 and 2 will also be created, remove unnesassary nodes with next line
  V(sg2)$color <- V(ug)$color[1:vcount(sg2)]
  sg2 <- induced.subgraph(sg2, which(degree(sg2,mode="total")>=1)) 
  nv <- vcount(sg2)
}
plot(sg2,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(sg2),asp=9/9,vertex.label.dist=1)
print(nv)


ug2 <- ug
col_edges <- rep("grey", length(E(ug2)))
col_edges[se] <- "red"
E(ug2)$color <- col_edges
plot(ug2, layout = layout_with_graphopt(ug1))

col_nodes <- rep("grey", length(V(ug2)));
col_nodes[e[se,1]] <- "red"
col_nodes[e[se,2]] <- "red"
V(ug2)$color <- col_nodes
set.seed(888)
plot(ug2,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(ug2),asp=9/9,vertex.label.dist=1)





### 3. Snowball Sampling
set.seed(525)
sv<- sample(1:vcount(ug), 5)
i = 1
count = 5

while (count < 500) {
  current <- sv[i]
  sv2 <- neighbors(ug, current)
  
  for (j in sv2) {
    if (is.element(j,sv) == FALSE) {
      sv <- c(sv, j)
      count = count+1
      if (count >= 500) {
        break
      }
  }
  }
  i = i + 1
}

sg3 <- induced.subgraph(ug, sv)
plot(sg3,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(sg3),asp=9/9,vertex.label.dist=1)

ug3 <- ug
V(ug3)$prop <- "pop" #Now all of the nodes are not chosen
V(ug3)[sv]$prop <- "sample" # Some of them are chosen as sample
V(ug3)$color <- ifelse(V(ug3)$prop == "sample", "red", "grey")

ev <- get.edges(ug3, E(ug3)); 
se <- V(ug3)[ev[,1]]$prop == "sample" & V(ug3)[ev[,2]]$prop == "sample"; 
E(ug3)$color <- ifelse(se, "red", "grey")
set.seed(888)
plot(ug3,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(ug3),asp=9/9,vertex.label.dist=1)




### 4. Respondent Driven Sampling

set.seed(525)
sv<- sample(1:vcount(ug), 5)
k = 3
count = 5
i = 1
se <- vector()
se <- na.omit(se)

while (count < 500) {
  current <- sv[i]
  sv2 <- neighbors(ug, current)
  ki <- 0 
  for (j in sv2) {
    if (is.element(j,sv) == FALSE) {
      sv <- c(sv, j)
      se <- rbind(se,c(current,j))
      count = count+1
      ki = ki+1
      if (count >= 500 || ki >= 3) {
        break
      }
    }
  }
  i = i + 1
}

sg4 <- graph.edgelist(se, directed=FALSE)
V(sg4)$color <- V(ug)$color[1:vcount(sg4)]
sg4 <- induced.subgraph(sg4, which(degree(sg4,mode="total")>=1))
plot(sg4,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(sg4),asp=9/9,vertex.label.dist=1)

seid <- get.edge.ids(ug, array(t(se)))

ug4 <- ug
V(ug4)$prop <- "pop" 
V(ug4)[sv]$prop <- "sample" # Some of them are chosen as sample
V(ug4)$color <- ifelse(V(ug4)$prop == "sample", "red", "grey")

E(ug4)$color <- "grey"
E(ug4)[seid]$color <- "red"
set.seed(888)
plot(ug4,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_graphopt(ug4),asp=9/9,vertex.label.dist=1)


ecount <- c(ecount(sg1),ecount(sg2),ecount(sg3),ecount(sg4))  
ncount <- c(vcount(sg1),vcount(sg2),vcount(sg3),vcount(sg4))  
density_table <-  c(edge_density(sg1),edge_density(sg2),edge_density(sg3),edge_density(sg4))
transitivity_table <- c(transitivity(sg1, type = "global"),transitivity(sg2, type = "global"),transitivity(sg3, type = "global"),transitivity(sg4, type = "global"))
is_connected_table <- c(is.connected(sg1),is.connected(sg2),is.connected(sg3),is.connected(sg4))

sampling_table <- cbind.data.frame(ncount,ecount,density_table,transitivity_table,is_connected_table)
sampling_table

groups(components(sg4))








############ 5

gc <- GiantComp

'''
ebc <- edge.betweenness.community(gc)
plot(ebc, GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))

cut <- cutat(ebc,2)
colors <- rainbow(2)
plot(gc, vertex.color=colors[cut],,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(gc))
'''

set.seed(525)
##edge betweenness
while (components(gc)$no <2) {
eb <- edge_betweenness(gc, e = E(gc), directed = FALSE)
maxid <- which(eb==max(eb))
print(maxid)
gc <- gc - E(gc)[maxid]
}

V(GiantComp)$color <- components(gc)$membership
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))

ebresult <- table(components(gc)$membership-1, V(gc)$value)
ebresult
err_rate0 = min(sum(ebresult[1,2],ebresult[2,1])/sum(ebresult), 1-sum(ebresult[1,2],ebresult[2,1])/sum(ebresult))
err_rate0


A = get.adjacency(GiantComp) #get the adjacency matrix

graphdist <- function(A){
  n = dim(A)[1]; p = dim(A)[2]
  distmatrix = matrix(0, nrow = n, ncol = n);
  for(i in 1:p){
    #print(i)
    for(j in 1:i){
      distmatrix[i,j] = sqrt(sum((A[i,-c(i,j)] - A[j,-c(i,j)])^2));
    }
  }
  return(distmatrix)
}

d = graphdist(A) #calculate the Euclidean distance
d = as.dist(d);

  
#par(mfrow = c(1,3))
clu1 = hclust(d, method = "complete")
clu2 = hclust(d, method = "single")
clu3 = hclust(d, method = "average")
plot(clu1)
plot(clu2)
plot(clu3)

memb1 <- cutree(clu1, k = 2)
result1 <- table(memb1-1, V(gc)$value)
err_rate1 = min(sum(result1[1,2],result1[2,1])/sum(result1), 1-sum(result1[1,2],result1[2,1])/sum(result1))

V(GiantComp)$color <- memb1
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))
result1
err_rate1


memb2 <- cutree(clu2, k = 2)
result2 <- table(memb2-1, V(gc)$value)
err_rate2 = min(sum(result2[1,2],result2[2,1])/sum(result2), 1-sum(result2[1,2],result2[2,1])/sum(result2))

V(GiantComp)$color <- memb2
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))
result2
err_rate2


memb3 <- cutree(clu3, k = 2)
result3 <- table(memb3-1, V(gc)$value)
err_rate3 = min(sum(result3[1,2],result3[2,1])/sum(result3), 1-sum(result3[1,2],result3[2,1])/sum(result3))

V(GiantComp)$color <- memb3
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))
result3
err_rate3


kc <- fastgreedy.community(GiantComp)
kc2 = cut_at(kc, no = 2); 

result4 <- table(kc2-1, V(gc)$value)
err_rate4 = min(sum(result4[1,2],result4[2,1])/sum(result4), 1-sum(result4[1,2],result4[2,1])/sum(result4))
V(GiantComp)$color <- kc2
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))
result4
err_rate4

modularity(GiantComp, kc2, weights = NULL)
modularity(GiantComp, V(GiantComp)$value+1, weights = NULL)



set.seed(525)
r <- eigen(A)$vector[,2]/eigen(A)$vector[,1]

km <- kmeans(r, centers = 2)

V(GiantComp)$color <- km$cluster
plot(GiantComp,edge.arrow.size=.3,vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))

result5 <- table(km$cluster-1, V(gc)$value)
err_rate5 = min(sum(result5[1,2],result5[2,1])/sum(result5), 1-sum(result5[1,2],result5[2,1])/sum(result5))
result5
err_rate5



'''
plot(eigen(A)$vector[,1],eigen(A)$vector[,2],col=V(gc)$value+1,pch=V(gc)$value+1)
legend(0.125,-0.1,paste(c("0","1")), col=1:2,pch = 1:2)


plot(eigen(A)$vector[,1],eigen(A)$vector[,2],col=(V(gc)$value == km$cluster-1) +1,pch=V(gc)$value+1)
plot(eigen(A)$vector[,1], eigen(A)$vector[,2]/eigen(A)$vector[,1], pch=V(gc)$value+1, col = (V(gc)$value == km$cluster-1) +1)
'''



el = get.edges(GiantComp,E(GiantComp))

#member <- km$cluster
member <- V(GiantComp)$value

c <- matrix(0,ecount(GiantComp),2)

for (i in 1:ecount(GiantComp)){
  c[i,1] = member[el[i,1]]
  c[i,2] = member[el[i,2]]
}

e00 <- sum((c[,1] == 0) * (c[,2] == 0))
e01 <- sum((c[,1] == 0) * (c[,2] == 1)) + sum((c[,1] == 1) * (c[,2] == 0))
e11 <- sum((c[,1] == 1) * (c[,2] == 1))
n0 <- sum(V(GiantComp)$value == 0)
n1 <- sum(V(GiantComp)$value == 1)

b00 <- e00/n0/(n0-1)*2
b11 <- e11/n1/(n1-1)*2
b01 <- e01/n0/n1
b00
b01
b11

set.seed(525)
V(GiantComp)$color <- vertex_attr(GiantComp,"value",index=V(GiantComp))+1
ecolor <- apply(c,1,function(c) {
  if ((c[1] == 0) && (c[2] == 0)) {1} 
  else if ((c[1] == 1) && (c[2] == 1)) {2} 
  else {"red"}}
  )

plot(GiantComp,edge.arrow.size=.3,edge.color =ecolor, vertex.size=2,vertex.label=NA,layout=layout_with_fr(GiantComp))







