library(ggplot2)
library(grid)
library(colorspace)

pdata <- read.csv("anticrispr_crispr_cas_res_pm_genome_gc_size_v2.csv")
#pdata <- cbind(pdata, data.frame(zxy = rep(0,dim(pdata)[1])))

ab_names <- names(pdata)[48:63]

cas_count <- apply(pdata[,28:44],1,sum) 
pdata <- cbind(pdata, data.frame(CAS_COUNT = cas_count))

acr_count <- apply(pdata[,c(5:20,24)],1,sum) 
pdata <- cbind(pdata, data.frame(ACR_COUNT = acr_count))

short_name <- function(str) {
	ls <- strsplit(str, " ")[[1]]
	
	if (ls[1] == "Candidatus") {
		strout <- paste(ls[1:3], collapse = " ")
		return(strout)
	} else if (ls[2] %in% c("sp.", "Sp.")) {
		strout <- paste(ls[1:3], collapse = " ")
		return(strout)
	} else if ((ls[1] == "Candidatus") && (ls[2] %in% c("sp.", "Sp."))) {
		strout <- paste(ls[1:4], collapse = " ")
		return(strout)
	} else {
		strout <- paste(ls[1:2], collapse = " ")
		return(strout)
	}
}

genera_name <- function(str) {
		ls <- strsplit(str, " ")[[1]]
		
		if (ls[1] == "Candidatus") {
			strout <- paste(ls[1:2], collapse = " ")
			return(strout)
		} else {
			strout <- ls[1]
			return(strout)
		}
}

do_oe_ab_acr <- function(ab, org) {
	x <- as.numeric(pdata$ACR_COUNT > 0)
	y <- as.numeric(pdata[[ab]])
	
	df <- data.frame(X = x, Y = y)
	
	if (org != "All") {
		inds <- grep(org, as.character(pdata$ORGANISM))
		df <- df[inds,]
	}
	
	n <- dim(df)[1]
	
	fexp <- n * mean(df$X) * mean(df$Y)
	fobs <- sum(df$X * df$Y)
	xfobs <- sum((1 - df$X) * df$Y)
	
	oe <- (fobs + 1e-6)  / (fexp + 1e-6)
	
	qt <- oe_simul_bt(as.numeric(df$X), as.numeric(df$Y))
	prop <- sum(df$X) / n
	
	cest <- "No Assoc."
	if (round(qt[1],5) > 0) {
		cest <- "Assoc."
	}
	
	if (round(qt[2],5) < 0) {
		cest <- "Assoc."
	}
	
	cat <- "(-0.01,0.01]"
	if (log(oe) > 0.01 && log(oe) <= 0.5) {
		cat <- "(0.01,0.5]"
	} else if (log(oe) > 0.5 && log(oe) <= 1.0) {
		cat <- "(0.5,1.0]"
	} else if (log(oe) > 1.0) {
		cat <- "1.0+"
	} else if (log(oe) > -0.5 && log(oe) <= -0.01) {
		cat <- "(-0.5,-0.01]"
	} else if (log(oe) > -1.0 && log(oe) <= -0.5) {
		cat <- "(-1.0,-0.5]"
	} else if (log(oe) <= -1.0) {
		cat <- "< -1.0"
	} else {
		cat <- "(-0.01,0.01]"
	}
	
	rsdf <- data.frame(AB = c(ab), ORG = c(org), N_GENOMES = c(n), 
	                   N_AB = c(sum(df$Y)), N_ACR = c(sum(df$X)), N_LACK_ACR = c(sum(1 - df$X)), 
	                   N_EXCLUDE_OBS = c(xfobs), N_COEXIST_OBS = c(fobs), N_COEXIST_EXP = c(fexp), 
	                   N_COEXIST_OBS_PSC = c(fobs + 1e-6), N_COEXIST_EXP_PSC = c(fexp + 1e-6), OE = c(oe), 
	                   LOG_OE = c(log(oe)), LOG_OE_CAT = c(cat), QTL_0005 = c(qt[1]), 
	                   QTL_9995 = c(qt[2]), SIG = c(cest), PROP = c(prop))
	print(c(as.character(ab),as.character(org)))
	return(rsdf)
}

do_oe_ab_crispr <- function(ab, org) {
	x <- as.numeric(pdata$CAS_COUNT > 0)
	y <- as.numeric(pdata[[ab]])

	df <- data.frame(X = 1 - x, Y = y)
	
	if (org != "All") {
		inds <- grep(org, as.character(pdata$ORGANISM))
		df <- df[inds,]
	}
	
	n <- dim(df)[1]
	
	fexp <- n * mean(df$X) * mean(df$Y)
	fobs <- sum(df$X * df$Y)
	xfobs <- sum((1 - df$X) * df$Y)
	
	oe <- (fobs + 1e-6)  / (fexp + 1e-6)
	
	qt <- oe_simul_bt(as.numeric(df$X), as.numeric(df$Y))
	prop <- sum(1 - df$X) / n
	
	cest <- "No Assoc."
	if ((round(qt[1],5) > 0) && (abs(log(oe)) > 0.2)) {
		cest <- "Assoc."
	}
	
	if ((round(qt[2],5) < 0) && (abs(log(oe)) > 0.2)) {
		cest <- "Assoc."
	}
	
	cat <- "(-0.2,0.2]"
	if (log(oe) > 0.2 && log(oe) <= 0.5) {
		cat <- "(0.2,0.5]"
	} else if (log(oe) > 0.5 && log(oe) <= 1.0) {
		cat <- "(0.5,1.0]"
	} else if (log(oe) > 1.0) {
		cat <- "1.0+"
	} else if (log(oe) > -0.5 && log(oe) <= -0.2) {
		cat <- "(-0.5,-0.2]"
	} else if (log(oe) > -1.0 && log(oe) <= -0.5) {
		cat <- "(-1.0,-0.5]"
	} else if (log(oe) <= -1.0) {
		cat <- "< -1.0"
	} else {
		cat <- "(-0.2,0.2]"
	}
	
	rsdf <- data.frame(AB = c(ab), ORG = c(org), N_GENOMES = c(n), 
	                   N_AB = c(sum(df$Y)), N_CR = c(sum(df$X)), N_LACK_CR = c(sum(df$X)), 
	                   N_COEXIST_OBS = c(xfobs), N_EXCLUDE_OBS = c(fobs), N_EXCLUDE_EXP = c(fexp), 
	                   N_EXCLUDE_OBS_PSC = c(fobs + 1e-6), N_EXCLUDE_EXP_PSC = c(fexp + 1e-6), OE = c(oe), 
	                   LOG_OE = c(log(oe)), LOG_OE_CAT = c(cat), QTL_0005 = c(qt[1]), 
	                   QTL_9995 = c(qt[2]), SIG = c(cest), PROP = c(prop))
	print(c(as.character(ab),as.character(org)))
	return(rsdf)
}

oe_simul <- function(vec1, vec2) {
	lst_rnd_vec_1 <- lapply(1:5000, FUN=function(x){rbinom(length(vec1), 1, mean(vec1))})
	lst_rnd_vec_2 <- lapply(1:5000, FUN=function(x){rbinom(length(vec2), 1, mean(vec2))})
	
	log_oe <- function(x, y) {
		n <- length(x)
		fexp <- n * mean(x) * mean(y)
		fobs <- sum(x * y)
		oe <- (fobs + 1e-6) / (fexp + 1e-6)
		log_oe <- log(oe)
		return(log_oe)
	}
	
	vals <- unlist(lapply(1:length(lst_rnd_vec_1), FUN = function(x){log_oe(lst_rnd_vec_1[[x]], lst_rnd_vec_2[[x]])}))
	qt <- as.numeric(quantile(vals, prob=c(0.0005,0.9995)))
	return(qt)
}

oe_simul_bt <- function(vec1, vec2) {
	lst_inds <- lapply(1:5000, FUN=function(x){sample(1:length(vec1), replace = TRUE)})
	
	log_oe <- function(x, y, inds) {
		xx <- x[inds]
		yy <- y[inds]
		n <- length(x)
		fexp <- n * mean(xx) * mean(yy)
		fobs <- sum(xx * yy)
		oe <- (fobs + 1e-6) / (fexp + 1e-6)
		log_oe <- log(oe)
		return(log_oe)
	}
	
	vals <- unlist(lapply(1:length(lst_inds), FUN = function(x){log_oe(vec1, vec2, lst_inds[[x]])}))
	qt <- as.numeric(quantile(vals, prob=c(0.0005,0.9995)))
	return(qt)
}

orgs_other <- c("Staphylococcus aureus", "Salmonella enterica", "Streptococcus pneumoniae",
                "Escherichia coli", "Mycobacterium tuberculosis", "Klebsiella pneumoniae",
                "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Mycobacterium abscessus",
                "Listeria monocytogenes", "Shigella sonnei", "Campylobacter jejuni",
                "Streptococcus suis", "Neisseria meningitidis", "Streptococcus agalactiae",
                "Vibrio parahaemolyticus", "Vibrio cholerae", "Enterococcus faecium",
                "Clostridioides difficile", "Burkholderia pseudomallei", "Helicobacter pylori",
                "Bordetella pertussis", "Enterobacter cloacae", "Legionella pneumophila",
                "Enterococcus faecalis", "Bacillus cereus", "Neisseria gonorrhoeae",
                "Staphylococcus epidermidis", "Streptococcus pyogenes", "Enterobacter hormaechei",
                "Serratia marcescens", "Campylobacter coli", "Yersinia pestis",
                "Leptospira interrogans", "Burkholderia ubonensis", "Pseudomonas syringae",
                "Burkholderia cenocepacia", "Streptococcus equi", "Ralstonia solanacearum")
                
orgs_other <- orgs_other[order(orgs_other)]                                
orgs <- c(orgs_other)
ab_names_rep <- rep(ab_names, length(orgs))
orgs_rep <- rep(orgs, each = length(ab_names))

org_set2 <- c("Pseudomonas aeruginosa", "Listeria monocytogenes", "Neisseria meningitidis", "Ralstonia solanacearum")
ab_names_rep_2 <- rep(ab_names, length(org_set2))
orgs_rep_2 <- rep(org_set2, each = length(ab_names))

org_set3 <- c("Staphylococcus", "Salmonella", "Streptococcus",
              "Escherichia", "Mycobacterium", "Klebsiella",
              "Acinetobacter", "Pseudomonas", "Mycobacterium",
              "Listeria", "Shigella", "Campylobacter",
              "Streptococcus", "Neisseria", "Streptococcus",
              "Vibrio", "Vibrio", "Enterococcus",
              "Clostridioides", "Burkholderia", "Helicobacter",
              "Bordetella", "Enterobacter", "Legionella",
              "Enterococcus", "Bacillus", "Neisseria",
              "Staphylococcus", "Streptococcus", "Enterobacter",
              "Serratia", "Campylobacter", "Yersinia",
              "Leptospira", "Burkholderia", "Pseudomonas",
              "Burkholderia", "Streptococcus", "Ralstonia")
org_set3 <- unique(org_set3)
ab_names_rep_3 <- rep(ab_names, length(org_set3))
orgs_rep_3 <- rep(org_set3, each = length(ab_names))

df_oe_org <- do.call(rbind, lapply(1:length(ab_names_rep_2), FUN=function(x){do_oe_ab_acr(ab_names_rep_2[x], orgs_rep_2[x])}))
df_oe_org_2 <- do.call(rbind, lapply(1:length(ab_names_rep), FUN=function(x){do_oe_ab_crispr(ab_names_rep[x], orgs_rep[x])}))
df_oe_org_3 <- do.call(rbind, lapply(1:length(ab_names_rep_3), FUN=function(x){do_oe_ab_crispr(ab_names_rep_3[x], orgs_rep_3[x])}))

df_oe_org_a <- df_oe_org[as.numeric(as.numeric(df_oe_org$PROP > 0.05) * as.numeric(df_oe_org$PROP < 0.95)) == 1,]
df_oe_org_2a <- df_oe_org_2[as.numeric(as.numeric(df_oe_org_2$PROP > 0.05) * as.numeric(df_oe_org_2$PROP < 0.95)) == 1,]
df_oe_org_3a <- df_oe_org_3[as.numeric(as.numeric(df_oe_org_3$PROP > 0.05) * as.numeric(df_oe_org_3$PROP < 0.95)) == 1,]

df_oe_org$LOG_OE_CAT <- factor(df_oe_org$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.01]", "(-0.01,0.01]", "(0.01,0.5]", "(0.5,1.0]", "1.0+")))
df_oe_org_2$LOG_OE_CAT <- factor(df_oe_org_2$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.2]", "(-0.2,0.2]", "(0.2,0.5]", "(0.5,1.0]", "1.0+")))
df_oe_org_3$LOG_OE_CAT <- factor(df_oe_org_3$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.2]", "(-0.2,0.2]", "(0.2,0.5]", "(0.5,1.0]", "1.0+")))

df_oe_org_a$LOG_OE_CAT <- factor(df_oe_org_a$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.01]", "(-0.01,0.01]", "(0.01,0.5]", "(0.5,1.0]", "1.0+")))
df_oe_org_2a$LOG_OE_CAT <- factor(df_oe_org_2a$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.2]", "(-0.2,0.2]", "(0.2,0.5]", "(0.5,1.0]", "1.0+")))
df_oe_org_3a$LOG_OE_CAT <- factor(df_oe_org_3a$LOG_OE_CAT, levels = rev(c("< -1.0", "(-1.0,-0.5]", "(-0.5,-0.2]", "(-0.2,0.2]", "(0.2,0.5]", "(0.5,1.0]", "1.0+")))

lst_cl <- diverge_hcl(7, c=100, l=c(50,90), power=1)

gg <- ggplot(data=df_oe_org_2a, aes(y=factor(ORG, levels = rev(unique(df_oe_org_2a$ORG))),x=factor(AB)))+geom_tile(aes(fill=LOG_OE_CAT,color=factor(SIG),width=0.95,height=0.95),size=0.4)
gg <- gg + labs(y="ORG", x="AB", fill="log(O/E)", color="Association")
gg <- gg + theme(axis.text.y=element_text(color="black",face="italic"),axis.text.x=element_text(angle=75,hjust=1.05,vjust=1.05,color="black"))
gg <- gg + xlab("")+ylab("")+ggtitle("")+guides(colour=FALSE)
gg <- gg + scale_color_manual(values=c("No Assoc."="white","Assoc."="black"))
gg <- gg + scale_fill_manual(values = rev(lst_cl))

gg1 <- ggplot(data=df_oe_org, aes(y=factor(ORG, levels = rev(unique(df_oe_org$ORG))),x=factor(AB)))+geom_tile(aes(fill=LOG_OE_CAT,color=factor(SIG),width=0.95,height=0.95),size=0.4)
gg1 <- gg1 + labs(y="ORG", x="AB", fill="log(O/E)", color="Association")
gg1 <- gg1 + theme(axis.text.y=element_text(color="black",face="italic"),axis.text.x=element_text(angle=75,hjust=1.05,vjust=1.05,color="black"))
gg1 <- gg1 + xlab("")+ylab("")+ggtitle("")+guides(colour=FALSE)
gg1 <- gg1 + scale_color_manual(values=c("No Assoc."="white","Assoc."="black"))
gg1 <- gg1 + scale_fill_manual(values = rev(lst_cl)[-1])

gg2 <- ggplot(data=df_oe_org_3a, aes(y=factor(ORG, levels = rev(unique(df_oe_org_3a$ORG))),x=factor(AB)))+geom_tile(aes(fill=LOG_OE_CAT,color=factor(SIG),width=0.95,height=0.95),size=0.4)
gg2 <- gg2 + labs(y="ORG", x="AB", fill="log(O/E)", color="Association")
gg2 <- gg2 + theme(axis.text.y=element_text(color="black",face="italic"),axis.text.x=element_text(angle=75,hjust=1.05,vjust=1.05,color="black"))
gg2 <- gg2 + xlab("")+ylab("")+ggtitle("")+guides(colour=FALSE)
gg2 <- gg2 + scale_color_manual(values=c("No Assoc."="white","Assoc."="black"))
gg2 <- gg2 + scale_fill_manual(values = rev(lst_cl))