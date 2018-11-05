library(ggplot2)
library(grid)
library(colorspace)

pdata <- read.csv("anticrispr_crispr_cas_res_pm_genome_gc_size_v2.csv")
df_num_orgs <- read.csv("num_organisms.csv")
df_num_orgs <- df_num_orgs[df_num_orgs$Freq > 9,]

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

do_fet_ab_acr <- function(ab, org) {
	x <- as.numeric(pdata$ACR_COUNT > 0)
	y <- as.numeric(pdata[[ab]])
	
	df <- data.frame(X = x, Y = y)
	
	if (org != "All") {
		inds <- grep(org, as.character(pdata$ORGANISM))
		df <- df[inds,]
	}
	
	n11 <- sum(df$X * df$Y)
	n10 <- sum(df$X * (1 - df$Y))
	n01 <- sum((1 -df$X) * df$Y)
	n00 <- sum((1 -df$X) * (1 - df$Y))
	tot <- sum(c(n11, n10, n01, n00))
	
	tab <- matrix(c(n11, n10, n01, n00), nrow=2)
	fet <- fisher.test(tab)
	
	or <- fet$estimate
	lor <- log(or)
	pval <- fet$p.value
	
	assoc_type <- c("Negative", "None", "Positive")
	assoc <- assoc_type[2 + (sign(lor) * as.numeric(pval < 1e-5))]
	dfout <- data.frame(AB = c(ab), ORG = c(org), N_GENOMES = c(tot),
	                    N_ARG_ACR = c(n11), N_NARG_ACR = c(n10), N_ARG_NACR = c(n01), N_NARG_NACR = c(n00),
	                    OR = c(or), LOR = c(lor), PVAL = c(pval), ASSOC = c(assoc))
	print(c(as.character(ab),as.character(org)))
	return(dfout)
}

do_fet_ab_crispr <- function(ab, org) {
	x <- as.numeric(pdata$CAS_COUNT > 0)
	y <- as.numeric(pdata[[ab]])
	
	df <- data.frame(X = x, Y = y)
	
	if (org != "All") {
		inds <- grep(org, as.character(pdata$ORGANISM))
		df <- df[inds,]
	}
	
	n11 <- sum(df$X * df$Y)
	n10 <- sum(df$X * (1 - df$Y))
	n01 <- sum((1 -df$X) * df$Y)
	n00 <- sum((1 -df$X) * (1 - df$Y))
	tot <- sum(c(n11, n10, n01, n00))
	
	tab <- matrix(c(n11, n10, n01, n00), nrow=2)
	fet <- fisher.test(tab)
	
	or <- fet$estimate
	lor <- log(or)
	pval <- fet$p.value
	
	assoc_type <- c("Negative", "None", "Positive")
	assoc <- assoc_type[2 + (sign(lor) * as.numeric(pval < 1e-5))]
	dfout <- data.frame(AB = c(ab), ORG = c(org), N_GENOMES = c(tot),
	                    N_ARG_CR = c(n11), N_NARG_CR = c(n10), N_ARG_NCR = c(n01), N_NARG_NCR = c(n00),
	                    OR = c(or), LOR = c(lor), PVAL = c(pval), ASSOC = c(assoc))
    print(c(as.character(ab),as.character(org)))
	return(dfout)
}

orgs <- unique(as.character(df_num_orgs$orgs_short))
ab_names_rep <- rep(ab_names, length(orgs))
orgs_rep <- rep(orgs, each = length(ab_names))

df_oe_org_all <- do.call(rbind, lapply(1:length(ab_names_rep), FUN=function(x){do_fet_ab_crispr(ab_names_rep[x], orgs_rep[x])}))

gg <- ggplot(data=df_oe_org_all, aes(x=factor(ORG, levels = unique(as.character(df_oe_org_all$ORG))[order(unique(as.character(df_oe_org_all$ORG)))]),y=factor(AB)))+geom_tile(aes(fill = factor(ASSOC, levels = c("Negative", "None", "Positive"))))
gg <- gg + labs(y="", x="", fill="Association (p < 1e-5)") + ggtitle("Coexistence of CRISPR-Cas and antibiotic resistance gene")
gg <- gg + theme(axis.text.y=element_text(color="black"),axis.text.x=element_text(hjust=1.05,angle=90,vjust=0.05,color="black",face="italic",size=4))
gg <- gg + scale_fill_manual(values = c("blue", "grey80", "red"))