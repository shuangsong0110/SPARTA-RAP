
Gene_Centric_Coding_Results_Summary2 <- function(agds_dir,gene_centric_coding_jobs_num,input_path,output_path,gene_centric_results_name,
                                                obj_nullmodel=NULL,known_loci=NULL,cMAC_cutoff=0,SNV_cutoff=10,
                                                method_cond=c("optimal","naive"),rare_maf_cutoff=0.01,
                                                QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                                Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                                Use_annotation_weights=FALSE,Annotation_name=NULL,
                                                alpha=2.5E-06,manhattan_plot=FALSE,QQ_plot=FALSE,
                                                cond_null_model_name=NULL,cond_null_model_dir=NULL,
                                                SPA_p_filter=FALSE,p_filter_cutoff=0.05){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## SPA status
	# if(!is.null(obj_nullmodel$use_SPA))
	# {
	# 	use_SPA <- obj_nullmodel$use_SPA
	# }else
	# {
	# 	use_SPA <- FALSE
	# }
	use_SPA=F
	#######################################################
	#     summarize unconditional analysis results
	#######################################################

	results_coding_genome <- c()

	for(kk in 1:gene_centric_coding_jobs_num)
	{
		#print(kk)
	  if(file.exists(paste0(input_path,gene_centric_results_name,"_",kk,".Rdata"))){
		results_coding <- get(load(paste0(input_path,gene_centric_results_name,"_",kk,".Rdata")))

		results_coding_genome <- c(results_coding_genome,results_coding)
	  }else{
	    print(kk)
	  }
	}

	results_plof_genome <- c()
	results_plof_ds_genome <- c()
	results_missense_genome <- c()
	results_disruptive_missense_genome <- c()
	results_synonymous_genome <- c()
#  rr <- unlist(lapply(results_coding_genome, function(x){
# if(is.null(x[1])){return(NA)}else{return(x[1])}
# }))
	for(kk in 1:length(results_coding_genome))
	#for(kk in 15726:15800)
	{
		results <- results_coding_genome[[kk]]

		if(is.null(results)==FALSE)
		{
			### plof
			if(results[3]=="plof")
			{
				results_plof_genome <- rbind(results_plof_genome,results)
			}
			### plof_ds
			if(results[3]=="plof_ds")
			{
				results_plof_ds_genome <- rbind(results_plof_ds_genome,results)
			}
			### missense
			if(results[3]=="missense")
			{
				results_missense_genome <- rbind(results_missense_genome,results)
			}
			### disruptive_missense
			if(results[3]=="disruptive_missense")
			{
				results_disruptive_missense_genome <- rbind(results_disruptive_missense_genome,results)
			}
			### synonymous
			if(results[3]=="synonymous")
			{
				results_synonymous_genome <- rbind(results_synonymous_genome,results)
			}
		}

		if(kk%%1000==0)
		{
			print(kk)
		}
	}

	###### cMAC_cutoff
	# plof
	results_plof_genome <- results_plof_genome[results_plof_genome[,"cMAC"]>cMAC_cutoff,]
	# plof + disruptive missense
	results_plof_ds_genome <- results_plof_ds_genome[results_plof_ds_genome[,"cMAC"]>cMAC_cutoff,]
	# missense
	results_missense_genome <- results_missense_genome[results_missense_genome[,"cMAC"]>cMAC_cutoff,]
	# disruptive missense
	results_disruptive_missense_genome <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"cMAC"]>cMAC_cutoff,]
	# synonymous
	results_synonymous_genome <- results_synonymous_genome[results_synonymous_genome[,"cMAC"]>cMAC_cutoff,]

	####### SNV_cutoff
	# plof
	results_plof_genome <- results_plof_genome[results_plof_genome[,"#SNV"]>SNV_cutoff,]
	# plof + disruptive missense
	results_plof_ds_genome <- results_plof_ds_genome[results_plof_ds_genome[,"#SNV"]>SNV_cutoff,]
	# missense
	results_missense_genome <- results_missense_genome[results_missense_genome[,"#SNV"]>SNV_cutoff,]
	# disruptive missense
	results_disruptive_missense_genome <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"#SNV"]>SNV_cutoff,]
	# synonymous
	results_synonymous_genome <- results_synonymous_genome[results_synonymous_genome[,"#SNV"]>SNV_cutoff,]
	
	### recalculate missense pvalue
	if(SNV_cutoff > 0)
	{
		genes_name_disruptive_missense <- as.vector(unlist(results_disruptive_missense_genome[,1])) 
		genes_name_missense <- as.vector(unlist(results_missense_genome[,1])) 

		recal_id <- (1:dim(results_missense_genome)[1])[(!genes_name_missense%in%genes_name_disruptive_missense)]

		for(kk in recal_id)
		{
			print(kk)
			results_m <- results_missense_genome[kk,]
			
			if(use_SPA)
			{
				## disruptive missense cMAC < cut_off, set p-value to 1 
				results_missense_genome[kk,"Burden(1,25)-Disruptive"] <- 1
				results_missense_genome[kk,"Burden(1,1)-Disruptive"] <- 1

				apc_num <- (length(results_m)-10)/2
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1))
				
				## calculate STAAR-B
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
				
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B"] <- 1
					}
				}
				
				results_missense_genome[kk,"STAAR-B"] <- results_m["STAAR-B"]
				
				## calculate STAAR-B(1,25)
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num)])
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B(1,25)"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B(1,25)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B(1,25)"] <- 1
					}
				}
				
				results_missense_genome[kk,"STAAR-B(1,25)"] <- results_m["STAAR-B(1,25)"]
				
				## calculate STAAR-B(1,1)
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1))])
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B(1,1)"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B(1,1)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B(1,1)"] <- 1
					}
				}
				
				results_missense_genome[kk,"STAAR-B(1,1)"] <- results_m["STAAR-B(1,1)"]
			}else
			{
				## disruptive missense cMAC < cut_off, set p-value to 1 
			  results_m["Burden(1,25)-Disruptive"] <- 1
			  results_m["Burden(1,1)-Disruptive"] <- 1
			  results_m["SKAT(1,25)-Disruptive"] <- 1
			  results_m["SKAT(1,1)-Disruptive"] <- 1			
			  results_m["ACAT-V(1,25)-Disruptive"] <- 1
			  results_m["ACAT-V(1,1)-Disruptive"] <- 1			
				
				apc_num <- (length(results_m)-19)/6
				
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1))

				results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
				results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num)]))
				results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1))]))
				results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1))]))
				results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1))]))
				results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1))]))
				results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1))]))
			}
			results_missense_genome[kk,] <- results_m
		}

	}

	###### whole-genome results
	# plof
	save(results_plof_genome,file=paste0(output_path,"plof.Rdata"))
	# plof + disruptive missense
	save(results_plof_ds_genome,file=paste0(output_path,"plof_ds.Rdata"))
	# missense
	save(results_missense_genome,file=paste0(output_path,"missense.Rdata"))
	# disruptive missense
	save(results_disruptive_missense_genome,file=paste0(output_path,"disruptive_missense.Rdata"))
	# synonymous
	save(results_synonymous_genome,file=paste0(output_path,"synonymous.Rdata"))

	if(F)
	{
		###### significant results
		# plof
		plof_sig <- results_plof_genome[results_plof_genome[,"STAAR-B"]<alpha,,drop=FALSE]
		write.csv(plof_sig,file=paste0(output_path,"plof_sig.csv"))

		# missense
		missense_sig <- results_missense_genome[results_missense_genome[,"STAAR-B"]<alpha,,drop=FALSE]
		write.csv(missense_sig,file=paste0(output_path,"missense_sig.csv"))

		# synonymous
		synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"STAAR-B"]<alpha,,drop=FALSE]
		write.csv(synonymous_sig,file=paste0(output_path,"synonymous_sig.csv"))

		# plof_ds
		plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-B"]<alpha,,drop=FALSE]
		write.csv(plof_ds_sig,file=paste0(output_path,"plof_ds_sig.csv"))

		# disruptive_missense
		disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-B"]<alpha,,drop=FALSE]
		write.csv(disruptive_missense_sig,file=paste0(output_path,"disruptive_missense_sig.csv"))

		# coding results
		coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")],
		                    missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
		coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
		coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
		coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])

		write.csv(coding_sig,file=paste0(output_path,"coding_sig.csv"))
	}else
	{
		# plof
		plof_sig <- results_plof_genome[results_plof_genome[,"STAAR-O"]<alpha,,drop=FALSE]
		write.csv(plof_sig,file=paste0(output_path,"plof_sig.csv"))

		# missense
		missense_sig <- results_missense_genome[results_missense_genome[,"STAAR-O"]<alpha,,drop=FALSE]
		write.csv(missense_sig,file=paste0(output_path,"missense_sig.csv"))

		# synonymous
		synonymous_sig <- results_synonymous_genome[results_synonymous_genome[,"STAAR-O"]<alpha,,drop=FALSE]
		write.csv(synonymous_sig,file=paste0(output_path,"synonymous_sig.csv"))

		# plof_ds
		plof_ds_sig <- results_plof_ds_genome[results_plof_ds_genome[,"STAAR-O"]<alpha,,drop=FALSE]
		write.csv(plof_ds_sig,file=paste0(output_path,"plof_ds_sig.csv"))

		# disruptive_missense
		disruptive_missense_sig <- results_disruptive_missense_genome[results_disruptive_missense_genome[,"STAAR-O"]<alpha,,drop=FALSE]
		write.csv(disruptive_missense_sig,file=paste0(output_path,"disruptive_missense_sig.csv"))

		# coding results
		coding_sig <- rbind(plof_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
		                    missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
		coding_sig <- rbind(coding_sig,synonymous_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
		coding_sig <- rbind(coding_sig,plof_ds_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
		coding_sig <- rbind(coding_sig,disruptive_missense_sig[,c("Gene name","Chr","Category","#SNV","cMAC","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])

		write.csv(coding_sig,file=paste0(output_path,"coding_sig.csv"))
	}

	#######################################################
	#     conditional analysis
	#######################################################

	if(length(known_loci)!=0)
	{
		if(!use_SPA)
		{
			### plof
			plof_sig_cond <- c()
			if(length(plof_sig)!=0)
			{
				if((class(plof_sig)[1]!="matrix")&(class(plof_sig)[1]!="data.frame"))
				{
					plof_sig <- matrix(plof_sig,nrow=1)
				}

				for(k in 1:dim(plof_sig)[1])
				{
					chr <- as.numeric(plof_sig[k,2])
					gene_name <- as.character(plof_sig[k,1])
					category <- as.character(plof_sig[k,3])

					gds.path <- agds_dir[chr]
					genofile <- seqOpen(gds.path)

					res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
					                                     obj_nullmodel=obj_nullmodel,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
					                                     variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
					                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

					plof_sig_cond <- rbind(plof_sig_cond,res_cond)

					seqClose(genofile)
				}
			}

			write.csv(plof_sig_cond,file=paste0(output_path,"plof_sig_cond.csv"))

			### missense
			missense_sig_cond <- c()
			if(length(missense_sig)!=0)
			{
				if((class(missense_sig)[1]!="matrix")&(class(missense_sig)[1]!="data.frame"))
				{
					missense_sig <- matrix(missense_sig,nrow=1)
				}

				for(k in 1:dim(missense_sig)[1])
				{
					chr <- as.numeric(missense_sig[k,2])
					gene_name <- as.character(missense_sig[k,1])
					category <- as.character(missense_sig[k,3])

					gds.path <- agds_dir[chr]
					genofile <- seqOpen(gds.path)

					res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
					                                     obj_nullmodel=obj_nullmodel,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
					                                     variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
					                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

					missense_sig_cond <- rbind(missense_sig_cond,res_cond)

					seqClose(genofile)
				}
			}
			write.csv(missense_sig_cond,file=paste0(output_path,"missense_sig_cond.csv"))

			### synonymous
			synonymous_sig_cond <- c()
			if(length(synonymous_sig)!=0)
			{
				if((class(synonymous_sig)[1]!="matrix")&(class(synonymous_sig)[1]!="data.frame"))
				{
					synonymous_sig <- matrix(synonymous_sig,nrow=1)
				}

				for(k in 1:dim(synonymous_sig)[1])
				{
					chr <- as.numeric(synonymous_sig[k,2])
					gene_name <- as.character(synonymous_sig[k,1])
					category <- as.character(synonymous_sig[k,3])

					gds.path <- agds_dir[chr]
					genofile <- seqOpen(gds.path)

					res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
					                                     obj_nullmodel=obj_nullmodel,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
					                                     variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
					                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

					synonymous_sig_cond <- rbind(synonymous_sig_cond,res_cond)

					seqClose(genofile)
				}
			}
			write.csv(synonymous_sig_cond,file=paste0(output_path,"synonymous_sig_cond.csv"))

			### plof_ds
			plof_ds_sig_cond <- c()
			if(length(plof_ds_sig)!=0)
			{
				if((class(plof_ds_sig)[1]!="matrix")&(class(plof_ds_sig)[1]!="data.frame"))
				{
					plof_ds_sig <- matrix(plof_ds_sig,nrow=1)
				}

				for(k in 1:dim(plof_ds_sig)[1])
				{
					chr <- as.numeric(plof_ds_sig[k,2])
					gene_name <- as.character(plof_ds_sig[k,1])
					category <- as.character(plof_ds_sig[k,3])

					gds.path <- agds_dir[chr]
					genofile <- seqOpen(gds.path)

					res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
					                                     obj_nullmodel=obj_nullmodel,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
					                                     variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
					                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

					plof_ds_sig_cond <- rbind(plof_ds_sig_cond,res_cond)

					seqClose(genofile)
				}
			}
			write.csv(plof_ds_sig_cond,file=paste0(output_path,"plof_ds_sig_cond.csv"))

			### disruptive_missense
			disruptive_missense_sig_cond <- c()
			if(length(disruptive_missense_sig)!=0)
			{
				if((class(disruptive_missense_sig)[1]!="matrix")&(class(disruptive_missense_sig)[1]!="data.frame"))
				{
					disruptive_missense_sig <- matrix(disruptive_missense_sig,nrow=1)
				}

				for(k in 1:dim(disruptive_missense_sig)[1])
				{
					chr <- as.numeric(disruptive_missense_sig[k,2])
					gene_name <- as.character(disruptive_missense_sig[k,1])
					category <- as.character(disruptive_missense_sig[k,3])

					gds.path <- agds_dir[chr]
					genofile <- seqOpen(gds.path)

					res_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
					                                     obj_nullmodel=obj_nullmodel,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
					                                     variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
					                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
					                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

					disruptive_missense_sig_cond <- rbind(disruptive_missense_sig_cond,res_cond)

					seqClose(genofile)
				}
			}
			write.csv(disruptive_missense_sig_cond,file=paste0(output_path,"disruptive_missense_sig_cond.csv"))

			## coding cond
			coding_sig_cond <- rbind(plof_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")],
			                         missense_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
			coding_sig_cond <- rbind(coding_sig_cond,synonymous_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
			coding_sig_cond <- rbind(coding_sig_cond,plof_ds_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])
			coding_sig_cond <- rbind(coding_sig_cond,disruptive_missense_sig_cond[,c("Gene name","Chr","Category","#SNV","SKAT(1,25)","Burden(1,1)","ACAT-V(1,25)","STAAR-O")])

			write.csv(coding_sig_cond,file=paste0(output_path,"coding_sig_cond.csv"))
		}else
		{
			### synonymous
			synonymous_sig_cond <- c()
			if(length(synonymous_sig)!=0)
			{
				if((class(synonymous_sig)[1]!="matrix")&(class(synonymous_sig)[1]!="data.frame"))
				{
					synonymous_sig <- matrix(synonymous_sig,nrow=1)
				}

				for(k in 1:dim(synonymous_sig)[1])
				{
					chr <- as.numeric(synonymous_sig[k,2])
					gene_name <- as.character(synonymous_sig[k,1])
					category <- as.character(synonymous_sig[k,3])

					if(file.exists(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))
					{
						gds.path <- agds_dir[chr]
						genofile <- seqOpen(gds.path)

						obj_nullmodel_cond <- get(load(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))

						res_cond <- Gene_Centric_Coding_cond_spa(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
						                                         obj_nullmodel=obj_nullmodel_cond,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
						                                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
						                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
						                                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff)
						seqClose(genofile)
					}else
					{
						res_cond <- synonymous_sig[k,,drop=FALSE]
						res_cond[1,3] <- "synonymous_cond"
					}

					synonymous_sig_cond <- rbind(synonymous_sig_cond,res_cond)
				}
			}
			write.csv(synonymous_sig_cond,file=paste0(output_path,"synonymous_sig_cond.csv"))

			### plof
			plof_sig_cond <- c()
			if(length(plof_sig)!=0)
			{
				if((class(plof_sig)[1]!="matrix")&(class(plof_sig)[1]!="data.frame"))
				{
					plof_sig <- matrix(plof_sig,nrow=1)
				}

				for(k in 1:dim(plof_sig)[1])
				{
					chr <- as.numeric(plof_sig[k,2])
					gene_name <- as.character(plof_sig[k,1])
					category <- as.character(plof_sig[k,3])

					if(file.exists(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))
					{
						gds.path <- agds_dir[chr]
						genofile <- seqOpen(gds.path)

						obj_nullmodel_cond <- get(load(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))

						res_cond <- Gene_Centric_Coding_cond_spa(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
						                                         obj_nullmodel=obj_nullmodel_cond,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
						                                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
						                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
						                                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff)
						seqClose(genofile)
					}else
					{
						res_cond <- plof_sig[k,,drop=FALSE]
						res_cond[1,3] <- "plof_cond"
					}

					plof_sig_cond <- rbind(plof_sig_cond,res_cond)
				}
			}
			write.csv(plof_sig_cond,file=paste0(output_path,"plof_sig_cond.csv"))

			### plof_ds
			plof_ds_sig_cond <- c()
			if(length(plof_ds_sig)!=0)
			{
				if((class(plof_ds_sig)[1]!="matrix")&(class(plof_ds_sig)[1]!="data.frame"))
				{
					plof_ds_sig <- matrix(plof_ds_sig,nrow=1)
				}

				for(k in 1:dim(plof_ds_sig)[1])
				{
					chr <- as.numeric(plof_ds_sig[k,2])
					gene_name <- as.character(plof_ds_sig[k,1])
					category <- as.character(plof_ds_sig[k,3])

					if(file.exists(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))
					{
						gds.path <- agds_dir[chr]
						genofile <- seqOpen(gds.path)

						obj_nullmodel_cond <- get(load(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))

						res_cond <- Gene_Centric_Coding_cond_spa(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
						                                         obj_nullmodel=obj_nullmodel_cond,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
						                                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
						                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
						                                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff)
						seqClose(genofile)
					}else
					{
						res_cond <- plof_ds_sig[k,,drop=FALSE]
						res_cond[1,3] <- "plof_ds_cond"
					}

					plof_ds_sig_cond <- rbind(plof_ds_sig_cond,res_cond)
				}
			}
			write.csv(plof_ds_sig_cond,file=paste0(output_path,"plof_ds_sig_cond.csv"))

			### missense
			missense_sig_cond <- c()
			if(length(missense_sig)!=0)
			{
				if((class(missense_sig)[1]!="matrix")&(class(missense_sig)[1]!="data.frame"))
				{
					missense_sig <- matrix(missense_sig,nrow=1)
				}

				for(k in 1:dim(missense_sig)[1])
				{
					chr <- as.numeric(missense_sig[k,2])
					gene_name <- as.character(missense_sig[k,1])
					category <- as.character(missense_sig[k,3])

					if(file.exists(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))
					{
						gds.path <- agds_dir[chr]
						genofile <- seqOpen(gds.path)

						obj_nullmodel_cond <- get(load(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))

						res_cond <- Gene_Centric_Coding_cond_spa(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
						                                         obj_nullmodel=obj_nullmodel_cond,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
						                                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
						                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
						                                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff)
						seqClose(genofile)
					}else
					{
						res_cond <- missense_sig[k,,drop=FALSE]
						res_cond[1,3] <- "missense_cond"
					}

					missense_sig_cond <- rbind(missense_sig_cond,res_cond)
				}
			}
			write.csv(missense_sig_cond,file=paste0(output_path,"missense_sig_cond.csv"))

			### disruptive_missense
			disruptive_missense_sig_cond <- c()
			if(length(disruptive_missense_sig)!=0)
			{
				if((class(disruptive_missense_sig)[1]!="matrix")&(class(disruptive_missense_sig)[1]!="data.frame"))
				{
					disruptive_missense_sig <- matrix(disruptive_missense_sig,nrow=1)
				}

				for(k in 1:dim(disruptive_missense_sig)[1])
				{
					chr <- as.numeric(disruptive_missense_sig[k,2])
					gene_name <- as.character(disruptive_missense_sig[k,1])
					category <- as.character(disruptive_missense_sig[k,3])

					if(file.exists(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))
					{
						gds.path <- agds_dir[chr]
						genofile <- seqOpen(gds.path)

						obj_nullmodel_cond <- get(load(paste0(cond_null_model_dir,cond_null_model_name,".chr",chr,".Rdata")))

						res_cond <- Gene_Centric_Coding_cond_spa(chr=chr,gene_name=gene_name,category=category,genofile=genofile,
						                                         obj_nullmodel=obj_nullmodel_cond,known_loci=known_loci,rare_maf_cutoff=rare_maf_cutoff,
						                                         variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,QC_label=QC_label,
						                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
						                                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff)
						seqClose(genofile)
					}else
					{
						res_cond <- disruptive_missense_sig[k,,drop=FALSE]
						res_cond[1,3] <- "disruptive_missense_cond"
					}

					disruptive_missense_sig_cond <- rbind(disruptive_missense_sig_cond,res_cond)
				}
			}
			write.csv(disruptive_missense_sig_cond,file=paste0(output_path,"disruptive_missense_sig_cond.csv"))

			# coding results
			coding_sig_cond <- rbind(plof_sig_cond[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")],
			                         missense_sig_cond[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
			coding_sig_cond <- rbind(coding_sig_cond,synonymous_sig_cond[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
			coding_sig_cond <- rbind(coding_sig_cond,plof_ds_sig_cond[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])
			coding_sig_cond <- rbind(coding_sig_cond,disruptive_missense_sig_cond[,c("Gene name","Chr","Category","#SNV","cMAC","Burden(1,1)","STAAR-B")])

			write.csv(coding_sig_cond,file=paste0(output_path,"coding_sig_cond.csv"))
		}
	}

	## manhattan plot
	if(manhattan_plot)
	{
		### plof
		results_STAAR <- results_plof_genome[,c(1,2,dim(results_plof_genome)[2])]

		results_m <- c()
		for(i in 1:dim(results_STAAR)[2])
		{
			results_m <- cbind(results_m,unlist(results_STAAR[,i]))
		}

		colnames(results_m) <- colnames(results_STAAR)
		results_m <- data.frame(results_m,stringsAsFactors = FALSE)
		results_m[,2] <- as.numeric(results_m[,2])
		results_m[,3] <- as.numeric(results_m[,3])

		genes_info_manhattan <- dplyr::left_join(genes_info,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))

		genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
		colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof"

		### plof_ds
		results_STAAR <- results_plof_ds_genome[,c(1,2,dim(results_plof_ds_genome)[2])]

		results_m <- c()
		for(i in 1:dim(results_STAAR)[2])
		{
			results_m <- cbind(results_m,unlist(results_STAAR[,i]))
		}

		colnames(results_m) <- colnames(results_STAAR)
		results_m <- data.frame(results_m,stringsAsFactors = FALSE)
		results_m[,2] <- as.numeric(results_m[,2])
		results_m[,3] <- as.numeric(results_m[,3])

		genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
		genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
		colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"

		### missense
		if(!use_SPA)
		{
			results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-6)]
		}else
		{
			results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-2)]
		}

		results_m <- c()
		for(i in 1:dim(results_STAAR)[2])
		{
			results_m <- cbind(results_m,unlist(results_STAAR[,i]))
		}

		colnames(results_m) <- colnames(results_STAAR)
		results_m <- data.frame(results_m,stringsAsFactors = FALSE)
		results_m[,2] <- as.numeric(results_m[,2])
		results_m[,3] <- as.numeric(results_m[,3])

		genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
		genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
		colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"

		### disruptive_missense
		results_STAAR <- results_disruptive_missense_genome[,c(1,2,dim(results_disruptive_missense_genome)[2])]

		results_m <- c()
		for(i in 1:dim(results_STAAR)[2])
		{
			results_m <- cbind(results_m,unlist(results_STAAR[,i]))
		}

		colnames(results_m) <- colnames(results_STAAR)
		results_m <- data.frame(results_m,stringsAsFactors = FALSE)
		results_m[,2] <- as.numeric(results_m[,2])
		results_m[,3] <- as.numeric(results_m[,3])

		genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
		genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
		colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"

		### synonymous
		results_STAAR <- results_synonymous_genome[,c(1,2,dim(results_synonymous_genome)[2])]

		results_m <- c()
		for(i in 1:dim(results_STAAR)[2])
		{
			results_m <- cbind(results_m,unlist(results_STAAR[,i]))
		}

		colnames(results_m) <- colnames(results_STAAR)
		results_m <- data.frame(results_m,stringsAsFactors = FALSE)
		results_m[,2] <- as.numeric(results_m[,2])
		results_m[,3] <- as.numeric(results_m[,3])

		genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
		genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
		colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"

		## ylim
		genes_info_manhattan[genes_info_manhattan==0] <- 1
		coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-4):dim(genes_info_manhattan)[2]])
		min_y <- ceiling(-log10(coding_minp)) + 1

		pch <- c(0,1,2,3,4)
		figure1 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof,sig.level=alpha,pch=pch[1],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

		figure2 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$plof_ds,sig.level=alpha,pch=pch[2],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

		figure3 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$missense,sig.level=alpha,pch=pch[3],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

		figure4 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$disruptive_missense,sig.level=alpha,pch=pch[4],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

		figure5 <- manhattan_plot(genes_info_manhattan[,2], (genes_info_manhattan[,3]+genes_info_manhattan[,4])/2, genes_info_manhattan$synonymous,sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                          auto.key=T,key=list(space="top", columns=5, title="Functional Category", cex.title=1, points=TRUE,pch=pch,text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))

		
		if(T){
		  genes_info_manhattan$minp <- pmin(genes_info_manhattan$plof,genes_info_manhattan$plof_ds,
		                                    genes_info_manhattan$missense,genes_info_manhattan$disruptive_missense,
		                                    genes_info_manhattan$synonymous)
		  
		  plot_info_sig <- genes_info_manhattan
		  plot_info_sig$label <- ''
		  plot_info_sig$label[which(genes_info_manhattan$minp<alpha)] <-plot_info_sig[which(genes_info_manhattan$minp<alpha),1]
		  
		  plot_info_sig <- plot_info_sig[order(plot_info_sig$minp, decreasing = F), ]
		}
		figure_test <- manhattan_plot2(plot_info_sig[,2], (plot_info_sig[,3]+plot_info_sig[,4])/2, 
		                               plot_info_sig$minp,
		                               sig.level=alpha,pch=pch[5],col = c("blue4", "orange3"),ylim=c(0,min_y),
		                               label=plot_info_sig$label,
		                               auto.key=T,key=list(space="top", columns=5, title="Functional Category", 
		                                                   cex.title=1, points=TRUE,pch=pch,
		                                                   text=list(c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"))))
		
		print("Manhattan plot")

		png(paste0(output_path,"gene_centric_coding_manhattan.png"), width = 9, height = 6, units = 'in', res = 600)

		print(figure1)
		print(figure_test,newpage = FALSE)
		print(figure2,newpage = FALSE)
		print(figure3,newpage = FALSE)
		print(figure4,newpage = FALSE)
		print(figure5,newpage = FALSE)

		dev.off()
	}

	## Q-Q plot
	if(QQ_plot)
	{
		if(!manhattan_plot)
		{
			### plof
			results_STAAR <- results_plof_genome[,c(1,2,dim(results_plof_genome)[2])]

			results_m <- c()
			for(i in 1:dim(results_STAAR)[2])
			{
				results_m <- cbind(results_m,unlist(results_STAAR[,i]))
			}

			colnames(results_m) <- colnames(results_STAAR)
			results_m <- data.frame(results_m,stringsAsFactors = FALSE)
			results_m[,2] <- as.numeric(results_m[,2])
			results_m[,3] <- as.numeric(results_m[,3])

			genes_info_manhattan <- dplyr::left_join(genes_info,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))

			genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
			colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof"

			### plof_ds
			results_STAAR <- results_plof_ds_genome[,c(1,2,dim(results_plof_ds_genome)[2])]

			results_m <- c()
			for(i in 1:dim(results_STAAR)[2])
			{
				results_m <- cbind(results_m,unlist(results_STAAR[,i]))
			}

			colnames(results_m) <- colnames(results_STAAR)
			results_m <- data.frame(results_m,stringsAsFactors = FALSE)
			results_m[,2] <- as.numeric(results_m[,2])
			results_m[,3] <- as.numeric(results_m[,3])

			genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
			genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
			colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "plof_ds"

			### missense
			if(!use_SPA)
			{
				results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-6)]
			}else
			{
				results_STAAR <- results_missense_genome[,c(1,2,dim(results_missense_genome)[2]-2)]
			}

			results_m <- c()
			for(i in 1:dim(results_STAAR)[2])
			{
				results_m <- cbind(results_m,unlist(results_STAAR[,i]))
			}

			colnames(results_m) <- colnames(results_STAAR)
			results_m <- data.frame(results_m,stringsAsFactors = FALSE)
			results_m[,2] <- as.numeric(results_m[,2])
			results_m[,3] <- as.numeric(results_m[,3])

			genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
			genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
			colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "missense"

			### disruptive_missense
			results_STAAR <- results_disruptive_missense_genome[,c(1,2,dim(results_disruptive_missense_genome)[2])]

			results_m <- c()
			for(i in 1:dim(results_STAAR)[2])
			{
				results_m <- cbind(results_m,unlist(results_STAAR[,i]))
			}

			colnames(results_m) <- colnames(results_STAAR)
			results_m <- data.frame(results_m,stringsAsFactors = FALSE)
			results_m[,2] <- as.numeric(results_m[,2])
			results_m[,3] <- as.numeric(results_m[,3])

			genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
			genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
			colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "disruptive_missense"

			### synonymous
			results_STAAR <- results_synonymous_genome[,c(1,2,dim(results_synonymous_genome)[2])]

			results_m <- c()
			for(i in 1:dim(results_STAAR)[2])
			{
				results_m <- cbind(results_m,unlist(results_STAAR[,i]))
			}

			colnames(results_m) <- colnames(results_STAAR)
			results_m <- data.frame(results_m,stringsAsFactors = FALSE)
			results_m[,2] <- as.numeric(results_m[,2])
			results_m[,3] <- as.numeric(results_m[,3])

			genes_info_manhattan <- dplyr::left_join(genes_info_manhattan,results_m,by=c("chromosome_name"="Chr","hgnc_symbol"="Gene.name"))
			genes_info_manhattan[is.na(genes_info_manhattan)] <- 1
			colnames(genes_info_manhattan)[dim(genes_info_manhattan)[2]] <- "synonymous"

			## ylim
			coding_minp <- min(genes_info_manhattan[,(dim(genes_info_manhattan)[2]-4):dim(genes_info_manhattan)[2]])
			min_y <- ceiling(-log10(coding_minp)) + 1
		}

		print("Q-Q plot")
		cex_point <- 1

		png(paste0(output_path,"gene_centric_coding_qqplot.png"), width = 8, height = 8, units = 'in', res = 600)

		### plof
		## remove unconverged p-values
		observed <- sort(genes_info_manhattan$plof[genes_info_manhattan$plof < 1])
		lobs <- -(log10(observed))

		expected <- c(1:length(observed))
		lexp <- -(log10(expected / (length(expected)+1)))

		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=0, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
		     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		### plof_ds
		## remove unconverged p-values
		observed <- sort(genes_info_manhattan$plof_ds[genes_info_manhattan$plof_ds < 1])

		lobs <- -(log10(observed))

		expected <- c(1:length(observed))
		lexp <- -(log10(expected / (length(expected)+1)))

		par(new=T)
		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=1, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
		     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		### missense
		## remove unconverged p-values
		observed <- sort(genes_info_manhattan$missense[genes_info_manhattan$missense < 1])

		lobs <- -(log10(observed))

		expected <- c(1:length(observed))
		lexp <- -(log10(expected / (length(expected)+1)))

		par(new=T)
		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=2, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
		     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		### disruptive_missense
		## remove unconverged p-values
		observed <- sort(genes_info_manhattan$disruptive_missense[genes_info_manhattan$disruptive_missense < 1])


		lobs <- -(log10(observed))

		expected <- c(1:length(observed))
		lexp <- -(log10(expected / (length(expected)+1)))

		par(new=T)
		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=3, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
		     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		### synonymous
		## remove unconverged p-values
		observed <- sort(genes_info_manhattan$synonymous[genes_info_manhattan$synonymous < 1])

		lobs <- -(log10(observed))

		expected <- c(1:length(observed))
		lexp <- -(log10(expected / (length(expected)+1)))

		par(new=T)
		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=4, cex=cex_point, xlim = c(0, 5), ylim = c(0, min_y),
		     xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		     font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		legend("topleft",legend=c("pLoF","pLoF+D","Missense","Disruptive Missense","Synonymous"),ncol=1,bty="o",box.lwd=1,pch=0:4,cex=1.5,text.font=2)

		dev.off()
	}

}




library(lattice)
manhattan_plot2<-function(chr, pos, pvalue,
                          sig.level=NA, annotate=NULL, label=NULL, ann.default=list(),
                          should.thin=T, thin.pos.places=2, thin.logp.places=2,
                          xlab="Chromosome", ylab=expression(-log[10](p-value)),
                          col=c("gray","darkgray"), panel.extra=NULL, pch=20,
                          use_logp=FALSE, cex=0.8, ...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  # Make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[, drop=T]
  }
  
  # Make sure positions are in kbp
  if (any(pos > 1e6)) pos <- pos / 1e6
  
  # Calculate absolute genomic position from relative chromosomal positions
  posmin <- tapply(pos, chr, min)
  posmax <- tapply(pos, chr, max)
  posshift <- head(c(0, cumsum(posmax)), -1)
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr]
  getGenPos <- function(cchr, cpos) {
    p <- posshift[as.character(cchr)] + cpos
    return(p)
  }
  
  # Parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default <- list(x="peak", y="peak", adj=NULL, pos=3, offset=0.5,
                        col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label <- function(rawval, groupname) {
    r <- list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval) >= 1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    }
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings <- vector("list", length(levels(grp)))
  ann.settings[[1]] <- list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings) > 1) {
    lcols <- trellis.par.get("superpose.symbol")$col
    lfills <- trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]] <- list(pch=pch,
                                col=lcols[(i-2) %% length(lcols) +1],
                                fill=lfills[(i-2) %% length(lfills) +1],
                                cex=cex, label=label.default)
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings) <- levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i > 1) {
      ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)
    }
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label,
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  
  # If label is provided, set it for annotation
  if (!is.null(label)) {
    if (length(label) != length(pvalue)) stop("Length of label must match pvalue")
  }
  
  # Reduce number of points plotted
  if (should.thin) {
    if (!use_logp) {
      thinned <- unique(data.frame(
        logp = round(-log10(pvalue), thin.logp.places),
        pos = round(genpos, thin.pos.places),
        chr = chr,
        grp = grp
      ))
    } else {
      thinned <- unique(data.frame(
        logp = round(pvalue, thin.logp.places),
        pos = round(genpos, thin.pos.places),
        chr = chr,
        grp = grp
      ))
    }
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    if (!use_logp) {
      logp <- -log10(pvalue)
    } else {
      logp <- pvalue
    }
  }
  rm(pos, pvalue)
  gc()
  
  # Custom axis to print chromosome names
  axis.chr <- function(side, ...) {
    if (side == "bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax + posmin) / 2 + posshift),
                 labels=levels(chr),
                 ticks=F, rot=0,
                 check.overlap=F)
    } else if (side == "top" || side == "right") {
      panel.axis(side=side, draw.labels=F, ticks=F)
    } else {
      axis.default(side=side, ...)
    }
  }
  
  # Make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr <- function(x, y, ...) {
    A <- list()
    maxy <- ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0))) + 0.5
    A$ylim <- c(0, maxy)
    A
  }
  
  xyplot(logp ~ genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings,
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if (!is.na(sig.level)) {
             # Add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=3, col="red")
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos)
           if (!is.null(panel.extra)) {
             panel.extra(x, y, getgenpos, ...)
           }
         },
         panel.groups=function(x, y, ..., subscripts, group.number) {
           A <- list(...)
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           # do.call("panel.xyplot", A)
           
           # Draw labels if provided
           print(table(label))
           # print(gs$label$show)
           if (!is.null(label)) {
             #   print(label[subscripts])
             panel.text(x, y, labels=label[subscripts], pos=3, offset=.9, cex=.85, col="blue")
           }
         },
         xlab=xlab, ylab=ylab,
         panel.extra=panel.extra, getgenpos=getGenPos, ...)
}

