version 1.0

workflow CallDenovo {
	input {
		File family_vcf
		File family_vcf_index

		String family_id
		String proband_id
		String father_id
		String mother_id
		String proband_sex

		File onekg_vcf
		File onekg_vcf_index

		File reference = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
		File reference_index = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
		File reference_dict = "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"
	}

	call CreatePed {
		input:
			family_id = family_id,
			proband_id = proband_id,
			mother_id = mother_id,
			father_id = father_id,
			proband_sex = proband_sex
	}

	call CalculateGenotypePosteriors {
		input:
			family_vcf = family_vcf,
			family_vcf_index = family_vcf_index,
			family_ped = CreatePed.family_ped,
			family_id = family_id,
			reference = reference,
			reference_index = reference_index,
			reference_dict = reference_dict,
			onekg_vcf = onekg_vcf,
			onekg_vcf_index = onekg_vcf_index
	}

	call VariantFiltration {
		input:
			posterior_vcf = CalculateGenotypePosteriors.posterior_vcf,
			posterior_vcf_index = CalculateGenotypePosteriors.posterior_vcf_index,
			family_id = family_id,
			reference = reference,
			reference_index = reference_index,
			reference_dict = reference_dict
	}

	call VariantAnnotator {
		input:
			filtered_vcf = VariantFiltration.filtered_vcf,
			filtered_vcf_index = VariantFiltration.filtered_vcf_index,
			family_ped = CreatePed.family_ped,
			family_id = family_id,
			reference = reference,
			reference_index = reference_index,
			reference_dict = reference_dict
	}

	call snpEff {
		input:
			denovo_vcf = VariantAnnotator.denovo_vcf,
			denovo_vcf_index = VariantAnnotator.denovo_vcf_index,
			family_id = family_id
	}

	call snpSift {
		input:
			annotated_vcf = snpEff.annotated_vcf,
			family_id = family_id
	}

	output {
		File high_confidence_impact_denovos = snpSift.high_impact_vcf
		File snpEff_summary = snpEff.snpEff_summary
	}
}

task CreatePed {
	input {
		String family_id
		String proband_id
		String mother_id
		String father_id
		String proband_sex
	}

	command <<<
		if [ "~{proband_sex}" = "F" ]; then
			PROBAND_SEX=2
		elif [ "~{proband_sex}" = "M" ]; then
			PROBAND_SEX=1
		else
			PROBAND_SEX=other
		fi

		echo -e "~{family_id}\t~{proband_id}\t~{father_id}\t~{mother_id}\t$PROBAND_SEX\t2" > ~{family_id}.ped
		echo -e "~{family_id}\t~{father_id}\t0\t0\t1\t1" >> ~{family_id}.ped
		echo -e "~{family_id}\t~{mother_id}\t0\t0\t2\t1" >> ~{family_id}.ped
	>>>

	output {
		File family_ped = "~{family_id}.ped"
	}

	runtime {
		docker: "ubuntu:xenial"
		cpu: 1
		memory: "3.75 GB"
		disks: "local-disk 10 HDD"
	}
}

task CalculateGenotypePosteriors {
	input {
		File family_vcf
		File family_vcf_index
		File family_ped
		String family_id

		File reference
		File reference_index
		File reference_dict

		File onekg_vcf
		File onekg_vcf_index
	}

	Int disk_size = ceil((size(family_vcf, "GB") + size(reference, "GB") + size(onekg_vcf, "GB")) * 2 + 50)

	command <<<
		java -jar /opt/extras/gatk/GenomeAnalysisTK.jar \
			-T CalculateGenotypePosteriors \
			-R ~{reference} \
			-supporting ~{onekg_vcf} \
			-V ~{family_vcf} \
			-pedValidationType SILENT \
			-ped ~{family_ped} \
			-o ~{family_id}.withPosteriors.vcf
	>>>

	output {
		File posterior_vcf = "~{family_id}.withPosteriors.vcf"
		File posterior_vcf_index = "~{family_id}.withPosteriors.vcf.idx"
	}

	runtime {
		docker: "gcr.io/cool-benefit-817/gatk:3.8.0"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " LOCAL"
	}
}

task VariantFiltration {
	input {
		File posterior_vcf
		File posterior_vcf_index
		String family_id

		File reference
		File reference_index
		File reference_dict
	}

	Int disk_size = ceil((size(posterior_vcf, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		java -jar /opt/extras/gatk/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R ~{reference} \
			--variant ~{posterior_vcf} \
			--genotypeFilterExpression "GQ < 20" \
			--genotypeFilterName "GQ20" \
			-o ~{family_id}.filtered.vcf 
	>>>

	output {
		File filtered_vcf = "~{family_id}.filtered.vcf"
		File filtered_vcf_index = "~{family_id}.filtered.vcf.idx"
	}

	runtime {
		docker: "gcr.io/cool-benefit-817/gatk:3.8.0"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " LOCAL"
	}
}

task VariantAnnotator {
	input {
		File filtered_vcf
		File filtered_vcf_index
		File family_ped
		String family_id

		File reference
		File reference_index
		File reference_dict
	}

	Int disk_size = ceil((size(filtered_vcf, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		java -jar /opt/extras/gatk/GenomeAnalysisTK.jar \
		-T VariantAnnotator \
		-R ~{reference} \
		-V ~{filtered_vcf} \
		-pedValidationType SILENT \
		-ped ~{family_ped} \
		-A PossibleDeNovo \
		-o ~{family_id}.denovo.vcf
	>>>

	output {
		File denovo_vcf = "~{family_id}.denovo.vcf"
		File denovo_vcf_index = "~{family_id}.denovo.vcf.idx"
	}

	runtime {
		docker: "gcr.io/cool-benefit-817/gatk:3.8.0"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " LOCAL"
	}
}

task snpEff {
	input {
		File denovo_vcf
		File denovo_vcf_index
		String family_id
	}

	Int disk_size = ceil(size(denovo_vcf, "GB") * 4 + 50)

	command <<<
		. snpEff_wrapper.sh \
		-C $CONFIG \
		-v ~{denovo_vcf} \
		-d GRCh38.86 \
		-S Homo_sapiens \
		-m Vertebrate_Mitochondrial

		# first filter only denovos
		java -jar /opt/snpEff/SnpSift.jar filter \
		"(exists hiConfDeNovo)" ~{denovo_vcf} \
		> ~{family_id}.denovo_only.vcf

		# annotate only high confidence de novo variants
		java -Xmx9G -jar /opt/snpEff/snpEff.jar \
		-c $CONFIG \
		$GENOME \
		~{family_id}.denovo_only.vcf \
		$ANNO_STRING \
		> ~{family_id}.denovo.annotated.vcf
	>>>

	output {
		File annotated_vcf = "~{family_id}.denovo.annotated.vcf"
		File snpEff_summary = "snpEff_summary.html"
	}

	runtime {
		docker: "gcr.io/cool-benefit-817/snpeff:4.3t"
		cpu: 2
		memory: "13 GB"
		disks: "local-disk " + disk_size + " HDD"
	}	
}

task snpSift {
	input {
		File annotated_vcf
		String family_id
	}

	Int disk_size = ceil(size(annotated_vcf, "GB") * 2 + 20)

	command <<<
		java -jar /opt/snpEff/SnpSift.jar filter \
		"(ANN[*].IMPACT = 'HIGH')" \
		~{annotated_vcf} \
		> ~{family_id}.denovo.annotated.high_impact.vcf
	>>>

	output {
		File high_impact_vcf = "~{family_id}.denovo.annotated.high_impact.vcf" 
	}

	runtime {
		docker: "gcr.io/cool-benefit-817/snpeff:4.3t"
		cpu: 2
		memory: "13 GB"
		disks: "local-disk " + disk_size + " HDD"
	}	
}
