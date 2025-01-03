version 1.0

### Define tasks
# task to run juicer pre
task juicer_pre {
	input {
		String fastqName
		String TopDir
		String GenomeAssembly
		String ChromSizes
		String? Enzyme
		Int? cpu = 64
		Int? mem = 150	
	}
	command {
		bash /opt/scripts/juicer.sh \
		-z ${GenomeAssembly} \
		-p ${ChromSizes} \
		~{"-s " + Enzyme} \
		-t ${cpu} \
		-d ${TopDir} \
		-D /opt
	}
	output {
		File hic = '${TopDir}/aligned/{fastqName}.hic'
	}
	runtime {
		cpu: "${cpu}"
		mem: "${mem}"
		docker: "faryabilab/juicer:v1"
	}
}

task hic2cool {
	input {
		File hic
		String TopDir
		String? resolutions = '0'
		Int? cpu = 12
		Int? mem = 64
		String fastqName
	}
	command {
		hic2cool convert \
		${hic} \
		-r ${resolutions} \
		-p ${cpu} \
		${TopDir}/${fastqName}.mcool
	}
	output {
		File mcool = '${TopDir}/${fastqName}.mcool'
	}
	runtime {
		cpu: "${cpu}"
		mem: "${mem}"
		docker: "faryabilab/hic2cool:v1"
	}
}

### Define workflow
workflow CooledJuicer {
	input {
		String TopDir
		String GenomeAssembly
		String ChromSizes
		String? Enzyme
		String resolutions=0
		File sampleList
	}	
	Array[Array[String]] samples = read_tsv(sampleList)
	scatter (sample in samples) {
		String sampleName = sample[0]
		call juicer_pre {
			input:
				fastqName=sampleName,
				TopDir=TopDir,
				GenomeAssembly=GenomeAssembly,
				ChromSizes=ChromSizes,
				Enzyme=Enzyme
		}
		call hic2cool {
			input:
				hic=juicer_pre.hic,
				TopDir=TopDir,
				resolutions=resolutions,
				fastqName=sampleName
		}
	}
	output {
		Array[File] hic = juicer_pre.hic
		Array[File] mcool = hic2cool.mcool
	}
}
