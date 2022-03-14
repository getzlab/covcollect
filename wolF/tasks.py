from wolf import Task

class Covcollect(Task):
    inputs = {
      "bam" : None,
      "bai" : None,
      "intervals" : None, # path to BED file or value for even bins
      "subset_chr" : -1,
      "subset_start" : 0,
      "subset_end" : 0
    }
    script = """
    /app/covcollect -b ${bam} -i ${intervals} -o coverage.bed -c ${subset_chr} -s ${subset_start} -e ${subset_end}
    """
    output_patterns = { "coverage" : "coverage.bed" }
    docker = "gcr.io/broad-getzlab-workflows/covcollect:fraglen_v192"
