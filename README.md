[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](ngstools.png){width=20%} - NGS-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/ngs-tools/wiki)** for examples and workflows.

### **run_freebayes.sh**

The bash file **[run_freebayes.sh](run_freebayes.sh)** calls varisant with freebayes from mappings and a reference genome.
```bash
# Usage: run_freebayes.sh -i <bam file> -r <fasta reference>
# script version 1.0, 2016_09_28
# [optional: -m <minmapq|20>]
# [optional: -q <minbaseq|20>]
# [optional: -F <minaltfrac|0.01>]
# [optional: -C <minaltcnt|10>]
```

### **makeENSref.sh**

The bash file **[makeENSref.sh](makeENSref.sh)** creates a series of reference files from ENSembl FTP downloads.
```bash
# Usage: makeENSref.sh
# -o <organism (default to <homo_sapiens>)> 
# -b <build number (default to <GRCh38>)> 
# -r <release number (default to 88)>
# script version 1.0, 2017_04_05
# [-h for this help]
```

### **gepard_plot.sh**

The bash file **[gepard_plot.sh](gepard_plot.sh)** creates a xy-plot from two related assemblies. The Java GUI tools does the same but this script is applicable to multiple inputs in batch. The original tool can be found at **https://github.com/univieCUBE/gepard**.
```bash
# Usage: gepard_plot.sh -x <reference assembly> -y <draft assembly> -p <path to gepard.jar and matrices>
# script version 1.0, 2017_04_21
# [optional: -o <result folder>]
# [optional: -w <word size:10>]
# [optional: -W <window size:0>]
# [optional: -h <this help text>]
```

### **mauve_reorder.sh**

The bash file **[mauve_reorder.sh](mauve_reorder.sh)** reorders a draft-assembly based on a reference assembly at CLI. The Java GUI tools does the same but this script is applicable to multiple inputs in batch. The original tool can be found at **http://darlinglab.org/mauve/download.html**
```bash
# Usage: mauve_reorder.sh -i <draft assembly> -r <reference assembly> -p <mauve path>
# script version 1.0, 2017_04_21
# [optional: -o <result folder>]
# [optional: -m <available RAM|1G>]
# [optional: -h <this help text>]
```
<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
